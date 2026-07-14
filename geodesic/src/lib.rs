//! Run orchestration for the `geodesic` binary (SAD.md §9.2): the pieces the
//! CLI drives but that also need to be callable from integration tests
//! (golden reference §13.7, determinism §13.5) without shelling out. The
//! binary crate is therefore lib + bin; `main.rs` only parses arguments and
//! calls in here. Kept in the binary crate, not `geodesic-engine`, because
//! the full loop needs `geodesic-io` (DCD + CSV writers) and the crate graph
//! (SAD.md §9.3) forbids the engine from depending on I/O.

use std::path::{Path, PathBuf};

use geodesic_core::{Axis, ComputeBackend, IoError, NumericalError, SimError, SimParams, SimState};
use geodesic_engine::constraint::{constrain_velocities, promote_hydrogen_bonds};
use geodesic_engine::cpu_backend::CpuBackend;
use geodesic_engine::force::{bonded, nonbonded};
use geodesic_engine::integrator::half_kick;
use geodesic_engine::neighbor;
use geodesic_io::config::{Backend, Config, DriftAction};
use geodesic_io::dcd::DcdWriter;
use geodesic_io::export::EnergyLogWriter;
use geodesic_io::{inpcrd, prmtop};

/// 1 amu·(Å/ps)² in kcal/mol is 1/20.455² — the same AKMA time constant the
/// DCD and inpcrd velocity units already use (geodesic-io), kept byte-identical
/// so kinetic energy and the file-format velocity conversions never disagree.
const KINETIC_ANG_PS_TO_KCAL: f64 = 1.0 / (20.455 * 20.455);

/// Defaults for the `energy` subcommand, which has only prmtop + inpcrd and no
/// config: standard AMBER cutoffs and a box large enough that no pair is ever
/// wrapped (the M1 fixtures are non-periodic GBSA systems). Printed alongside
/// the result so the number is reproducible.
const ENERGY_DEFAULT_R_SWITCH: f64 = 10.0;
const ENERGY_DEFAULT_R_CUTOFF: f64 = 12.0;
const ENERGY_DEFAULT_R_SKIN: f64 = 14.0;
const ENERGY_DEFAULT_BOX: f64 = 1000.0;

/// Per-term potential energy breakdown from the `energy` subcommand.
pub struct EnergyReport {
    pub bond: f64,
    pub angle: f64,
    pub dihedral: f64,
    pub nonbonded: f64,
    pub total: f64,
    pub r_switch: f64,
    pub r_cutoff: f64,
}

/// What a completed `run` produced, for the CLI to report and for tests to
/// assert on.
pub struct RunSummary {
    pub n_steps: u64,
    pub n_frames: usize,
    pub final_potential: f64,
    pub final_kinetic: f64,
    pub trajectory: PathBuf,
    pub energy_log: PathBuf,
}

fn read_file(path: &Path) -> Result<String, SimError> {
    std::fs::read_to_string(path).map_err(|e| SimError::Io(IoError { path: path.to_path_buf(), source: e }))
}

/// Resolves an input/output path from config against the config file's own
/// directory; an absolute path passes through unchanged.
fn resolve(base: &Path, path: &Path) -> PathBuf {
    base.join(path)
}

fn kinetic_energy(state: &SimState, mass: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (i, &m) in mass.iter().enumerate() {
        let v2 = state.vel_x[i] * state.vel_x[i]
            + state.vel_y[i] * state.vel_y[i]
            + state.vel_z[i] * state.vel_z[i];
        sum += m * v2;
    }
    0.5 * sum * KINETIC_ANG_PS_TO_KCAL
}

fn check_forces_finite(state: &SimState, step: u64) -> Result<(), NumericalError> {
    for i in 0..state.force_x.len() {
        if !state.force_x[i].is_finite() {
            return Err(NumericalError::NanInForce { step, atom: i, component: Axis::X });
        }
        if !state.force_y[i].is_finite() {
            return Err(NumericalError::NanInForce { step, atom: i, component: Axis::Y });
        }
        if !state.force_z[i].is_finite() {
            return Err(NumericalError::NanInForce { step, atom: i, component: Axis::Z });
        }
    }
    Ok(())
}

fn check_positions_finite(state: &SimState, step: u64) -> Result<(), NumericalError> {
    for i in 0..state.pos_x.len() {
        if !state.pos_x[i].is_finite() {
            return Err(NumericalError::NanInPos { step, atom: i, component: Axis::X });
        }
        if !state.pos_y[i].is_finite() {
            return Err(NumericalError::NanInPos { step, atom: i, component: Axis::Y });
        }
        if !state.pos_z[i].is_finite() {
            return Err(NumericalError::NanInPos { step, atom: i, component: Axis::Z });
        }
    }
    Ok(())
}

/// Evaluates the full force-field potential (bond + angle + dihedral + LJ) at
/// the coordinates in the inpcrd, with no constraint promotion, for the
/// `geodesic energy` convenience subcommand (SAD.md §9.2, §12.3). This is the
/// number a user copies into `config.toml`'s `total_energy` before a run.
pub fn energy_from_files(prmtop_path: &Path, inpcrd_path: &Path) -> Result<EnergyReport, SimError> {
    let prmtop_text = read_file(prmtop_path)?;
    let (atoms, topology) = prmtop::parse(&prmtop_text)?;
    let n_atoms = atoms.mass.len();
    let inpcrd_text = read_file(inpcrd_path)?;
    let state = inpcrd::parse(&inpcrd_text, n_atoms, false)?;

    let mut fx = vec![0.0; n_atoms];
    let mut fy = vec![0.0; n_atoms];
    let mut fz = vec![0.0; n_atoms];

    let bond = bonded::compute_bond_forces(&state, &topology, &mut fx, &mut fy, &mut fz);
    let angle = bonded::compute_angle_forces(&state, &topology, &mut fx, &mut fy, &mut fz);
    let dihedral = bonded::compute_dihedral_forces(&state, &topology, &mut fx, &mut fy, &mut fz);

    let params = SimParams {
        n_atoms,
        n_steps: 0,
        dt: 0.0,
        box_size: [ENERGY_DEFAULT_BOX; 3],
        periodic: false,
        r_cutoff: ENERGY_DEFAULT_R_CUTOFF,
        r_skin: ENERGY_DEFAULT_R_SKIN,
        r_switch: ENERGY_DEFAULT_R_SWITCH,
        max_constr_iter: 100,
        constr_tol: 1e-10,
        frame_interval: 1,
        n_threads: 1,
        total_energy: 0.0,
    };
    let mut list_state = SimState::new(n_atoms);
    list_state.pos_x = state.pos_x.clone();
    list_state.pos_y = state.pos_y.clone();
    list_state.pos_z = state.pos_z.clone();
    let list = neighbor::build(&mut list_state, &params, &topology);
    let nonbonded = nonbonded::compute_pair_forces(
        &state,
        &atoms,
        &list.pair_i,
        &list.pair_j,
        list.r_cutoff,
        list.r_switch,
        params.box_size,
        &mut fx,
        &mut fy,
        &mut fz,
    );

    let total = bond + angle + dihedral + nonbonded;
    Ok(EnergyReport {
        bond,
        angle,
        dihedral,
        nonbonded,
        total,
        r_switch: ENERGY_DEFAULT_R_SWITCH,
        r_cutoff: ENERGY_DEFAULT_R_CUTOFF,
    })
}

/// Parses `config.toml`, resolving all file paths against the config file's
/// directory, then runs the simulation to completion (SAD.md §9.2).
pub fn run_from_config_file(config_path: &Path) -> Result<RunSummary, SimError> {
    let config_text = read_file(config_path)?;
    let config = Config::from_toml_str(&config_text)?;
    let base = config_path.parent().unwrap_or_else(|| Path::new("."));
    run(config, base)
}

/// The BAB + RATTLE step loop (SAD.md §2.3, §12.5). Sequencing per step:
/// half_kick -> geodesic_drift -> rebuild-if-stale -> compute_forces ->
/// half_kick -> constrain_velocities, with an initial constrain_velocities
/// before the loop so E(0) is measured against a constraint-consistent
/// velocity (see memory.md's energy-conservation lesson).
fn run(config: Config, base: &Path) -> Result<RunSummary, SimError> {
    if config.backend != Backend::Cpu {
        return Err(SimError::Config(geodesic_core::ConfigError::InvalidValue {
            key: "run.backend".to_string(),
            value: format!("{:?}", config.backend).to_lowercase(),
            reason: "only the \"cpu\" backend is available in this build (gpu lands in v0.5, hybrid in v0.6)".to_string(),
        }));
    }

    let prmtop_path = resolve(base, &config.prmtop);
    let inpcrd_path = resolve(base, &config.inpcrd);
    let trajectory_path = resolve(base, &config.trajectory);
    let energy_log_path = resolve(base, &config.energy_log);

    let prmtop_text = read_file(&prmtop_path)?;
    let (atoms, mut topology) = prmtop::parse(&prmtop_text)?;
    let n_atoms = atoms.mass.len();
    let inpcrd_text = read_file(&inpcrd_path)?;
    let mut state = inpcrd::parse(&inpcrd_text, n_atoms, false)?;

    if config.constrain_hydrogens {
        promote_hydrogen_bonds(&mut topology, &atoms);
    }

    let dt = config.dt;
    let n_steps = config.n_steps;
    let frame_interval = config.frame_interval as u64;
    let max_iter = config.max_constr_iter;
    let tol = config.constr_tol;

    let params = SimParams {
        n_atoms,
        n_steps,
        dt,
        box_size: config.box_size,
        periodic: config.periodic,
        r_cutoff: config.r_cutoff,
        r_skin: config.r_skin,
        r_switch: config.r_switch,
        max_constr_iter: max_iter,
        constr_tol: tol,
        frame_interval: config.frame_interval,
        n_threads: config.n_threads,
        total_energy: config.total_energy,
    };

    let mut backend: Box<dyn ComputeBackend> = Box::new(CpuBackend::new(atoms, topology, &params));

    let mut dcd = DcdWriter::create(&trajectory_path, n_atoms, config.frame_interval, dt)?;
    let mut csv = EnergyLogWriter::create(&energy_log_path, n_atoms)?;

    let mut f32x = vec![0.0f32; n_atoms];
    let mut f32y = vec![0.0f32; n_atoms];
    let mut f32z = vec![0.0f32; n_atoms];

    // Initial force evaluation, so the first half_kick uses forces at the
    // starting configuration, and the initial velocity projection onto the
    // constraint tangent space (memory.md).
    backend.build_neighbor_list(&mut state, &params);
    {
        let f = backend.compute_forces(&state);
        state.force_x.copy_from_slice(&f.fx);
        state.force_y.copy_from_slice(&f.fy);
        state.force_z.copy_from_slice(&f.fz);
    }
    state.potential_energy = backend.potential_energy();
    check_forces_finite(&state, 0)?;
    constrain_velocities(
        backend.topology(),
        backend.atoms(),
        &state.pos_x,
        &state.pos_y,
        &state.pos_z,
        &mut state.vel_x,
        &mut state.vel_y,
        &mut state.vel_z,
        max_iter,
        tol,
        0,
    )?;
    state.kinetic_energy = kinetic_energy(&state, &backend.atoms().mass);

    let mut n_frames = 0usize;
    write_frame(&mut dcd, &state, &mut f32x, &mut f32y, &mut f32z)?;
    csv.write_row(0, 0.0, state.potential_energy, state.kinetic_energy)?;
    n_frames += 1;

    let e0 = state.potential_energy + state.kinetic_energy;

    for _ in 0..n_steps {
        half_kick(&mut state, backend.atoms(), dt / 2.0);
        backend.geodesic_drift(&mut state, dt)?;
        check_positions_finite(&state, state.step)?;

        if backend.needs_rebuild(&state) {
            backend.build_neighbor_list(&mut state, &params);
        }
        {
            let f = backend.compute_forces(&state);
            state.force_x.copy_from_slice(&f.fx);
            state.force_y.copy_from_slice(&f.fy);
            state.force_z.copy_from_slice(&f.fz);
        }
        state.potential_energy = backend.potential_energy();
        check_forces_finite(&state, state.step)?;

        half_kick(&mut state, backend.atoms(), dt / 2.0);
        constrain_velocities(
            backend.topology(),
            backend.atoms(),
            &state.pos_x,
            &state.pos_y,
            &state.pos_z,
            &mut state.vel_x,
            &mut state.vel_y,
            &mut state.vel_z,
            max_iter,
            tol,
            state.step,
        )?;
        state.step += 1;

        if state.step % frame_interval == 0 {
            state.kinetic_energy = kinetic_energy(&state, &backend.atoms().mass);
            let time_ps = state.step as f64 * dt;
            write_frame(&mut dcd, &state, &mut f32x, &mut f32y, &mut f32z)?;
            csv.write_row(state.step, time_ps, state.potential_energy, state.kinetic_energy)?;
            n_frames += 1;
            check_energy_drift(&config, &state, e0, state.step)?;
        }
    }

    dcd.close()?;

    Ok(RunSummary {
        n_steps,
        n_frames,
        final_potential: state.potential_energy,
        final_kinetic: state.kinetic_energy,
        trajectory: trajectory_path,
        energy_log: energy_log_path,
    })
}

fn write_frame(
    dcd: &mut DcdWriter,
    state: &SimState,
    f32x: &mut [f32],
    f32y: &mut [f32],
    f32z: &mut [f32],
) -> Result<(), SimError> {
    for i in 0..state.pos_x.len() {
        f32x[i] = state.pos_x[i] as f32;
        f32y[i] = state.pos_y[i] as f32;
        f32z[i] = state.pos_z[i] as f32;
    }
    dcd.write_frame(f32x, f32y, f32z)?;
    Ok(())
}

/// SAD.md §12.4: warn (or stop, if configured) when |E(t) - E(0)| exceeds the
/// threshold. Absent a `[monitoring]` block, energy is not monitored.
fn check_energy_drift(config: &Config, state: &SimState, e0: f64, step: u64) -> Result<(), SimError> {
    let Some(m) = &config.monitoring else { return Ok(()) };
    let e = state.potential_energy + state.kinetic_energy;
    let drift = (e - e0).abs();
    if drift <= m.energy_drift_threshold_kcal {
        return Ok(());
    }
    match m.energy_drift_action {
        DriftAction::Warn => {
            eprintln!(
                "warning: energy drift {drift:.4} kcal/mol exceeds threshold {:.4} at step {step}",
                m.energy_drift_threshold_kcal
            );
            Ok(())
        }
        DriftAction::Stop => Err(SimError::Numerical(NumericalError::EnergyDrift {
            step,
            drift_kcal: drift,
            threshold_kcal: m.energy_drift_threshold_kcal,
        })),
    }
}
