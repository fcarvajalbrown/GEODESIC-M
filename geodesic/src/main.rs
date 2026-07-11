use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand};

/// GEODESIC-M: deterministic, physics-based exploration of protein
/// conformational manifolds (SAD.md §1). M1 ships the CPU headless engine.
#[derive(Parser)]
#[command(name = "geodesic", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Evaluate the initial potential energy of a system, to fill
    /// config.toml's `total_energy` before a run (SAD.md §9.2).
    Energy {
        /// AMBER prmtop (topology + force field).
        prmtop: PathBuf,
        /// AMBER inpcrd (initial coordinates).
        inpcrd: PathBuf,
    },
    /// Run a simulation from a config.toml, writing a DCD trajectory and an
    /// energy CSV (SAD.md §9.2).
    Run {
        /// TOML configuration file.
        config: PathBuf,
    },
}

fn main() -> ExitCode {
    // Pin OpenBLAS to one thread at process start, before any crate that might
    // touch LAPACK initializes (SAD.md §11.3) -- harmless in the M1 CPU build,
    // required for a deterministic PSL eigenvalue solve once topology lands.
    std::env::set_var("OPENBLAS_NUM_THREADS", "1");

    let cli = Cli::parse();
    match run(cli.command) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(command: Command) -> Result<(), geodesic_core::SimError> {
    match command {
        Command::Energy { prmtop, inpcrd } => {
            let r = geodesic::energy_from_files(&prmtop, &inpcrd)?;
            println!("potential energy (kcal/mol), r_switch={:.1} A, r_cutoff={:.1} A:", r.r_switch, r.r_cutoff);
            println!("  bond       {:>16.4}", r.bond);
            println!("  angle      {:>16.4}", r.angle);
            println!("  dihedral   {:>16.4}", r.dihedral);
            println!("  non-bonded {:>16.4}", r.nonbonded);
            println!("  total      {:>16.4}", r.total);
            Ok(())
        }
        Command::Run { config } => {
            let summary = geodesic::run_from_config_file(&config)?;
            println!(
                "run complete: {} steps, {} frames written",
                summary.n_steps, summary.n_frames
            );
            println!("  trajectory {}", summary.trajectory.display());
            println!("  energy log {}", summary.energy_log.display());
            println!(
                "  final V={:.4} kcal/mol  KE={:.4} kcal/mol",
                summary.final_potential, summary.final_kinetic
            );
            Ok(())
        }
    }
}
