use crate::fixed_width::chunk_fields;
use geodesic_core::{ConfigError, SimState};

/// AMBER restart/inpcrd velocities are stored in units of Å per (1/20.455 ps)
/// — multiply by this factor to recover Å/ps (matches the AKMA time unit
/// used elsewhere in the AMBER/CHARMM ecosystem: 1/20.455 ps ≈ 48.888 fs).
const AMBER_VELOCITY_TO_ANG_PER_PS: f64 = 20.455;

/// `has_box` disambiguates the optional trailing velocity block from the
/// optional trailing periodic-box line — both are 1 line of 6F12.7 for a
/// small enough system, so line count alone can't tell them apart.
pub fn parse(text: &str, expected_n_atoms: usize, has_box: bool) -> Result<SimState, ConfigError> {
    let all_lines: Vec<&str> = text.lines().collect();
    if all_lines.len() < 2 {
        return Err(ConfigError::InvalidValue {
            key: "inpcrd".to_string(),
            value: String::new(),
            reason: "file has fewer than 2 lines (missing title or atom-count header)".to_string(),
        });
    }

    let header = all_lines[1];
    let natom_field = chunk_fields(header, 5, 1)[0];
    let natom: usize = natom_field.trim().parse().map_err(|_| ConfigError::InvalidValue {
        key: "inpcrd header".to_string(),
        value: natom_field.to_string(),
        reason: "expected an atom count".to_string(),
    })?;

    if natom != expected_n_atoms {
        return Err(ConfigError::PhysicallyInvalid {
            description: format!(
                "inpcrd declares {natom} atoms but prmtop declares {expected_n_atoms}"
            ),
        });
    }

    let n_coord_values = natom * 3;
    let n_coord_lines = n_coord_values.div_ceil(6);
    let data_lines = &all_lines[2..];
    if data_lines.len() < n_coord_lines {
        return Err(ConfigError::InvalidValue {
            key: "inpcrd".to_string(),
            value: String::new(),
            reason: format!(
                "expected {n_coord_lines} coordinate line(s) for {natom} atoms, found {}",
                data_lines.len()
            ),
        });
    }

    let coord_text: String = data_lines[..n_coord_lines].concat();
    let coords = parse_f64_block(&coord_text, n_coord_values, "inpcrd coordinates")?;

    let mut state = SimState::new(natom);
    for i in 0..natom {
        state.pos_x[i] = coords[3 * i];
        state.pos_y[i] = coords[3 * i + 1];
        state.pos_z[i] = coords[3 * i + 2];
    }

    let remaining = &data_lines[n_coord_lines..];
    let n_box_lines = usize::from(has_box);
    let velocity_lines_available = remaining.len().saturating_sub(n_box_lines);

    if velocity_lines_available >= n_coord_lines && n_coord_lines > 0 {
        let vel_text: String = remaining[..n_coord_lines].concat();
        let vels = parse_f64_block(&vel_text, n_coord_values, "inpcrd velocities")?;
        for i in 0..natom {
            state.vel_x[i] = vels[3 * i] * AMBER_VELOCITY_TO_ANG_PER_PS;
            state.vel_y[i] = vels[3 * i + 1] * AMBER_VELOCITY_TO_ANG_PER_PS;
            state.vel_z[i] = vels[3 * i + 2] * AMBER_VELOCITY_TO_ANG_PER_PS;
        }
    }

    Ok(state)
}

fn parse_f64_block(text: &str, count: usize, context: &str) -> Result<Vec<f64>, ConfigError> {
    chunk_fields(text, 12, count)
        .into_iter()
        .map(|field| {
            field.trim().parse::<f64>().map_err(|_| ConfigError::InvalidValue {
                key: context.to_string(),
                value: field.to_string(),
                reason: "expected a float".to_string(),
            })
        })
        .collect()
}
