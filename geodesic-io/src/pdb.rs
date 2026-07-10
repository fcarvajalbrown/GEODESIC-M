use geodesic_core::{AtomMeta, ConfigError, Element, IoError, SimState};
use std::io::Write;
use std::path::Path;

fn io_err(path: &Path, source: std::io::Error) -> IoError {
    IoError {
        path: path.to_path_buf(),
        source,
    }
}

fn infer_element(name: &str) -> Element {
    for c in name.chars() {
        if c.is_ascii_alphabetic() {
            return match c.to_ascii_uppercase() {
                'H' => Element::H,
                'C' => Element::C,
                'N' => Element::N,
                'O' => Element::O,
                'S' => Element::S,
                _ => Element::Unknown,
            };
        }
    }
    Element::Unknown
}

fn field(line: &str, start: usize, end: usize) -> &str {
    let end = end.min(line.len());
    if start >= line.len() {
        ""
    } else {
        &line[start..end]
    }
}

/// Secondary PDB-only input mode (SAD.md §10.1): positions and metadata
/// only. There is no bonded topology, no LJ/mass parameters, and no
/// force-field parameter table here — a bare PDB doesn't carry them, and
/// inventing default per-element values would risk silently wrong
/// physics. Callers combine this with their own parameters (test
/// fixtures, coarse-grained models) rather than treating it as a
/// complete simulation input.
pub fn parse_positions(text: &str) -> Result<(SimState, Vec<AtomMeta>), ConfigError> {
    let mut pos_x = Vec::new();
    let mut pos_y = Vec::new();
    let mut pos_z = Vec::new();
    let mut meta = Vec::new();

    for line in text.lines() {
        let record = field(line, 0, 6).trim();
        if record != "ATOM" && record != "HETATM" {
            continue;
        }

        let name = field(line, 12, 16).trim();
        let res_name = field(line, 17, 20).trim();
        let _ = res_name; // not represented in AtomMeta; kept for future use
        let chain = field(line, 21, 22).bytes().next().unwrap_or(b' ');
        let res_seq = field(line, 22, 26).trim();
        let x = field(line, 30, 38).trim();
        let y = field(line, 38, 46).trim();
        let z = field(line, 46, 54).trim();
        let element_col = field(line, 76, 78).trim();

        let parse_coord = |s: &str, axis: &str| -> Result<f64, ConfigError> {
            s.parse::<f64>().map_err(|_| ConfigError::InvalidValue {
                key: format!("PDB {axis} coordinate"),
                value: s.to_string(),
                reason: "expected a float".to_string(),
            })
        };
        pos_x.push(parse_coord(x, "x")?);
        pos_y.push(parse_coord(y, "y")?);
        pos_z.push(parse_coord(z, "z")?);

        let residue_id: u32 = res_seq.parse().unwrap_or(0);
        let element = if element_col.is_empty() {
            infer_element(name)
        } else {
            infer_element(element_col)
        };
        let mut atom_name = [b' '; 4];
        let name_bytes = name.as_bytes();
        let n = name_bytes.len().min(4);
        atom_name[..n].copy_from_slice(&name_bytes[..n]);

        meta.push(AtomMeta {
            element,
            residue_id,
            atom_name,
            chain_id: chain,
        });
    }

    let natom = pos_x.len();
    let mut state = SimState::new(natom);
    state.pos_x = pos_x;
    state.pos_y = pos_y;
    state.pos_z = pos_z;

    Ok((state, meta))
}

/// On-demand single-frame snapshot (SAD.md §10.5). B-factor column holds
/// the PSL flexibility score if available, 0.00 otherwise, so VMD/PyMOL
/// can color by flexibility without post-processing.
pub fn write_snapshot(
    path: &Path,
    pos_x: &[f64],
    pos_y: &[f64],
    pos_z: &[f64],
    meta: &[AtomMeta],
    flexibility: Option<&[f32]>,
) -> Result<(), IoError> {
    let mut file = std::fs::File::create(path).map_err(|e| io_err(path, e))?;
    for i in 0..pos_x.len() {
        let name = String::from_utf8_lossy(&meta[i].atom_name);
        let name = name.trim_end();
        let chain = meta[i].chain_id as char;
        let bfactor = flexibility.map(|f| f[i]).unwrap_or(0.0);
        let element = element_symbol(meta[i].element);
        writeln!(
            file,
            "ATOM  {:>5} {:<4} {:<3} {}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            i + 1,
            name,
            "RES",
            chain,
            meta[i].residue_id,
            pos_x[i],
            pos_y[i],
            pos_z[i],
            1.00,
            bfactor,
            element,
        )
        .map_err(|e| io_err(path, e))?;
    }
    writeln!(file, "END").map_err(|e| io_err(path, e))?;
    Ok(())
}

fn element_symbol(e: Element) -> &'static str {
    match e {
        Element::H => "H",
        Element::C => "C",
        Element::N => "N",
        Element::O => "O",
        Element::S => "S",
        Element::Unknown => "X",
    }
}
