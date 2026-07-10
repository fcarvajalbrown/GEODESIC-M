use crate::fixed_width::chunk_fields;
use geodesic_core::{AtomData, AtomMeta, BondedTopology, ConfigError, Element};
use std::collections::HashMap;

/// Raw %FLAG sections, keyed by flag name, holding the concatenated
/// fixed-width data text (newlines stripped) for that section.
struct Sections(HashMap<String, String>);

impl Sections {
    fn parse(text: &str) -> Sections {
        let mut sections: HashMap<String, String> = HashMap::new();
        let mut current: Option<String> = None;
        let mut skip_format_line = false;

        for line in text.lines() {
            if let Some(flag) = line.strip_prefix("%FLAG ") {
                let flag = flag.trim().to_string();
                sections.entry(flag.clone()).or_default();
                current = Some(flag);
                skip_format_line = true;
                continue;
            }
            if line.starts_with("%FORMAT") {
                skip_format_line = false;
                continue;
            }
            if line.starts_with('%') {
                // %VERSION, %COMMENT, or any other non-data directive
                continue;
            }
            if skip_format_line {
                // data line arrived before a %FORMAT line was seen; treat as data anyway
                skip_format_line = false;
            }
            if let Some(flag) = &current {
                sections.entry(flag.clone()).or_default().push_str(line);
            }
        }

        Sections(sections)
    }

    fn raw(&self, flag: &str) -> Result<&str, ConfigError> {
        self.0
            .get(flag)
            .map(|s| s.as_str())
            .ok_or_else(|| ConfigError::MissingRequired(format!("%FLAG {flag}")))
    }

    /// Fixed-width integer fields, `width` characters each.
    fn ints(&self, flag: &str, width: usize, count: usize) -> Result<Vec<i64>, ConfigError> {
        let raw = self.raw(flag)?;
        chunk_fields(raw, width, count)
            .into_iter()
            .map(|field| {
                field.trim().parse::<i64>().map_err(|_| ConfigError::InvalidValue {
                    key: format!("%FLAG {flag}"),
                    value: field.to_string(),
                    reason: "expected an integer".to_string(),
                })
            })
            .collect()
    }

    /// Fixed-width float fields, `width` characters each.
    fn floats(&self, flag: &str, width: usize, count: usize) -> Result<Vec<f64>, ConfigError> {
        let raw = self.raw(flag)?;
        chunk_fields(raw, width, count)
            .into_iter()
            .map(|field| {
                field.trim().parse::<f64>().map_err(|_| ConfigError::InvalidValue {
                    key: format!("%FLAG {flag}"),
                    value: field.to_string(),
                    reason: "expected a float".to_string(),
                })
            })
            .collect()
    }

    /// Fixed-width 4-character name fields.
    fn names4(&self, flag: &str, count: usize) -> Result<Vec<[u8; 4]>, ConfigError> {
        let raw = self.raw(flag)?;
        chunk_fields(raw, 4, count)
            .into_iter()
            .map(|field| {
                let bytes = field.as_bytes();
                let mut name = [b' '; 4];
                let n = bytes.len().min(4);
                name[..n].copy_from_slice(&bytes[..n]);
                Ok(name)
            })
            .collect()
    }
}

fn infer_element(name: &[u8; 4]) -> Element {
    for &b in name.iter() {
        let c = b as char;
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

fn residue_id_for_atom(atom_idx_0based: usize, residue_pointer: &[i64]) -> u32 {
    let atom_1based = atom_idx_0based + 1;
    let mut res_id = 0usize;
    for (r, &start) in residue_pointer.iter().enumerate() {
        if (start as usize) <= atom_1based {
            res_id = r;
        } else {
            break;
        }
    }
    res_id as u32
}

/// Bonds/angles/dihedrals in prmtop encode atom indices as 3x the 0-based
/// coordinate offset (recovered via integer division by 3), and negative
/// signs on the 3rd/4th dihedral atom flag "skip 1-4" / "improper" — M1
/// doesn't compute a separate scaled 1-4 term (see BondedTopology::excl_i;
/// 1-4 pairs are fully excluded like 1-2/1-3), and impropers use the same
/// cosine-series formula as propers, so only the unsigned index is kept.
fn parse_bonds(raw: &[i64]) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut i = Vec::new();
    let mut j = Vec::new();
    let mut ty = Vec::new();
    for c in raw.chunks(3) {
        i.push((c[0] / 3) as u32);
        j.push((c[1] / 3) as u32);
        ty.push((c[2] - 1) as u32);
    }
    (i, j, ty)
}

fn parse_angles(raw: &[i64]) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut i = Vec::new();
    let mut j = Vec::new();
    let mut k = Vec::new();
    let mut ty = Vec::new();
    for c in raw.chunks(4) {
        i.push((c[0] / 3) as u32);
        j.push((c[1] / 3) as u32);
        k.push((c[2] / 3) as u32);
        ty.push((c[3] - 1) as u32);
    }
    (i, j, k, ty)
}

type DihedralIndices = (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>);

fn parse_dihedrals(raw: &[i64]) -> DihedralIndices {
    let mut i = Vec::new();
    let mut j = Vec::new();
    let mut k = Vec::new();
    let mut l = Vec::new();
    let mut ty = Vec::new();
    for c in raw.chunks(5) {
        i.push((c[0].abs() / 3) as u32);
        j.push((c[1].abs() / 3) as u32);
        k.push((c[2].abs() / 3) as u32);
        l.push((c[3].abs() / 3) as u32);
        ty.push((c[4] - 1) as u32);
    }
    (i, j, k, l, ty)
}

pub fn parse(text: &str) -> Result<(AtomData, BondedTopology), ConfigError> {
    let sections = Sections::parse(text);

    let pointers = sections.ints("POINTERS", 8, 20)?;
    let natom = pointers[0] as usize;
    let ntypes = pointers[1] as usize;
    let nbonh = pointers[2] as usize;
    let ntheth = pointers[4] as usize;
    let nphih = pointers[6] as usize;
    let nnb = pointers[10] as usize;
    let nres = pointers[11] as usize;
    let nbona = pointers[12] as usize;
    let ntheta = pointers[13] as usize;
    let nphia = pointers[14] as usize;
    let numbnd = pointers[15] as usize;
    let numang = pointers[16] as usize;
    let nptra = pointers[17] as usize;

    let names = sections.names4("ATOM_NAME", natom)?;
    let charge = sections.floats("CHARGE", 16, natom)?;
    let mass = sections.floats("MASS", 16, natom)?;
    let atom_type_index = sections.ints("ATOM_TYPE_INDEX", 8, natom)?;
    let number_excluded = sections.ints("NUMBER_EXCLUDED_ATOMS", 8, natom)?;
    let nonbonded_parm_index = sections.ints("NONBONDED_PARM_INDEX", 8, ntypes * ntypes)?;
    let n_lj_types = ntypes * (ntypes + 1) / 2;
    let lj_acoef = sections.floats("LENNARD_JONES_ACOEF", 16, n_lj_types)?;
    let lj_bcoef = sections.floats("LENNARD_JONES_BCOEF", 16, n_lj_types)?;
    let residue_pointer = sections.ints("RESIDUE_POINTER", 8, nres)?;
    let excluded_atoms_list = sections.ints("EXCLUDED_ATOMS_LIST", 8, nnb)?;

    let mut epsilon = Vec::with_capacity(natom);
    let mut sigma = Vec::with_capacity(natom);
    for (i, &ty_raw) in atom_type_index.iter().enumerate() {
        let ty = ty_raw as usize;
        let ico = nonbonded_parm_index[ntypes * (ty - 1) + ty - 1];
        if ico <= 0 {
            return Err(ConfigError::PhysicallyInvalid {
                description: format!(
                    "atom {i} uses a 10-12 H-bond LJ type (NONBONDED_PARM_INDEX={ico}); not supported"
                ),
            });
        }
        let idx = (ico - 1) as usize;
        let (a, b) = (lj_acoef[idx], lj_bcoef[idx]);
        if a == 0.0 || b == 0.0 {
            epsilon.push(0.0);
            sigma.push(0.0);
        } else {
            sigma.push((a / b).powf(1.0 / 6.0));
            epsilon.push((b * b) / (4.0 * a));
        }
    }

    let mut meta = Vec::with_capacity(natom);
    for (i, &name) in names.iter().enumerate() {
        meta.push(AtomMeta {
            element: infer_element(&name),
            residue_id: residue_id_for_atom(i, &residue_pointer),
            atom_name: name,
            chain_id: 0,
        });
    }

    let atom_data = AtomData {
        epsilon,
        sigma,
        mass,
        charge,
        meta,
    };

    let bonds_h = sections.ints("BONDS_INC_HYDROGEN", 8, 3 * nbonh)?;
    let bonds_a = sections.ints("BONDS_WITHOUT_HYDROGEN", 8, 3 * nbona)?;
    let angles_h = sections.ints("ANGLES_INC_HYDROGEN", 8, 4 * ntheth)?;
    let angles_a = sections.ints("ANGLES_WITHOUT_HYDROGEN", 8, 4 * ntheta)?;
    let dihed_h = sections.ints("DIHEDRALS_INC_HYDROGEN", 8, 5 * nphih)?;
    let dihed_a = sections.ints("DIHEDRALS_WITHOUT_HYDROGEN", 8, 5 * nphia)?;

    let bond_force_k = sections.floats("BOND_FORCE_CONSTANT", 16, numbnd)?;
    let bond_r0 = sections.floats("BOND_EQUIL_VALUE", 16, numbnd)?;
    let angle_force_k = sections.floats("ANGLE_FORCE_CONSTANT", 16, numang)?;
    let angle_th0 = sections.floats("ANGLE_EQUIL_VALUE", 16, numang)?;
    let dihed_force_k = sections.floats("DIHEDRAL_FORCE_CONSTANT", 16, nptra)?;
    let dihed_period = sections.floats("DIHEDRAL_PERIODICITY", 16, nptra)?;
    let dihed_phase = sections.floats("DIHEDRAL_PHASE", 16, nptra)?;

    let (mut bond_i, mut bond_j, bond_ty1) = parse_bonds(&bonds_h);
    let (bi2, bj2, bond_ty2) = parse_bonds(&bonds_a);
    bond_i.extend(bi2);
    bond_j.extend(bj2);
    let bond_ty: Vec<u32> = bond_ty1.into_iter().chain(bond_ty2).collect();
    let bond_k: Vec<f64> = bond_ty.iter().map(|&t| bond_force_k[t as usize]).collect();
    let bond_r0v: Vec<f64> = bond_ty.iter().map(|&t| bond_r0[t as usize]).collect();

    let (mut angle_i, mut angle_j, mut angle_k, angle_ty1) = parse_angles(&angles_h);
    let (ai2, aj2, ak2, angle_ty2) = parse_angles(&angles_a);
    angle_i.extend(ai2);
    angle_j.extend(aj2);
    angle_k.extend(ak2);
    let angle_ty: Vec<u32> = angle_ty1.into_iter().chain(angle_ty2).collect();
    let angle_kth: Vec<f64> = angle_ty.iter().map(|&t| angle_force_k[t as usize]).collect();
    let angle_th0v: Vec<f64> = angle_ty.iter().map(|&t| angle_th0[t as usize]).collect();

    let (mut dihed_i, mut dihed_j, mut dihed_k, mut dihed_l, dihed_ty1) = parse_dihedrals(&dihed_h);
    let (di2, dj2, dk2, dl2, dihed_ty2) = parse_dihedrals(&dihed_a);
    dihed_i.extend(di2);
    dihed_j.extend(dj2);
    dihed_k.extend(dk2);
    dihed_l.extend(dl2);
    let dihed_ty: Vec<u32> = dihed_ty1.into_iter().chain(dihed_ty2).collect();
    let dihed_kphi: Vec<f64> = dihed_ty.iter().map(|&t| dihed_force_k[t as usize]).collect();
    let dihed_n: Vec<u32> = dihed_ty
        .iter()
        .map(|&t| dihed_period[t as usize].round() as u32)
        .collect();
    let dihed_delta: Vec<f64> = dihed_ty.iter().map(|&t| dihed_phase[t as usize]).collect();

    let mut excl_i = Vec::new();
    let mut excl_j = Vec::new();
    let mut cursor = 0usize;
    for (i, &count_raw) in number_excluded.iter().enumerate() {
        let count = count_raw as usize;
        for _ in 0..count {
            let raw = excluded_atoms_list[cursor];
            cursor += 1;
            if raw == 0 {
                continue;
            }
            excl_i.push(i as u32);
            excl_j.push((raw - 1) as u32);
        }
    }

    let bonded = BondedTopology {
        bond_i,
        bond_j,
        bond_k,
        bond_r0: bond_r0v,
        angle_i,
        angle_j,
        angle_k,
        angle_kth,
        angle_th0: angle_th0v,
        dihed_i,
        dihed_j,
        dihed_k,
        dihed_l,
        dihed_kphi,
        dihed_n,
        dihed_delta,
        constr_i: Vec::new(),
        constr_j: Vec::new(),
        constr_dsq: Vec::new(),
        excl_i,
        excl_j,
    };

    Ok((atom_data, bonded))
}
