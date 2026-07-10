use geodesic_core::IoError;
use serde::Serialize;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Molar gas constant R in kcal/(mol·K) — CODATA R = 8.314462618 J/(mol·K)
/// converted via 1 kcal = 4184 J.
const KB_KCAL_PER_MOL_K: f64 = 0.0019872041;

fn io_err(path: &Path, source: std::io::Error) -> IoError {
    IoError {
        path: path.to_path_buf(),
        source,
    }
}

fn temperature_kelvin(kinetic_kcal: f64, n_atoms: usize) -> f64 {
    2.0 * kinetic_kcal / (3.0 * n_atoms as f64 * KB_KCAL_PER_MOL_K)
}

/// Appends one CSV line per call, flushing immediately so the file is
/// safe to tail during a run (SAD.md §10.4). Never buffers rows in memory.
pub struct EnergyLogWriter {
    file: File,
    path: PathBuf,
    n_atoms: usize,
}

impl EnergyLogWriter {
    pub fn create(path: &Path, n_atoms: usize) -> Result<Self, IoError> {
        let mut file = File::create(path).map_err(|e| io_err(path, e))?;
        writeln!(file, "step,time_ps,potential_kcal,kinetic_kcal,total_kcal,temperature_K")
            .map_err(|e| io_err(path, e))?;
        Ok(EnergyLogWriter {
            file,
            path: path.to_path_buf(),
            n_atoms,
        })
    }

    pub fn write_row(&mut self, step: u64, time_ps: f64, potential_kcal: f64, kinetic_kcal: f64) -> Result<(), IoError> {
        let total_kcal = potential_kcal + kinetic_kcal;
        let temperature_k = temperature_kelvin(kinetic_kcal, self.n_atoms);
        writeln!(
            self.file,
            "{step},{time_ps:.3},{potential_kcal:.3},{kinetic_kcal:.3},{total_kcal:.3},{temperature_k:.2}"
        )
        .map_err(|e| io_err(&self.path, e))?;
        self.file.flush().map_err(|e| io_err(&self.path, e))
    }
}

#[derive(Serialize)]
pub struct BarcodeMetadata {
    pub n_atoms: usize,
    pub n_frames: usize,
    pub frame_interval: u32,
}

#[derive(Serialize)]
pub struct BarcodeEntry {
    pub dim: u8,
    pub birth: f64,
    /// Infinite bars are encoded as -1.0, never as f64::INFINITY (not
    /// representable in JSON) — construct with `BarcodeEntry::infinite`.
    pub death: f64,
}

impl BarcodeEntry {
    pub fn finite(dim: u8, birth: f64, death: f64) -> Self {
        BarcodeEntry { dim, birth, death }
    }

    pub fn infinite(dim: u8, birth: f64) -> Self {
        BarcodeEntry { dim, birth, death: -1.0 }
    }
}

#[derive(Serialize)]
pub struct Barcode {
    pub metadata: BarcodeMetadata,
    pub barcode: Vec<BarcodeEntry>,
}

pub fn write_barcode(path: &Path, barcode: &Barcode) -> Result<(), IoError> {
    let json = serde_json::to_string_pretty(barcode)
        .map_err(|e| io_err(path, std::io::Error::new(std::io::ErrorKind::InvalidData, e)))?;
    std::fs::write(path, json).map_err(|e| io_err(path, e))
}
