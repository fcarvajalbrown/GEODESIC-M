use geodesic_core::IoError;
use std::fs::File;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// 1 AKMA = 1/20.455 ps (cross-checked against the inpcrd velocity unit,
/// which uses the same underlying AMBER/CHARMM time constant).
const AKMA_PER_PS: f64 = 20.455;

/// Byte offset of the NSET field within the 84-byte header block —
/// 4-byte record marker + 4-byte "CORD" magic, then ICNTRL[0].
const NSET_OFFSET: u64 = 8;

/// Writes CHARMM/NAMD-style DCD trajectory files (single-precision X/Y/Z,
/// no per-frame unit cell block). Frames are appended one at a time;
/// the frame count (NSET) is a placeholder until `close()` patches it.
pub struct DcdWriter {
    file: File,
    path: PathBuf,
    natom: usize,
    n_frames_written: u32,
    closed: bool,
}

impl DcdWriter {
    pub fn create(path: &Path, natom: usize, frame_interval: u32, dt_ps: f64) -> Result<Self, IoError> {
        let mut file = File::create(path).map_err(|e| io_err(path, e))?;
        write_header(&mut file, natom, frame_interval, dt_ps).map_err(|e| io_err(path, e))?;
        Ok(DcdWriter {
            file,
            path: path.to_path_buf(),
            natom,
            n_frames_written: 0,
            closed: false,
        })
    }

    pub fn write_frame(&mut self, pos_x: &[f32], pos_y: &[f32], pos_z: &[f32]) -> Result<(), IoError> {
        if pos_x.len() != self.natom || pos_y.len() != self.natom || pos_z.len() != self.natom {
            return Err(io_err(
                &self.path,
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!(
                        "frame has {}/{}/{} x/y/z values, writer expects {} atoms",
                        pos_x.len(),
                        pos_y.len(),
                        pos_z.len(),
                        self.natom
                    ),
                ),
            ));
        }
        write_frame_records(&mut self.file, pos_x, pos_y, pos_z).map_err(|e| io_err(&self.path, e))?;
        self.n_frames_written += 1;
        Ok(())
    }

    pub fn close(mut self) -> Result<(), IoError> {
        self.patch_nset()
    }

    fn patch_nset(&mut self) -> Result<(), IoError> {
        if self.closed {
            return Ok(());
        }
        self.file
            .seek(SeekFrom::Start(NSET_OFFSET))
            .map_err(|e| io_err(&self.path, e))?;
        self.file
            .write_all(&(self.n_frames_written as i32).to_le_bytes())
            .map_err(|e| io_err(&self.path, e))?;
        self.file.flush().map_err(|e| io_err(&self.path, e))?;
        self.closed = true;
        Ok(())
    }
}

impl Drop for DcdWriter {
    fn drop(&mut self) {
        // Best-effort safety net if `close()` was never called explicitly —
        // errors here can't be surfaced, `close()` is the way to observe them.
        let _ = self.patch_nset();
    }
}

fn io_err(path: &Path, source: std::io::Error) -> IoError {
    IoError {
        path: path.to_path_buf(),
        source,
    }
}

fn write_header(file: &mut File, natom: usize, frame_interval: u32, dt_ps: f64) -> std::io::Result<()> {
    let delta = (dt_ps * AKMA_PER_PS) as f32;

    file.write_all(&84i32.to_le_bytes())?;
    file.write_all(b"CORD")?;
    for i in 0..20 {
        match i {
            0 => file.write_all(&0i32.to_le_bytes())?, // NSET, patched on close
            1 => file.write_all(&0i32.to_le_bytes())?, // ISTART
            2 => file.write_all(&(frame_interval as i32).to_le_bytes())?, // NSAVC
            8 => file.write_all(&0i32.to_le_bytes())?, // NAMNF: 0 fixed atoms
            9 => file.write_all(&delta.to_le_bytes())?, // DELTA, AKMA units, f32
            19 => file.write_all(&24i32.to_le_bytes())?, // CHARMM version marker
            _ => file.write_all(&0i32.to_le_bytes())?,
        }
    }
    file.write_all(&84i32.to_le_bytes())?;

    let mut title = [b' '; 80];
    let msg = b"Created by GEODESIC-M";
    title[..msg.len()].copy_from_slice(msg);
    file.write_all(&84i32.to_le_bytes())?; // 4 (NTITLE) + 80 (one title line)
    file.write_all(&1i32.to_le_bytes())?;
    file.write_all(&title)?;
    file.write_all(&84i32.to_le_bytes())?;

    file.write_all(&4i32.to_le_bytes())?;
    file.write_all(&(natom as i32).to_le_bytes())?;
    file.write_all(&4i32.to_le_bytes())?;

    Ok(())
}

fn write_frame_records(file: &mut File, x: &[f32], y: &[f32], z: &[f32]) -> std::io::Result<()> {
    for arr in [x, y, z] {
        let nbytes = (arr.len() * 4) as i32;
        file.write_all(&nbytes.to_le_bytes())?;
        for &v in arr {
            file.write_all(&v.to_le_bytes())?;
        }
        file.write_all(&nbytes.to_le_bytes())?;
    }
    Ok(())
}
