use std::io::Read;

// Minimal manual DCD reader, independent of the writer implementation,
// used only to verify dcd.rs actually produces a structurally valid file.
struct RawDcd {
    icntrl: [i32; 20],
    natom: i32,
    frames: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)>,
}

fn read_i32(bytes: &[u8], pos: &mut usize) -> i32 {
    let v = i32::from_le_bytes(bytes[*pos..*pos + 4].try_into().unwrap());
    *pos += 4;
    v
}

fn read_f32_array(bytes: &[u8], pos: &mut usize) -> Vec<f32> {
    let marker = read_i32(bytes, pos) as usize;
    let n = marker / 4;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        out.push(f32::from_le_bytes(bytes[*pos..*pos + 4].try_into().unwrap()));
        *pos += 4;
    }
    let end_marker = read_i32(bytes, pos) as usize;
    assert_eq!(marker, end_marker, "record marker mismatch");
    out
}

fn parse_raw_dcd(bytes: &[u8]) -> RawDcd {
    let mut pos = 0usize;
    let open = read_i32(bytes, &mut pos);
    assert_eq!(open, 84);
    assert_eq!(&bytes[pos..pos + 4], b"CORD");
    pos += 4;
    let mut icntrl = [0i32; 20];
    for slot in icntrl.iter_mut() {
        *slot = read_i32(bytes, &mut pos);
    }
    let close = read_i32(bytes, &mut pos);
    assert_eq!(close, 84);

    let title_open = read_i32(bytes, &mut pos);
    let ntitle = read_i32(bytes, &mut pos);
    pos += (ntitle as usize) * 80;
    let title_close = read_i32(bytes, &mut pos);
    assert_eq!(title_open, title_close);

    let natom_open = read_i32(bytes, &mut pos);
    assert_eq!(natom_open, 4);
    let natom = read_i32(bytes, &mut pos);
    let natom_close = read_i32(bytes, &mut pos);
    assert_eq!(natom_close, 4);

    let mut frames = Vec::new();
    while pos < bytes.len() {
        let x = read_f32_array(bytes, &mut pos);
        let y = read_f32_array(bytes, &mut pos);
        let z = read_f32_array(bytes, &mut pos);
        frames.push((x, y, z));
    }

    RawDcd { icntrl, natom, frames }
}

#[test]
fn header_and_frames_round_trip() {
    let path = std::env::temp_dir().join("geodesic_dcd_test.dcd");
    let frame_interval = 500u32;
    let dt_ps = 0.004;

    let mut writer = geodesic_io::dcd::DcdWriter::create(&path, 2, frame_interval, dt_ps).unwrap();
    writer.write_frame(&[0.0, 3.0], &[0.0, 0.0], &[0.0, 0.0]).unwrap();
    writer.write_frame(&[0.1, 3.1], &[0.0, 0.0], &[0.0, 0.0]).unwrap();
    writer.close().unwrap();

    let mut bytes = Vec::new();
    std::fs::File::open(&path).unwrap().read_to_end(&mut bytes).unwrap();
    let dcd = parse_raw_dcd(&bytes);

    assert_eq!(dcd.icntrl[0], 2, "NSET should equal frames written");
    assert_eq!(dcd.icntrl[1], 0, "ISTART");
    assert_eq!(dcd.icntrl[2], frame_interval as i32, "NSAVC");
    let expected_delta = (dt_ps * 20.455) as f32;
    let got_delta = f32::from_le_bytes(dcd.icntrl[9].to_le_bytes());
    assert!((got_delta - expected_delta).abs() < 1e-6);
    assert_eq!(dcd.natom, 2);

    assert_eq!(dcd.frames.len(), 2);
    assert_eq!(dcd.frames[0].0, vec![0.0, 3.0]);
    assert_eq!(dcd.frames[1].0, vec![0.1, 3.1]);

    std::fs::remove_file(&path).ok();
}
