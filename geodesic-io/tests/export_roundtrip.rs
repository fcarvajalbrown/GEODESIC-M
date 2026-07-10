use geodesic_io::export::{Barcode, BarcodeEntry, BarcodeMetadata, EnergyLogWriter};

#[test]
fn energy_log_format_and_temperature() {
    let path = std::env::temp_dir().join("geodesic_energy_test.csv");
    {
        let mut writer = EnergyLogWriter::create(&path, 100).unwrap();
        // Ek = 1.5 * N * kB * T  =>  T = 2*Ek / (3*N*kB); pick Ek so T comes out round
        // 3*100*0.0019872041 = 0.59616123; Ek = 0.59616123 * 150 => T = 300 K
        writer.write_row(0, 0.0, -100.0, 0.59616123 * 150.0).unwrap();
    }
    let contents = std::fs::read_to_string(&path).unwrap();
    let mut lines = contents.lines();
    assert_eq!(
        lines.next().unwrap(),
        "step,time_ps,potential_kcal,kinetic_kcal,total_kcal,temperature_K"
    );
    let row = lines.next().unwrap();
    let fields: Vec<&str> = row.split(',').collect();
    assert_eq!(fields[0], "0");
    assert_eq!(fields[1], "0.000");
    assert_eq!(fields[5], "300.00");

    std::fs::remove_file(&path).ok();
}

#[test]
fn barcode_json_schema_and_infinite_encoding() {
    let path = std::env::temp_dir().join("geodesic_barcode_test.json");
    let barcode = Barcode {
        metadata: BarcodeMetadata {
            n_atoms: 2,
            n_frames: 10,
            frame_interval: 500,
        },
        barcode: vec![
            BarcodeEntry::finite(0, 0.0, 142.5),
            BarcodeEntry::infinite(1, 18.3),
        ],
    };
    geodesic_io::export::write_barcode(&path, &barcode).unwrap();

    let text = std::fs::read_to_string(&path).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(parsed["metadata"]["n_atoms"], 2);
    assert_eq!(parsed["barcode"][0]["death"], 142.5);
    assert_eq!(parsed["barcode"][1]["death"], -1.0);

    std::fs::remove_file(&path).ok();
}
