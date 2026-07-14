use geodesic_gpu::device;

#[test]
fn context_creation_is_infallible_or_cleanly_absent() {
    match device::try_new() {
        Ok(_) => {}
        Err(geodesic_core::BackendError::NoAdapter) => {
            eprintln!("skipping: no GPU adapter (DX12/Vulkan) available");
        }
        Err(e) => panic!("unexpected backend error: {e}"),
    }
}
