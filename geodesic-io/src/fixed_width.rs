/// Slices `raw` into `count` fields of exactly `width` characters each.
/// Fortran fixed-width records can abut with no separating space, so
/// whitespace-splitting is not safe here.
pub(crate) fn chunk_fields(raw: &str, width: usize, count: usize) -> Vec<&str> {
    let bytes = raw.as_bytes();
    let mut fields = Vec::with_capacity(count);
    for i in 0..count {
        let start = i * width;
        let end = (start + width).min(bytes.len());
        let field = if start < bytes.len() { &raw[start..end] } else { "" };
        fields.push(field);
    }
    fields
}
