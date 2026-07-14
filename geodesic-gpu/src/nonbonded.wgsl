struct Params {
  n_atoms: u32,
  r_cutoff: f32,
  r_switch: f32,
  _pad0: f32,
  box_size: vec3<f32>,
  _pad1: f32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> positions: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> sigma: array<f32>;
@group(0) @binding(3) var<storage, read> epsilon: array<f32>;
@group(0) @binding(4) var<storage, read> offsets: array<u32>;
@group(0) @binding(5) var<storage, read> neighbors: array<u32>;
@group(0) @binding(6) var<storage, read_write> out_force: array<vec4<f32>>;

fn min_image(d: f32, box_len: f32) -> f32 {
  return d - box_len * round(d / box_len);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let i = gid.x;
  if (i >= params.n_atoms) { return; }
  let pi = positions[i].xyz;
  let si = sigma[i];
  let ei = epsilon[i];
  let rc2 = params.r_cutoff * params.r_cutoff;
  var acc = vec3<f32>(0.0, 0.0, 0.0);
  var energy = 0.0;
  let start = offsets[i];
  let end = offsets[i + 1u];
  for (var k = start; k < end; k = k + 1u) {
    let j = neighbors[k];
    let pj = positions[j].xyz;
    let dx = min_image(pj.x - pi.x, params.box_size.x);
    let dy = min_image(pj.y - pi.y, params.box_size.y);
    let dz = min_image(pj.z - pi.z, params.box_size.z);
    let r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > rc2 || r2 == 0.0) { continue; }
    let r = sqrt(r2);
    let sig = 0.5 * (si + sigma[j]);
    let eps = sqrt(ei * epsilon[j]);
    if (eps == 0.0) { continue; }
    let sr = sig / r;
    let sr2 = sr * sr;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    let v_lj = 4.0 * eps * (sr12 - sr6);
    let f_lj = 24.0 * eps / r * (2.0 * sr12 - sr6);
    var v = v_lj;
    var f_radial = f_lj;
    if (r > params.r_switch) {
      let denom = params.r_cutoff - params.r_switch;
      let u = (r - params.r_switch) / denom;
      let u2 = u * u;
      let s = 1.0 - 10.0 * u2 * u + 15.0 * u2 * u2 - 6.0 * u2 * u2 * u;
      let ds_dr = -30.0 * u2 * (1.0 - u) * (1.0 - u) / denom;
      v = v_lj * s;
      f_radial = f_lj * s - v_lj * ds_dr;
    }
    energy = energy + v;
    let inv_r = 1.0 / r;
    acc.x = acc.x - f_radial * dx * inv_r;
    acc.y = acc.y - f_radial * dy * inv_r;
    acc.z = acc.z - f_radial * dz * inv_r;
  }
  out_force[i] = vec4<f32>(acc.x, acc.y, acc.z, energy);
}
