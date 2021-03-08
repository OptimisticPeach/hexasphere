use crate::Vec3;

pub fn geometric_slerp<V: Vec3>(a: V, b: V, p: f32) -> V {
    let angle = a.dot(b).acos();

    let sin = angle.sin().recip();
    a * (((1.0 - p) * angle).sin() * sin) + b * ((p * angle).sin() * sin)
}

pub fn geometric_slerp_half<V: Vec3>(a: V, b: V) -> V {
    (a + b) * (2.0 * (1.0 + a.dot(b))).sqrt().recip()
}

pub fn geometric_slerp_multiple<V: Vec3>(a: V, b: V, indices: &[u32], points: &mut [V]) {
    let angle = a.dot(b).acos();
    let sin = angle.sin().recip();

    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] =
            a * (((1.0 - percent) * angle).sin() * sin) + b * ((percent * angle).sin() * sin);
    }
}
