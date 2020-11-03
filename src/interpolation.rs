use glam::f32::Vec3A;

///
/// Implements spherical interpolation along the great arc created by
/// the initial points. This returns a new point `p` percent of the way
/// along that arc.
///
/// Note: `a` and `b` should both be normalized for normalized results.
///
pub fn geometric_slerp(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
    let angle = a.dot(b).acos();

    let sin = angle.sin().recip();
    a * (((1.0 - p) * angle).sin() * sin) + b * ((p * angle).sin() * sin)
}

///
/// This is an optimization for the `geometric_slerp` in the case where `p`
/// is `0.5` or 50%.
///
/// Note: `a` and `b` should both be normalized for normalized results.
///
pub fn geometric_slerp_half(a: Vec3A, b: Vec3A) -> Vec3A {
    let angle = a.dot(b).acos();

    let sin_denom = angle.sin().recip();
    let sin_numer = (angle * 0.5).sin();

    (a + b) * sin_denom * sin_numer
}

///
/// This is an optimization for the case where multiple points require the
/// calculation of varying values of `p` for the same start and end points.
///
/// See the intended use in [`BaseShape::interpolate_multiple`].
///
/// Note: `a` and `b` should both be normalized for normalized results.
///
pub fn geometric_slerp_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
    let angle = a.dot(b).acos();
    let sin = angle.sin().recip();

    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] =
            a * (((1.0 - percent) * angle).sin() * sin) + b * ((percent * angle).sin() * sin);
    }
}

///
/// Performs normalized linear interpolation. This creates distortion when
/// compared with spherical interpolation along an arc, however this is most
/// likely faster, as though this avoids expensive sin and acos calculations.
///
pub fn normalized_lerp(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
    ((1.0 - p) * a + p * b).normalize()
}

///
/// This is an optimization of `normalized_lerp` which avoids a multiplication.
///
pub fn normalized_lerp_half(a: Vec3A, b: Vec3A) -> Vec3A {
    ((a + b) * 0.5).normalize()
}

///
/// This is provided as a plug in for people who need it, but this implements
/// essentially the same algorithm as `BaseShape` would without ever being
/// reimplemented.
///
pub fn normalized_lerp_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] = ((1.0 - percent) * a + percent * b).normalize();
    }
}

///
/// Simple linear interpolation. No weirdness here.
///
pub fn lerp(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
    (1.0 - p) * a + p * b
}

///
/// Gives the average of the two points.
///
pub fn lerp_half(a: Vec3A, b: Vec3A) -> Vec3A {
    (a + b) * 0.5
}

///
/// This is provided as a plug in for people who need it, but this implements
/// essentially the same algorithm as `BaseShape` would without ever being
/// reimplemented.
///
pub fn lerp_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] = (1.0 - percent) * a + percent * b;
    }
}
