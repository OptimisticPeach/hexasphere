#[cfg(feature = "shape-extras")]
use super::EquilateralBaseShape;
use super::{interpolation, BaseShape, Subdivided, Triangle};
use alloc::{boxed::Box, vec::Vec};
use glam::Vec3A;

///
/// Implements an icosahedron as the base shape.
///
/// - 12 vertices
/// - 20 faces
/// - 30 edges
///
/// This shape has the best results for a sphere.
///
/// The resulting smaller triangles are close to being
/// equilateral, so if one draws lines from the center
/// of the each triangle to the middle of the each edge
/// then the result will be 12 pentagons and many hexagons.
///
#[derive(Default, Copy, Clone, Debug)]
pub struct IcoSphereBase;

impl BaseShape for IcoSphereBase {
    #[inline]
    fn initial_points(&self) -> Vec<Vec3A> {
        consts::icosphere::INITIAL_POINTS.to_vec()
    }

    #[inline]
    fn triangles(&self) -> Box<[Triangle]> {
        Box::new(consts::icosphere::TRIANGLES)
    }
    const EDGES: usize = consts::icosphere::EDGES;

    #[inline]
    fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        interpolation::geometric_slerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
        interpolation::geometric_slerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        interpolation::geometric_slerp_multiple(a, b, indices, points);
    }
}

#[cfg(feature = "shape-extras")]
impl EquilateralBaseShape for IcoSphereBase {
    #[inline]
    fn triangle_normals() -> &'static [Vec3A] {
        &consts::icosphere::TRIANGLE_NORMALS
    }

    #[inline]
    fn triangle_min_dot() -> f32 {
        consts::icosphere::MIN_NORMAL_DOT
    }
}

///
/// Icosphere.
///
/// See [`IcoSphereBase`].
///
pub type IcoSphere<T> = Subdivided<T, IcoSphereBase>;

impl<T> IcoSphere<T> {
    ///
    /// Calculate distance from the center of a shape (pentagon or hexagon)
    /// to one of the vertices of the shape.
    ///
    /// In other words, the radius of the circumscribed circle.
    ///
    pub fn radius_shapes(&self) -> f32 {
        let subdivisions = self.subdivisions as f32 + 1.0;
        const DEFAULT_ANGLE: f32 = 1.10714871779409085306;
        let angle = DEFAULT_ANGLE / subdivisions;
        crate::math::sin(angle * 0.5) * 2.0
    }
}

///
/// Implements the same shape as [`IcoSphereBase`], however
/// it uses normalized linear interpolation, rather than
/// geometric spherical linear interpolation. (`nlerp` over `slerp`).
///
#[derive(Default, Copy, Clone, Debug)]
pub struct NormIcoSphereBase;

impl BaseShape for NormIcoSphereBase {
    #[inline]
    fn initial_points(&self) -> Vec<Vec3A> {
        consts::icosphere::INITIAL_POINTS.to_vec()
    }

    #[inline]
    fn triangles(&self) -> Box<[Triangle]> {
        Box::new(consts::icosphere::TRIANGLES)
    }
    const EDGES: usize = consts::icosphere::EDGES;

    #[inline]
    fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        interpolation::normalized_lerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
        interpolation::normalized_lerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        interpolation::normalized_lerp_multiple(a, b, indices, points);
    }
}

#[cfg(feature = "shape-extras")]
impl EquilateralBaseShape for NormIcoSphereBase {
    #[inline]
    fn triangle_normals() -> &'static [Vec3A] {
        &consts::icosphere::TRIANGLE_NORMALS
    }

    #[inline]
    fn triangle_min_dot() -> f32 {
        consts::icosphere::MIN_NORMAL_DOT
    }
}

///
/// Normalized Icosphere.
///
/// See [`NormIcoSphereBase`].
///
pub type NormIcoSphere<T> = Subdivided<T, NormIcoSphereBase>;

///
/// Implements a tetrahedron as the base shape.
///
/// - 4 vertices
/// - 4 faces
/// - 6 edges
///
/// This shape provides somewhat skewed results for a
/// sphere, especially at lower subdivisions.
/// I recommend that subdivisions of higher than 10
/// be used for acceptable results.
///
#[derive(Default, Copy, Clone, Debug)]
pub struct TetraSphereBase;

impl BaseShape for TetraSphereBase {
    #[inline]
    fn initial_points(&self) -> Vec<Vec3A> {
        consts::tetrasphere::INITIAL_POINTS.to_vec()
    }

    #[inline]
    fn triangles(&self) -> Box<[Triangle]> {
        Box::new(consts::tetrasphere::TRIANGLES)
    }
    const EDGES: usize = consts::tetrasphere::EDGES;

    #[inline]
    fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        interpolation::geometric_slerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
        interpolation::geometric_slerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        interpolation::geometric_slerp_multiple(a, b, indices, points);
    }
}

#[cfg(feature = "shape-extras")]
impl EquilateralBaseShape for TetraSphereBase {
    #[inline]
    fn triangle_normals() -> &'static [Vec3A] {
        &consts::tetrasphere::TRIANGLE_NORMALS
    }

    #[inline]
    fn triangle_min_dot() -> f32 {
        consts::tetrasphere::MIN_NORMAL_DOT
    }
}

///
/// Tetrasphere (Sphere from tetrahedron).
///
/// See [`TetraSphereBase`].
///
pub type TetraSphere<T> = Subdivided<T, TetraSphereBase>;

///
/// Implements a single triangle as the base shape.
///
/// - 3 vertices
/// - 1 face
/// - 3 edges
///
/// This is a triangle on the XZ plane. The circumscribed
/// circle on the triangle has radius 1.0.
///
#[derive(Default, Copy, Clone, Debug)]
pub struct TriangleBase;

impl BaseShape for TriangleBase {
    #[inline]
    fn initial_points(&self) -> Vec<Vec3A> {
        consts::triangle::INITIAL_POINTS.to_vec()
    }

    #[inline]
    fn triangles(&self) -> Box<[Triangle]> {
        core::slice::from_ref(&consts::triangle::TRIANGLE)
            .to_vec()
            .into()
    }
    const EDGES: usize = consts::triangle::EDGES;

    #[inline]
    fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        interpolation::lerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
        interpolation::lerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        interpolation::lerp_multiple(a, b, indices, points);
    }
}

#[cfg(feature = "shape-extras")]
impl EquilateralBaseShape for TriangleBase {
    #[inline]
    fn triangle_normals() -> &'static [Vec3A] {
        core::slice::from_ref(&consts::triangle::TRIANGLE_NORMAL)
    }

    #[inline]
    fn triangle_min_dot() -> f32 {
        -1.0
    }
}

///
/// A triangle.
///
/// See [`TriangleBase`].
///
pub type TrianglePlane<T> = Subdivided<T, TriangleBase>;

///
/// Implements a square as the base shape.
///
/// - 4 vertices
/// - 2 faces
/// - 5 edges
///
/// This is a square on the XZ plane.
///
#[derive(Default, Copy, Clone, Debug)]
pub struct SquareBase;

impl BaseShape for SquareBase {
    #[inline]
    fn initial_points(&self) -> Vec<Vec3A> {
        consts::square::INITIAL_POINTS.to_vec()
    }

    #[inline]
    fn triangles(&self) -> Box<[Triangle]> {
        consts::square::TRIANGLES.to_vec().into()
    }
    const EDGES: usize = consts::square::EDGES;

    #[inline]
    fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        interpolation::lerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
        interpolation::lerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        interpolation::lerp_multiple(a, b, indices, points);
    }
}

///
/// A square.
///
/// See [`SquareBase`].
///
pub type SquarePlane<T> = Subdivided<T, SquareBase>;

///
/// Implements a cube as the base shape.
///
/// - 8 vertices
/// - 12 faces (2 triangles per face makes 12 technically)
/// - 18 edges
///
/// This is a cube where half the diagonal is 1.0. This is to
/// enable this to be used in making a sphere.
///
#[derive(Default, Copy, Clone, Debug)]
pub struct CubeBase;

impl BaseShape for CubeBase {
    #[inline]
    fn initial_points(&self) -> Vec<Vec3A> {
        consts::cube::INITIAL_POINTS.to_vec()
    }

    #[inline]
    fn triangles(&self) -> Box<[Triangle]> {
        consts::cube::TRIANGLES.to_vec().into()
    }
    const EDGES: usize = consts::cube::EDGES;

    #[inline]
    fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        interpolation::geometric_slerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
        interpolation::geometric_slerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        interpolation::geometric_slerp_multiple(a, b, indices, points);
    }
}

///
/// A cube sphere.
///
/// See [`CubeBase`].
///
pub type CubeSphere<T> = Subdivided<T, CubeBase>;

///
/// Constant values for the shapes provided by this library.
///
mod consts {
    pub mod square {
        use crate::{Triangle, TriangleContents};
        use glam::Vec3A;

        #[rustfmt::skip]
        pub(crate) const INITIAL_POINTS: [Vec3A; 4] = [
            Vec3A::new( 1.0, 0.0,  1.0),
            Vec3A::new( 1.0, 0.0, -1.0),
            Vec3A::new(-1.0, 0.0, -1.0),
            Vec3A::new(-1.0, 0.0,  1.0),
        ];

        pub const TRIANGLES: [Triangle; 2] = [
            Triangle {
                a: 0,
                b: 1,
                c: 2,

                ab_edge: 1,
                bc_edge: 2,
                ca_edge: 0,
                ab_forward: true,
                bc_forward: true,
                ca_forward: true,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 0,
                b: 2,
                c: 3,

                ab_edge: 0,
                bc_edge: 3,
                ca_edge: 4,
                ab_forward: true,
                bc_forward: true,
                ca_forward: true,
                contents: TriangleContents::None,
            },
        ];

        pub const EDGES: usize = 5;
    }
    pub mod triangle {
        use crate::{Triangle, TriangleContents};
        use constgebra::const_soft_float::soft_f32::SoftF32;
        use glam::Vec3A;

        pub(crate) static INITIAL_POINTS: [Vec3A; 3] = [
            Vec3A::new(
                -SoftF32(3.0f32).sqrt().div(SoftF32(2.0)).to_f32(),
                0.0,
                -0.5,
            ),
            Vec3A::new(SoftF32(3.0f32).sqrt().div(SoftF32(2.0)).to_f32(), 0.0, -0.5),
            Vec3A::new(0.0, 0.0, 1.0),
        ];

        #[cfg(feature = "shape-extras")]
        pub(crate) const TRIANGLE_NORMAL: Vec3A = Vec3A::Y;

        pub const TRIANGLE: Triangle = Triangle {
            a: 2,
            b: 1,
            c: 0,

            ab_edge: 0,
            bc_edge: 1,
            ca_edge: 2,
            ab_forward: true,
            bc_forward: true,
            ca_forward: true,
            contents: TriangleContents::None,
        };

        pub const EDGES: usize = 3;
    }
    pub mod tetrasphere {
        use crate::{Triangle, TriangleContents};
        use constgebra::const_soft_float::soft_f32::SoftF32;
        #[cfg(feature = "shape-extras")]
        use constgebra::{const_soft_float::soft_f64::SoftF64, CVector, Operation};
        use glam::Vec3A;

        pub const TRIANGLES: [Triangle; 4] = [
            Triangle {
                a: 0,
                b: 1,
                c: 2,

                ab_edge: 0,
                bc_edge: 1,
                ca_edge: 2,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 0,
                b: 3,
                c: 1,

                ab_edge: 3,
                bc_edge: 4,
                ca_edge: 0,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 1,
                b: 3,
                c: 2,

                ab_edge: 4,
                bc_edge: 5,
                ca_edge: 1,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 0,
                b: 2,
                c: 3,

                ab_edge: 2,
                bc_edge: 5,
                ca_edge: 3,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
        ];
        pub const EDGES: usize = 6;

        #[cfg(feature = "shape-extras")]
        pub(super) const fn normal<const I: usize>(
            triangles: &[Triangle],
            initial_points: &[Vec3A],
        ) -> Vec3A {
            const fn f32_arr_to_f64<const N: usize>(arr: [f32; N]) -> [f64; N] {
                let mut ret = [0.0_f64; N];
                let mut i = 0;
                while i < N {
                    ret[i] = arr[i] as f64;
                    i += 1;
                }
                ret
            }
            const fn f64_arr_to_f32<const N: usize>(arr: [f64; N]) -> [f32; N] {
                let mut ret = [0.0_f32; N];
                let mut i = 0;
                while i < N {
                    ret[i] = arr[i] as f32;
                    i += 1;
                }
                ret
            }

            let triangle = &triangles[I];

            let p1 = CVector::new_vector(f32_arr_to_f64(
                initial_points[triangle.a as usize].to_array(),
            ));
            let p2 = CVector::new_vector(f32_arr_to_f64(
                initial_points[triangle.b as usize].to_array(),
            ));
            let p3 = CVector::new_vector(f32_arr_to_f64(
                initial_points[triangle.c as usize].to_array(),
            ));

            let normal = p1.add(p2).add(p3);
            let len = {
                let [x, y, z] = normal.finish_vector();
                let x = SoftF64(x);
                let y = SoftF64(y);
                let z = SoftF64(z);

                x.mul(x).add(y.mul(y)).add(z.mul(z)).sqrt()
            };
            let normal = normal.apply_each(Operation::Div(len.0));

            return Vec3A::from_array(f64_arr_to_f32(normal.finish_vector()));
        }

        pub(crate) const INITIAL_POINTS: [Vec3A; 4] = [
            Vec3A::new(
                SoftF32(8.0f32).sqrt().div(SoftF32(3.0)).to_f32(),
                SoftF32(-1.0).div(SoftF32(3.0)).to_f32(),
                0.0,
            ),
            Vec3A::new(
                -SoftF32(2.0f32).sqrt().div(SoftF32(3.0)).to_f32(),
                SoftF32(-1.0).div(SoftF32(3.0)).to_f32(),
                SoftF32(2.0f32).div(SoftF32(3.0)).sqrt().to_f32(),
            ),
            Vec3A::new(
                -SoftF32(2.0f32).sqrt().div(SoftF32(3.0)).to_f32(),
                SoftF32(-1.0).div(SoftF32(3.0)).to_f32(),
                SoftF32(2.0f32).div(SoftF32(3.0)).sqrt().neg().to_f32(),
            ),
            Vec3A::new(0.0, 1.0, 0.0),
        ];

        #[cfg(feature = "shape-extras")]
        pub(crate) const MIN_NORMAL_DOT: f32 = SoftF32(7.0f32).sqrt().div(SoftF32(3.0)).to_f32();

        #[cfg(feature = "shape-extras")]
        pub(crate) const TRIANGLE_NORMALS: [Vec3A; 4] = [
            normal::<0>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<2>(&TRIANGLES, &INITIAL_POINTS),
            normal::<3>(&TRIANGLES, &INITIAL_POINTS),
        ];
    }
    pub mod cube {
        use crate::{Triangle, TriangleContents};
        use constgebra::const_soft_float::soft_f32::SoftF32;
        use glam::Vec3A;

        #[rustfmt::skip]
        pub(crate) static INITIAL_POINTS: [Vec3A; 8] = {
            let val = SoftF32(1.0).div(SoftF32(3.0f32).sqrt()).to_f32();
            [
                Vec3A::new(-val, -val, -val),
                Vec3A::new( val, -val, -val),
                Vec3A::new(-val,  val, -val),
                Vec3A::new( val,  val, -val),

                Vec3A::new(-val, -val,  val),
                Vec3A::new( val, -val,  val),
                Vec3A::new(-val,  val,  val),
                Vec3A::new( val,  val,  val),
            ]
        };

        pub const TRIANGLES: &[Triangle; 12] = &[
            // Back
            Triangle {
                a: 0,
                b: 2,
                c: 3,

                ab_edge: 1,
                bc_edge: 2,
                ca_edge: 12,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 0,
                b: 3,
                c: 1,

                ab_edge: 12,
                bc_edge: 3,
                ca_edge: 0,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            // Top
            Triangle {
                a: 2,
                b: 7,
                c: 3,

                ab_edge: 14,
                bc_edge: 6,
                ca_edge: 2,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 2,
                b: 6,
                c: 7,

                ab_edge: 5,
                bc_edge: 10,
                ca_edge: 14,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            // Left
            Triangle {
                a: 4,
                b: 2,
                c: 0,

                ab_edge: 13,
                bc_edge: 1,
                ca_edge: 4,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 4,
                b: 6,
                c: 2,

                ab_edge: 9,
                bc_edge: 5,
                ca_edge: 13,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            // Bottom
            Triangle {
                a: 1,
                b: 5,
                c: 4,

                ab_edge: 7,
                bc_edge: 8,
                ca_edge: 16,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 1,
                b: 4,
                c: 0,

                ab_edge: 16,
                bc_edge: 4,
                ca_edge: 0,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            // Right
            Triangle {
                a: 1,
                b: 7,
                c: 5,

                ab_edge: 15,
                bc_edge: 11,
                ca_edge: 7,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 1,
                b: 3,
                c: 7,

                ab_edge: 3,
                bc_edge: 6,
                ca_edge: 15,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            // Front
            Triangle {
                a: 5,
                b: 6,
                c: 4,

                ab_edge: 17,
                bc_edge: 9,
                ca_edge: 8,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 5,
                b: 7,
                c: 6,

                ab_edge: 11,
                bc_edge: 10,
                ca_edge: 17,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
        ];

        pub const EDGES: usize = 18;
    }
    pub mod icosphere {
        use crate::{Triangle, TriangleContents};
        #[cfg(feature = "shape-extras")]
        use constgebra::const_soft_float::soft_f32::SoftF32;
        use glam::Vec3A;

        #[cfg(feature = "shape-extras")]
        use super::tetrasphere::normal;

        pub(crate) const INITIAL_POINTS: [Vec3A; 12] = [
            // North Pole
            Vec3A::Y,
            // Top Ring
            Vec3A::new(
                0.89442719099991585541,
                0.44721359549995792770,
                0.00000000000000000000,
            ),
            Vec3A::new(
                0.27639320225002106390,
                0.44721359549995792770,
                0.85065080835203987775,
            ),
            Vec3A::new(
                -0.72360679774997882507,
                0.44721359549995792770,
                0.52573111211913370333,
            ),
            Vec3A::new(
                -0.72360679774997904712,
                0.44721359549995792770,
                -0.52573111211913348129,
            ),
            Vec3A::new(
                0.27639320225002084186,
                0.44721359549995792770,
                -0.85065080835203998877,
            ),
            // Bottom Ring
            Vec3A::new(
                0.72360679774997871405,
                -0.44721359549995792770,
                -0.52573111211913392538,
            ),
            Vec3A::new(
                0.72360679774997904712,
                -0.44721359549995792770,
                0.52573111211913337026,
            ),
            Vec3A::new(
                -0.27639320225002073084,
                -0.44721359549995792770,
                0.85065080835203998877,
            ),
            Vec3A::new(
                -0.89442719099991585541,
                -0.44721359549995792770,
                0.00000000000000000000,
            ),
            Vec3A::new(
                -0.27639320225002139697,
                -0.44721359549995792770,
                -0.85065080835203976672,
            ),
            // South Pole
            Vec3A::NEG_Y,
        ];

        #[cfg(feature = "shape-extras")]
        pub(crate) static TRIANGLE_NORMALS: [Vec3A; 20] = [
            normal::<0>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<2>(&TRIANGLES, &INITIAL_POINTS),
            normal::<3>(&TRIANGLES, &INITIAL_POINTS),
            normal::<4>(&TRIANGLES, &INITIAL_POINTS),
            normal::<5>(&TRIANGLES, &INITIAL_POINTS),
            normal::<6>(&TRIANGLES, &INITIAL_POINTS),
            normal::<7>(&TRIANGLES, &INITIAL_POINTS),
            normal::<8>(&TRIANGLES, &INITIAL_POINTS),
            normal::<9>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
            normal::<1>(&TRIANGLES, &INITIAL_POINTS),
        ];

        #[cfg(feature = "shape-extras")]
        pub(crate) static MIN_NORMAL_DOT: f32 = ((SoftF32(1.0f32).div(SoftF32(30.0)))
            .mul(SoftF32(25.0).add(SoftF32(5.0_f32)).sqrt()))
        .sqrt()
        .to_f32();

        pub const TRIANGLES: [Triangle; 20] = [
            // Top
            Triangle {
                a: 0,
                b: 2,
                c: 1,

                ab_edge: 0,
                bc_edge: 5,
                ca_edge: 4,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //0
            Triangle {
                a: 0,
                b: 3,
                c: 2,

                ab_edge: 1,
                bc_edge: 6,
                ca_edge: 0,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //1
            Triangle {
                a: 0,
                b: 4,
                c: 3,

                ab_edge: 2,
                bc_edge: 7,
                ca_edge: 1,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //2
            Triangle {
                a: 0,
                b: 5,
                c: 4,

                ab_edge: 3,
                bc_edge: 8,
                ca_edge: 2,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //3
            Triangle {
                a: 0,
                b: 1,
                c: 5,

                ab_edge: 4,
                bc_edge: 9,
                ca_edge: 3,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //4
            // First ring
            Triangle {
                a: 5,
                b: 1,
                c: 6,

                ab_edge: 9,
                bc_edge: 10,
                ca_edge: 15,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //5
            Triangle {
                a: 1,
                b: 2,
                c: 7,

                ab_edge: 5,
                bc_edge: 11,
                ca_edge: 16,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //6
            Triangle {
                a: 2,
                b: 3,
                c: 8,

                ab_edge: 6,
                bc_edge: 12,
                ca_edge: 17,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //7
            Triangle {
                a: 3,
                b: 4,
                c: 9,

                ab_edge: 7,
                bc_edge: 13,
                ca_edge: 18,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //8
            Triangle {
                a: 4,
                b: 5,
                c: 10,

                ab_edge: 8,
                bc_edge: 14,
                ca_edge: 19,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //9
            // Second ring
            Triangle {
                a: 5,
                b: 6,
                c: 10,

                ab_edge: 15,
                bc_edge: 20,
                ca_edge: 14,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //10
            Triangle {
                a: 1,
                b: 7,
                c: 6,

                ab_edge: 16,
                bc_edge: 21,
                ca_edge: 10,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //11
            Triangle {
                a: 2,
                b: 8,
                c: 7,

                ab_edge: 17,
                bc_edge: 22,
                ca_edge: 11,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //12
            Triangle {
                a: 3,
                b: 9,
                c: 8,

                ab_edge: 18,
                bc_edge: 23,
                ca_edge: 12,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //13
            Triangle {
                a: 4,
                b: 10,
                c: 9,

                ab_edge: 19,
                bc_edge: 24,
                ca_edge: 13,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //14
            // Bottom
            Triangle {
                a: 10,
                b: 6,
                c: 11,

                ab_edge: 20,
                bc_edge: 26,
                ca_edge: 25,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //15
            Triangle {
                a: 6,
                b: 7,
                c: 11,

                ab_edge: 21,
                bc_edge: 27,
                ca_edge: 26,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //16
            Triangle {
                a: 7,
                b: 8,
                c: 11,

                ab_edge: 22,
                bc_edge: 28,
                ca_edge: 27,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //17
            Triangle {
                a: 8,
                b: 9,
                c: 11,

                ab_edge: 23,
                bc_edge: 29,
                ca_edge: 28,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //18
            Triangle {
                a: 9,
                b: 10,
                c: 11,

                ab_edge: 24,
                bc_edge: 25,
                ca_edge: 29,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            }, //19
        ];

        pub const EDGES: usize = 30;
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn assert_normal_numbers() {
        use super::consts;
        let all_source_points = [
            &consts::cube::INITIAL_POINTS[..],
            &consts::icosphere::INITIAL_POINTS[..],
            &consts::square::INITIAL_POINTS[..],
            &consts::tetrasphere::INITIAL_POINTS[..],
            &consts::triangle::INITIAL_POINTS[..],
        ];

        for point in all_source_points.iter().flat_map(|x| x.iter()) {
            assert!(point.is_finite());
        }
    }
}
