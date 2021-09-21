use super::{
    BaseShape,
    Triangle,
    interpolation,
    EquilateralBaseShape,
    Subdivided,

};
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
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A])  {
        interpolation::geometric_slerp_multiple(a, b, indices, points);
    }
}

impl EquilateralBaseShape for IcoSphereBase {
    #[inline]
    fn triangle_normals() -> &'static [Vec3A] {
        &*consts::icosphere::TRIANGLE_NORMALS
    }

    #[inline]
    fn triangle_min_dot() -> f32 {
        *consts::icosphere::MIN_NORMAL_DOT
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

        (angle * 0.5).sin() * 2.0
    }
}

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

impl EquilateralBaseShape for TetraSphereBase {
    #[inline]
    fn triangle_normals() -> &'static [Vec3A] {
        &*consts::tetrasphere::TRIANGLE_NORMALS
    }

    #[inline]
    fn triangle_min_dot() -> f32 {
        *consts::tetrasphere::MIN_NORMAL_DOT
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
        core::slice::from_ref(&consts::triangle::TRIANGLE).to_vec().into()
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

impl EquilateralBaseShape for TriangleBase {
    #[inline]
    fn triangle_normals() -> &'static [Vec3A] {
        core::slice::from_ref(&*consts::triangle::TRIANGLE_NORMAL)
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

        lazy_static::lazy_static! {
            pub(crate) static ref INITIAL_POINTS: [Vec3A; 4] = [
                Vec3A::new(1.0, 0.0, 1.0),
                Vec3A::new(1.0, 0.0, -1.0),
                Vec3A::new(-1.0, 0.0, -1.0),
                Vec3A::new(-1.0, 0.0, 1.0),
            ];

            pub(crate) static ref TRIANGLE_NORMALS: [Vec3A; 2] = [
                Vec3A::new(0.0, 1.0, 0.0),
                Vec3A::new(0.0, 1.0, 0.0),
            ];
        }

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
        use glam::Vec3A;

        lazy_static::lazy_static! {
            pub(crate) static ref INITIAL_POINTS: [Vec3A; 3] = [
                Vec3A::new(-3.0f32.sqrt() / 2.0, 0.0, -0.5),
                Vec3A::new( 3.0f32.sqrt() / 2.0, 0.0, -0.5),
                Vec3A::new(0.0, 0.0, 1.0),
            ];

            pub(crate) static ref TRIANGLE_NORMAL: Vec3A = Vec3A::new(0.0, 1.0, 0.0);
        }

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
        use glam::Vec3A;

        lazy_static::lazy_static! {
            pub(crate) static ref INITIAL_POINTS: [Vec3A; 4] = [
                Vec3A::new((8.0f32).sqrt() / 3.0, -1.0 / 3.0, 0.0),
                Vec3A::new(-(2.0f32).sqrt() / 3.0, -1.0 / 3.0, (2.0f32 / 3.0).sqrt()),
                Vec3A::new(-(2.0f32).sqrt() / 3.0, -1.0 / 3.0, -(2.0f32 / 3.0).sqrt()),
                Vec3A::new(0.0, 1.0, 0.0),
            ];

            pub(crate) static ref TRIANGLE_NORMALS: [Vec3A; 4] = {
                fn normal(triangle: usize) -> Vec3A {
                    (
                        INITIAL_POINTS[TRIANGLES[triangle].a as usize] +
                        INITIAL_POINTS[TRIANGLES[triangle].b as usize] +
                        INITIAL_POINTS[TRIANGLES[triangle].c as usize]
                    ) / 3.0
                }

                [
                    normal(0),
                    normal(1),
                    normal(2),
                    normal(3),
                ]
            };

            pub(crate) static ref MIN_NORMAL_DOT: f32 = 7.0f32.sqrt() / 3.0;
        }

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

            }
        ];
        pub const EDGES: usize = 6;
    }
    pub mod cube {
        use crate::{Triangle, TriangleContents};
        use glam::Vec3A;

        lazy_static::lazy_static! {
            pub(crate) static ref INITIAL_POINTS: [Vec3A; 8] = {
                let val = (3.0f32).sqrt().recip();
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
        }

        pub const TRIANGLES: [Triangle; 12] = [
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
        use glam::Vec3A;

        lazy_static::lazy_static! {
            pub(crate) static ref INITIAL_POINTS: [Vec3A; 12] = [
                Vec3A::from(RAW_POINTS[0]),
                Vec3A::from(RAW_POINTS[1]),
                Vec3A::from(RAW_POINTS[2]),
                Vec3A::from(RAW_POINTS[3]),
                Vec3A::from(RAW_POINTS[4]),
                Vec3A::from(RAW_POINTS[5]),
                Vec3A::from(RAW_POINTS[6]),
                Vec3A::from(RAW_POINTS[7]),
                Vec3A::from(RAW_POINTS[8]),
                Vec3A::from(RAW_POINTS[9]),
                Vec3A::from(RAW_POINTS[10]),
                Vec3A::from(RAW_POINTS[11]),
            ];

            pub(crate) static ref TRIANGLE_NORMALS: [Vec3A; 20] = {
                fn normal(triangle: usize) -> Vec3A {
                    (
                        INITIAL_POINTS[TRIANGLES[triangle].a as usize] +
                        INITIAL_POINTS[TRIANGLES[triangle].b as usize] +
                        INITIAL_POINTS[TRIANGLES[triangle].c as usize]
                    ) / 3.0
                }

                [
                    normal(0),
                    normal(1),
                    normal(2),
                    normal(3),
                    normal(4),
                    normal(5),
                    normal(6),
                    normal(7),
                    normal(8),
                    normal(9),
                    normal(10),
                    normal(11),
                    normal(12),
                    normal(13),
                    normal(14),
                    normal(15),
                    normal(16),
                    normal(17),
                    normal(18),
                    normal(19),
                ]
            };

            pub(crate) static ref MIN_NORMAL_DOT: f32 = ((1.0f32 / 30.0) * (25.0 + 5.0_f32.sqrt())).sqrt();
        }

        const RAW_POINTS: [[f32; 3]; 12] = [
            // North Pole
            [0.0, 1.0, 0.0],
            // Top Ring
            [
                0.89442719099991585541,
                0.44721359549995792770,
                0.00000000000000000000,
            ],
            [
                0.27639320225002106390,
                0.44721359549995792770,
                0.85065080835203987775,
            ],
            [
                -0.72360679774997882507,
                0.44721359549995792770,
                0.52573111211913370333,
            ],
            [
                -0.72360679774997904712,
                0.44721359549995792770,
                -0.52573111211913348129,
            ],
            [
                0.27639320225002084186,
                0.44721359549995792770,
                -0.85065080835203998877,
            ],
            // Bottom Ring
            [
                0.72360679774997871405,
                -0.44721359549995792770,
                -0.52573111211913392538,
            ],
            [
                0.72360679774997904712,
                -0.44721359549995792770,
                0.52573111211913337026,
            ],
            [
                -0.27639320225002073084,
                -0.44721359549995792770,
                0.85065080835203998877,
            ],
            [
                -0.89442719099991585541,
                -0.44721359549995792770,
                0.00000000000000000000,
            ],
            [
                -0.27639320225002139697,
                -0.44721359549995792770,
                -0.85065080835203976672,
            ],
            // South Pole
            [0.0, -1.0, 0.0],
        ];

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
