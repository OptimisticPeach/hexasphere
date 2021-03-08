use super::{
    BaseShape,
    Triangle,
    interpolation,
    Subdivided,
    Vec3,
};

#[derive(Default, Copy, Clone, Debug)]
pub struct IcoSphereBase;

impl<V: Vec3> BaseShape<V> for IcoSphereBase {
    #[inline]
    fn initial_points(&self) -> Vec<V> {
        consts::icosphere::INITIAL_POINTS.iter().map(|x| V::from_arr3(*x)).collect()
    }

    #[inline]
    fn triangles(&self) -> Box<[Triangle]> {
        Box::new(consts::icosphere::TRIANGLES)
    }
    const EDGES: usize = consts::icosphere::EDGES;

    #[inline]
    fn interpolate(&self, a: V, b: V, p: f32) -> V {
        interpolation::geometric_slerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(&self, a: V, b: V) -> V {
        interpolation::geometric_slerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(&self, a: V, b: V, indices: &[u32], points: &mut [V])  {
        interpolation::geometric_slerp_multiple(a, b, indices, points);
    }
}

pub type IcoSphere<T, V> = Subdivided<T, V, IcoSphereBase>;

mod consts {
    pub mod icosphere {
        use crate::{Triangle, TriangleContents};

        pub const INITIAL_POINTS: [[f32; 3]; 12] = [
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
