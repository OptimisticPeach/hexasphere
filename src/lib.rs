//!
//! Library for hexagonally tiling a sphere.
//!
//! Note this is *not* an explicit rewrite of `hexasphere.js`.
//!

use glam::Vec3A;
use std::ops::Index;

use Slice::*;

///
/// Contents of one of the main triangular faces.
///
/// This is *not* the entire face, it is the points which do
/// not include the exterior edges and vertices since those are
/// shared.
///
#[derive(Clone, Debug)]
enum TriangleContents {
    ///
    /// Nothing inside the triangle: subdivisions 0 and 1
    ///
    None,
    ///
    /// One point inside the triangle: subdivision 2
    ///
    One(u32),
    ///
    /// Three points inside the triangle: subdivision 3
    ///
    Three { a: u32, b: u32, c: u32 },
    ///
    /// Six points inside the triangle: subdivision 4
    ///
    Six {
        a: u32,
        b: u32,
        c: u32,
        ab: u32,
        bc: u32,
        ca: u32,
    },
    ///
    /// More than six points inside the triangle: subdivision > 4
    ///
    More {
        a: u32,
        b: u32,
        c: u32,
        // Separated into three `my_side_length` segments
        // to save on extra allocations.
        sides: Vec<u32>,
        my_side_length: u32,
        ///
        /// Contents of the inner triangle.
        ///
        // TODO: Make `Triangle` contain a `Vec<TriangleContents>`.
        contents: Box<TriangleContents>,
    },
}

impl TriangleContents {
    ///
    /// Creates a `None` variant.
    ///
    pub fn none() -> Self {
        Self::None
    }

    ///
    /// Creates a `One` by interpolating two values.
    ///
    fn one(ab: Slice<u32>, bc: Slice<u32>, points: &mut Vec<Vec3A>) -> Self {
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), 2);
        let p1 = points[ab[0] as usize];
        let p2 = points[bc[1] as usize];
        let index = points.len() as u32;
        points.push(geometric_slerp_half(p1, p2));
        TriangleContents::One(index)
    }

    ///
    /// Creates a `Three` variant from a `One` variant.
    ///
    fn three(&mut self, ab: Slice<u32>, bc: Slice<u32>, ca: Slice<u32>, points: &mut Vec<Vec3A>) {
        use TriangleContents::*;

        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert_eq!(ab.len(), 3);

        match self {
            &mut One(x) => {
                let ab = points[ab[1] as usize];
                let bc = points[bc[1] as usize];
                let ca = points[ca[1] as usize];

                let a = geometric_slerp_half(ab, ca);
                let b = geometric_slerp_half(bc, ab);
                let c = geometric_slerp_half(ca, bc);

                points.extend_from_slice(&[b, c]);
                points[x as usize] = a;
                *self = Three {
                    a: x,
                    b: points.len() as u32 - 2,
                    c: points.len() as u32 - 1,
                };
            }
            _ => panic!("Self is {:?} while it should be One", self),
        }
    }

    ///
    /// Creates a `Six` variant from a `Three` variant.
    ///
    fn six(&mut self, ab: Slice<u32>, bc: Slice<u32>, ca: Slice<u32>, points: &mut Vec<Vec3A>) {
        use TriangleContents::*;

        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert_eq!(ab.len(), 4);

        match self {
            &mut Three {
                a: a_index,
                b: b_index,
                c: c_index,
            } => {
                let aba = points[ab[1] as usize];
                let abb = points[ab[2] as usize];
                let bcb = points[bc[1] as usize];
                let bcc = points[bc[2] as usize];
                let cac = points[ca[1] as usize];
                let caa = points[ca[2] as usize];

                let a = geometric_slerp_half(aba, caa);
                let b = geometric_slerp_half(abb, bcb);
                let c = geometric_slerp_half(bcc, cac);

                let ab = geometric_slerp_half(a, b);
                let bc = geometric_slerp_half(b, c);
                let ca = geometric_slerp_half(c, a);

                points[a_index as usize] = a;
                points[b_index as usize] = b;
                points[c_index as usize] = c;
                points.extend_from_slice(&[ab, bc, ca]);

                *self = Six {
                    a: a_index,
                    b: b_index,
                    c: c_index,
                    ab: points.len() as u32 - 3,
                    bc: points.len() as u32 - 2,
                    ca: points.len() as u32 - 1,
                };
            }
            _ => panic!("Found {:?} whereas a Three was expected", self),
        }
    }

    ///
    /// Subdivides this given the surrounding points.
    ///
    pub fn subdivide(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut Vec<Vec3A>,
    ) {
        use TriangleContents::*;
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert!(ab.len() >= 2);
        match self {
            None => *self = Self::one(ab, bc, points),
            One(_) => self.three(ab, bc, ca, points),
            Three { .. } => self.six(ab, bc, ca, points),
            &mut Six {
                a,
                b,
                c,
                ab: ab_idx,
                bc: bc_idx,
                ca: ca_idx,
            } => {
                *self = More {
                    a,
                    b,
                    c,
                    sides: vec![ab_idx, bc_idx, ca_idx],
                    my_side_length: 1,
                    contents: Box::new(Self::none()),
                };
                self.subdivide(ab, bc, ca, points);
            }
            &mut More {
                a: a_idx,
                b: b_idx,
                c: c_idx,
                ref mut sides,
                ref mut contents,
                ref mut my_side_length,
            } => {
                points.extend_from_slice(&[Vec3A::zero(), Vec3A::zero(), Vec3A::zero()]);
                let len = points.len() as u32;
                sides.extend_from_slice(&[len - 3, len - 2, len - 1]);
                *my_side_length += 1;
                let side_length = *my_side_length as usize;

                let outer_len = ab.len();

                let aba = points[ab[1] as usize];
                let abb = points[ab[outer_len - 2] as usize];
                let bcb = points[bc[1] as usize];
                let bcc = points[bc[outer_len - 2] as usize];
                let cac = points[ca[1] as usize];
                let caa = points[ca[outer_len - 2] as usize];

                points[a_idx as usize] = geometric_slerp_half(aba, caa);
                points[b_idx as usize] = geometric_slerp_half(abb, bcb);
                points[c_idx as usize] = geometric_slerp_half(bcc, cac);

                let ab = &sides[0..side_length];
                let bc = &sides[side_length..side_length * 2];
                let ca = &sides[side_length * 2..];

                geometric_slerp_multiple(
                    points[a_idx as usize],
                    points[b_idx as usize],
                    ab,
                    points,
                );
                geometric_slerp_multiple(
                    points[b_idx as usize],
                    points[c_idx as usize],
                    bc,
                    points,
                );
                geometric_slerp_multiple(
                    points[c_idx as usize],
                    points[a_idx as usize],
                    ca,
                    points,
                );

                contents.subdivide(Forward(ab), Forward(bc), Forward(ca), points);
            }
        }
    }

    pub fn idx_ab(&self, idx: usize) -> u32 {
        use TriangleContents::*;
        match self {
            None => panic!("Invalid Index, len is 0, but got {}", idx),
            One(x) => {
                if idx != 0 {
                    panic!("Invalid Index, len is 1, but got {}", idx);
                } else {
                    *x
                }
            }
            Three { a, b, .. } => *[a, b][idx],
            Six { a, b, ab, .. } => *[a, ab, b][idx],
            &More {
                a,
                b,
                ref sides,
                my_side_length,
                ..
            } => match idx {
                0 => a,
                x if (1..(my_side_length as usize + 1)).contains(&x) => sides[x - 1],
                x if x == my_side_length as usize + 1 => b,
                _ => panic!(
                    "Invalid Index, len is {}, but got {}",
                    my_side_length + 2,
                    idx
                ),
            },
        }
    }
    pub fn idx_bc(&self, idx: usize) -> u32 {
        use TriangleContents::*;
        match self {
            None => panic!("Invalid Index, len is 0, but got {}", idx),
            One(x) => {
                if idx != 0 {
                    panic!("Invalid Index, len is 1, but got {}", idx);
                } else {
                    *x
                }
            }
            Three { c, b, .. } => *[b, c][idx],
            Six { b, c, bc, .. } => *[b, bc, c][idx],
            &More {
                b,
                c,
                ref sides,
                my_side_length,
                ..
            } => match idx {
                0 => b,
                x if (1..(my_side_length as usize + 1)).contains(&x) => {
                    sides[my_side_length as usize + (x - 1)]
                }
                x if x == my_side_length as usize + 1 => c,
                _ => panic!(
                    "Invalid Index, len is {}, but got {}",
                    my_side_length + 2,
                    idx
                ),
            },
        }
    }
    pub fn idx_ca(&self, idx: usize) -> u32 {
        use TriangleContents::*;
        match self {
            None => panic!("Invalid Index, len is 0, but got {}", idx),
            One(x) => {
                if idx != 0 {
                    panic!("Invalid Index, len is 1, but got {}", idx);
                } else {
                    *x
                }
            }
            Three { c, a, .. } => *[c, a][idx],
            Six { c, a, ca, .. } => *[c, ca, a][idx],
            &More {
                c,
                a,
                ref sides,
                my_side_length,
                ..
            } => match idx {
                0 => c,
                x if (1..(my_side_length as usize + 1)).contains(&x) => {
                    sides[my_side_length as usize * 2 + x - 1]
                }
                x if x == my_side_length as usize + 1 => a,
                _ => panic!(
                    "Invalid Index, len is {}, but got {}",
                    my_side_length + 2,
                    idx
                ),
            },
        }
    }

    pub fn add_indices(&self, buffer: &mut Vec<u32>) {
        use TriangleContents::*;
        match self {
            None | One(_) => {}
            &Three { a, b, c } => buffer.extend_from_slice(&[a, b, c]),
            &Six {
                a,
                b,
                c,
                ab,
                bc,
                ca,
            } => {
                buffer.extend_from_slice(&[a, ab, ca]);
                buffer.extend_from_slice(&[ab, b, bc]);
                buffer.extend_from_slice(&[bc, c, ca]);

                buffer.extend_from_slice(&[ab, bc, ca]);
            }
            &More {
                a,
                b,
                c,
                ref sides,
                my_side_length,
                ref contents,
            } => {
                let my_side_length = my_side_length as usize;
                let ab = &sides[0..my_side_length];
                let bc = &sides[my_side_length..my_side_length * 2];
                let ca = &sides[my_side_length * 2..];

                // Contents are always stored forward.
                add_indices_triangular(
                    a,
                    b,
                    c,
                    Forward(ab),
                    Forward(bc),
                    Forward(ca),
                    &**contents,
                    buffer,
                );
                contents.add_indices(buffer);
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
enum Slice<'a, T> {
    Forward(&'a [T]),
    Backward(&'a [T]),
}

impl<'a, T> Slice<'a, T> {
    fn len(&self) -> usize {
        match self {
            &Forward(x) | &Backward(x) => x.len(),
        }
    }
}

impl<'a, T> Index<usize> for Slice<'a, T> {
    type Output = <[T] as Index<usize>>::Output;

    fn index(&self, idx: usize) -> &Self::Output {
        match self {
            Forward(x) => x.index(idx),
            Backward(x) => x.index((x.len() - 1) - idx),
        }
    }
}

struct Triangle {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub ab: usize,
    pub bc: usize,
    pub ca: usize,
    pub ab_forward: bool,
    pub bc_forward: bool,
    pub ca_forward: bool,

    pub contents: TriangleContents,
}

impl Default for Triangle {
    fn default() -> Self {
        Self {
            a: 0,
            b: 0,
            c: 0,
            ab: 0,
            bc: 0,
            ca: 0,
            ab_forward: false,
            bc_forward: false,
            ca_forward: false,
            contents: TriangleContents::None,
        }
    }
}

impl Triangle {
    fn subdivide_edges<'a>(
        &'a mut self,
        edges: &mut [(Vec<u32>, bool); 30],
        points: &mut Vec<Vec3A>,
    ) -> usize {
        let mut divide = |p1: u32, p2: u32, edge_idx: usize, forward: &mut bool| {
            if !edges[edge_idx].1 {
                edges[edge_idx].0.push(points.len() as u32);
                points.push(Vec3A::zero());

                geometric_slerp_multiple(
                    points[p1 as usize],
                    points[p2 as usize],
                    &edges[edge_idx].0,
                    points,
                );
                edges[edge_idx].1 = true;
                *forward = true;
            } else {
                *forward = false;
            }
        };

        divide(self.a, self.b, self.ab, &mut self.ab_forward);
        divide(self.b, self.c, self.bc, &mut self.bc_forward);
        divide(self.c, self.a, self.ca, &mut self.ca_forward);

        edges[self.ab].0.len()
    }

    pub fn subdivide(&mut self, edges: &mut [(Vec<u32>, bool); 30], points: &mut Vec<Vec3A>) {
        let side_length = self.subdivide_edges(edges, points) + 1;

        if side_length > 2 {
            let ab = if self.ab_forward {
                Forward(&edges[self.ab].0)
            } else {
                Backward(&edges[self.ab].0)
            };
            let bc = if self.bc_forward {
                Forward(&edges[self.bc].0)
            } else {
                Backward(&edges[self.bc].0)
            };
            let ca = if self.ca_forward {
                Forward(&edges[self.ca].0)
            } else {
                Backward(&edges[self.ca].0)
            };
            self.contents.subdivide(ab, bc, ca, points);
        }
    }

    pub fn add_indices(&self, buffer: &mut Vec<u32>, edges: &[(Vec<u32>, bool); 30]) {
        let ab = if self.ab_forward {
            Forward(&edges[self.ab].0)
        } else {
            Backward(&edges[self.ab].0)
        };
        let bc = if self.bc_forward {
            Forward(&edges[self.bc].0)
        } else {
            Backward(&edges[self.bc].0)
        };
        let ca = if self.ca_forward {
            Forward(&edges[self.ca].0)
        } else {
            Backward(&edges[self.ca].0)
        };

        add_indices_triangular(self.a, self.b, self.c, ab, bc, ca, &self.contents, buffer);

        self.contents.add_indices(buffer);
    }
}

///
/// Hexagonally tiled sphere. Created with the `new` method.
///
/// Points are represented as the triangular subdivisions of
/// the faces of an icosahedron, the dual of which creates a
/// hexagonally tiled sphere with 12 pentagons. The number of
/// subdivisions is exactly equal to the number of hexagons
/// between adjacent pentagons.
///
/// Points are stored in normalized space, on the surface of a
/// unit sphere. The initial points are spherically interpolated
/// to preserve accuracy, and not distort due to projections
/// such as linear interpolation and then normalization.
///
pub struct Hexasphere<T> {
    points: Vec<Vec3A>,
    data: Vec<T>,
    triangles: [Triangle; 20],
    shared_edges: [(Vec<u32>, bool); 30],
    subdivisions: usize,
}

impl<T> Hexasphere<T> {
    ///
    /// Creates a hexagonally tiled sphere.
    ///
    /// The number of hexagons between pentagons is `subdivisions`, and
    /// `generator` permits you to create a unique value for every point
    /// on the sphere.
    ///
    pub fn new(subdivisions: usize, generator: impl FnMut(Vec3A) -> T) -> Self {
        let mut this = Self {
            points: vec![
                // North Pole
                Vec3A::new(0.0, 1.0, 0.0),
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
                    0.00000000000000032861,
                ),
                Vec3A::new(
                    -0.27639320225002139697,
                    -0.44721359549995792770,
                    -0.85065080835203976672,
                ),
                // South Pole
                Vec3A::new(0.0, -1.0, 0.0),
            ],
            shared_edges: [
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
                (Vec::new(), true),
            ],
            triangles: [
                // Top
                Triangle {
                    a: 0,
                    b: 2,
                    c: 1,

                    ab: 0,
                    bc: 5,
                    ca: 4,
                    ..Default::default()
                }, //0
                Triangle {
                    a: 0,
                    b: 3,
                    c: 2,

                    ab: 1,
                    bc: 6,
                    ca: 0,
                    ..Default::default()
                }, //1
                Triangle {
                    a: 0,
                    b: 4,
                    c: 3,

                    ab: 2,
                    bc: 7,
                    ca: 1,
                    ..Default::default()
                }, //2
                Triangle {
                    a: 0,
                    b: 5,
                    c: 4,

                    ab: 3,
                    bc: 8,
                    ca: 2,
                    ..Default::default()
                }, //3
                Triangle {
                    a: 0,
                    b: 1,
                    c: 5,

                    ab: 4,
                    bc: 9,
                    ca: 3,
                    ..Default::default()
                }, //4
                // First ring
                Triangle {
                    a: 5,
                    b: 1,
                    c: 6,

                    ab: 9,
                    bc: 10,
                    ca: 15,
                    ..Default::default()
                }, //5
                Triangle {
                    a: 1,
                    b: 2,
                    c: 7,

                    ab: 5,
                    bc: 11,
                    ca: 16,
                    ..Default::default()
                }, //6
                Triangle {
                    a: 2,
                    b: 3,
                    c: 8,

                    ab: 6,
                    bc: 12,
                    ca: 17,
                    ..Default::default()
                }, //7
                Triangle {
                    a: 3,
                    b: 4,
                    c: 9,

                    ab: 7,
                    bc: 13,
                    ca: 18,
                    ..Default::default()
                }, //8
                Triangle {
                    a: 4,
                    b: 5,
                    c: 10,

                    ab: 8,
                    bc: 14,
                    ca: 19,
                    ..Default::default()
                }, //9
                // Second ring
                Triangle {
                    a: 5,
                    b: 6,
                    c: 10,

                    ab: 15,
                    bc: 20,
                    ca: 14,
                    ..Default::default()
                }, //10
                Triangle {
                    a: 1,
                    b: 7,
                    c: 6,

                    ab: 16,
                    bc: 21,
                    ca: 10,
                    ..Default::default()
                }, //11
                Triangle {
                    a: 2,
                    b: 8,
                    c: 7,

                    ab: 17,
                    bc: 22,
                    ca: 11,
                    ..Default::default()
                }, //12
                Triangle {
                    a: 3,
                    b: 9,
                    c: 8,

                    ab: 18,
                    bc: 23,
                    ca: 12,
                    ..Default::default()
                }, //13
                Triangle {
                    a: 4,
                    b: 10,
                    c: 9,

                    ab: 19,
                    bc: 24,
                    ca: 13,
                    ..Default::default()
                }, //14
                // Bottom
                Triangle {
                    a: 10,
                    b: 6,
                    c: 11,

                    ab: 20,
                    bc: 26,
                    ca: 25,
                    ..Default::default()
                }, //15
                Triangle {
                    a: 6,
                    b: 7,
                    c: 11,

                    ab: 21,
                    bc: 27,
                    ca: 26,
                    ..Default::default()
                }, //16
                Triangle {
                    a: 7,
                    b: 8,
                    c: 11,

                    ab: 22,
                    bc: 28,
                    ca: 27,
                    ..Default::default()
                }, //17
                Triangle {
                    a: 8,
                    b: 9,
                    c: 11,

                    ab: 23,
                    bc: 29,
                    ca: 28,
                    ..Default::default()
                }, //18
                Triangle {
                    a: 9,
                    b: 10,
                    c: 11,

                    ab: 24,
                    bc: 25,
                    ca: 29,
                    ..Default::default()
                }, //19
            ],
            subdivisions: 1,
            data: vec![],
        };

        for _ in 0..subdivisions {
            this.subdivide();
        }

        this.data = this.points.iter().copied().map(generator).collect();

        this
    }

    fn subdivide(&mut self) {
        for (_, done) in &mut self.shared_edges {
            *done = false;
        }

        for triangle in &mut self.triangles {
            triangle.subdivide(&mut self.shared_edges, &mut self.points);
        }
    }

    pub fn raw_points(&self) -> &[Vec3A] {
        &self.points
    }

    ///
    /// Appends the indices for the triangle into `buffer`.
    ///
    pub fn get_indices(&self, triangle: usize, buffer: &mut Vec<u32>) {
        self.triangles[triangle].add_indices(buffer, &self.shared_edges);
    }

    pub fn subdivisions(&self) -> usize {
        self.subdivisions
    }

    pub fn raw_data(&self) -> &[T] {
        &self.data
    }

    ///
    /// Calculate the number of indices which each main
    /// triangle will add to the vertex buffer.
    ///
    /// # Equation
    ///
    /// ```text
    /// (subdivisions + 1)²
    /// ```
    ///
    pub fn indices_per_main_triangle(&self) -> usize {
        (self.subdivisions + 1) * (self.subdivisions + 1)
    }

    ///
    /// Calculate the number of vertices contained within
    /// each main triangle including the vertices which are
    /// shared with another main triangle.
    ///
    /// # Equation
    ///
    /// ```text
    /// (subdivisions + 1) * (subdivisions + 2) / 2
    /// ```
    ///
    pub fn vertices_per_main_triangle_shared(&self) -> usize {
        (self.subdivisions + 1) * (self.subdivisions + 2) / 2
    }

    ///
    /// Calculate the number of vertices contained within each
    /// main triangle excluding the ones that are shared with
    /// other main triangles.
    ///
    /// # Equation
    ///
    /// ```text
    /// {
    /// { subdivisions < 2  : 0
    /// {
    /// { subdivisions >= 2 : (subdivisions - 1) * subdivisions / 2
    /// {
    /// ```
    ///
    pub fn vertices_per_main_triangle_unique(&self) -> usize {
        if self.subdivisions < 2 {
            return 0;
        }
        (self.subdivisions - 1) * self.subdivisions / 2
    }

    ///
    /// Calculate the number of vertices along the edges
    /// of the main triangles and the vertices of the main
    /// triangles.
    ///
    /// # Equation
    ///
    /// ```text
    /// subdivisions * 30 + 12
    /// ```
    ///
    pub fn shared_vertices(&self) -> usize {
        self.subdivisions * 30 + 12
    }
}

/// Note: `a` and `b` should both be normalized for normalized results.
#[allow(dead_code)]
fn geometric_slerp(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
    let angle = a.dot(b).acos();

    let sin = angle.sin().recip();
    a * (((1.0 - p) * angle).sin() * sin) + b * ((p * angle).sin() * sin)
}

fn geometric_slerp_half(a: Vec3A, b: Vec3A) -> Vec3A {
    let angle = a.dot(b).acos();

    let sin_denom = angle.sin().recip();
    let sin_numer = (angle * 0.5).sin();

    (a + b) * sin_denom * sin_numer
}

/// Note: `a` and `b` should both be normalized for normalized results.
fn geometric_slerp_multiple<'a>(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
    let angle = a.dot(b).acos();
    let sin = angle.sin().recip();

    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] =
            a * (((1.0 - percent) * angle).sin() * sin) + b * ((percent * angle).sin() * sin);
    }
}

fn add_indices_triangular(
    a: u32,
    b: u32,
    c: u32,
    ab: Slice<u32>,
    bc: Slice<u32>,
    ca: Slice<u32>,
    contents: &TriangleContents,
    buffer: &mut Vec<u32>,
) {
    let subdivisions = ab.len();
    if subdivisions == 0 {
        buffer.extend_from_slice(&[a, b, c]);
        return;
    } else if subdivisions == 1 {
        buffer.extend_from_slice(&[a, ab[0], ca[0]]);
        buffer.extend_from_slice(&[b, bc[0], ab[0]]);
        buffer.extend_from_slice(&[c, ca[0], bc[0]]);
        buffer.extend_from_slice(&[ab[0], bc[0], ca[0]]);
        return;
    } else if subdivisions == 2 {
        buffer.extend_from_slice(&[a, ab[0], ca[1]]);
        buffer.extend_from_slice(&[b, bc[0], ab[1]]);
        buffer.extend_from_slice(&[c, ca[0], bc[1]]);

        buffer.extend_from_slice(&[ab[1], contents.idx_ab(0), ab[0]]);
        buffer.extend_from_slice(&[bc[1], contents.idx_ab(0), bc[0]]);
        buffer.extend_from_slice(&[ca[1], contents.idx_ab(0), ca[0]]);

        buffer.extend_from_slice(&[ab[0], contents.idx_ab(0), ca[1]]);
        buffer.extend_from_slice(&[bc[0], contents.idx_ab(0), ab[1]]);
        buffer.extend_from_slice(&[ca[0], contents.idx_ab(0), bc[1]]);
        return;
    }

    let last_idx = ab.len() - 1;

    buffer.extend_from_slice(&[a, ab[0], ca[last_idx]]);
    buffer.extend_from_slice(&[b, bc[0], ab[last_idx]]);
    buffer.extend_from_slice(&[c, ca[0], bc[last_idx]]);

    buffer.extend_from_slice(&[ab[0], contents.idx_ab(0), ca[last_idx]]);
    buffer.extend_from_slice(&[bc[0], contents.idx_bc(0), ab[last_idx]]);
    buffer.extend_from_slice(&[ca[0], contents.idx_ca(0), bc[last_idx]]);

    for i in 0..last_idx - 1 {
        // Exclude special case: last_idx - 1.
        // AB
        buffer.extend_from_slice(&[ab[i], ab[i + 1], contents.idx_ab(i)]);
        buffer.extend_from_slice(&[ab[i + 1], contents.idx_ab(i + 1), contents.idx_ab(i)]);
        // BC
        buffer.extend_from_slice(&[bc[i], bc[i + 1], contents.idx_bc(i)]);
        buffer.extend_from_slice(&[bc[i + 1], contents.idx_bc(i + 1), contents.idx_bc(i)]);
        // CA
        buffer.extend_from_slice(&[ca[i], ca[i + 1], contents.idx_ca(i)]);
        buffer.extend_from_slice(&[ca[i + 1], contents.idx_ca(i + 1), contents.idx_ca(i)]);
    }

    // Deal with special case: last_idx - 1
    buffer.extend_from_slice(&[
        ab[last_idx],
        contents.idx_ab(last_idx - 1),
        ab[last_idx - 1],
    ]);

    buffer.extend_from_slice(&[
        bc[last_idx],
        contents.idx_bc(last_idx - 1),
        bc[last_idx - 1],
    ]);

    buffer.extend_from_slice(&[
        ca[last_idx],
        contents.idx_ca(last_idx - 1),
        ca[last_idx - 1],
    ]);
}

#[cfg(test)]
mod tests {
    use crate::Hexasphere;
    use crate::Slice::Forward;
    use glam::Vec3A;

    // Starting points aren't _quite_ precise enough to use `f32::EPSILON`.
    const EPSILON: f32 = 0.0000002;

    #[test]
    fn slerp_one() {
        use super::geometric_slerp_half;
        let p1 = Vec3A::new(0.360492952832, 0.932761936915, 0.0);
        let p2 = Vec3A::new(0.975897449331, 0.218229623081, 0.0);

        let expected = Vec3A::new(0.757709663147, 0.652591806854, 0.0);

        let result = geometric_slerp_half(p1, p2);

        assert!((expected - result).length() <= EPSILON);

        // Another test case
        let p1 = Vec3A::new(-0.24953852315, 0.0, 0.968364872073);
        let p2 = Vec3A::new(-0.948416666565, 0.0, 0.317026539239);

        let expected = Vec3A::new(-0.681787771301, 0.0, 0.731550022148);

        let result = geometric_slerp_half(p1, p2);

        assert!((expected - result).length() <= EPSILON);
    }

    #[test]
    fn slerp_many() {
        use super::geometric_slerp_multiple;

        let p1 = Vec3A::new(0.0, -0.885330189449, 0.464962854054);
        let p2 = Vec3A::new(0.0, 0.946042343528, 0.324043028395);

        let expected = Vec3A::new(0.0, 0.0767208624118, 0.997052611085);

        let mut result = Vec3A::zero();
        geometric_slerp_multiple(p1, p2, &[0], std::slice::from_mut(&mut result));

        assert!((expected - result).length() <= EPSILON);

        let p1 = Vec3A::new(0.876621956288, 0.0, 0.481179743707);
        let p2 = Vec3A::new(-0.391617625614, 0.0, -0.920128053756);

        let expected = [
            Vec3A::new(0.999975758841, 0.0, 0.00696288230076),
            Vec3A::new(0.883237589397, 0.0, -0.468925751774),
            Vec3A::new(0.554436024709, 0.0, -0.83222634812),
            Vec3A::new(0.0925155945469, 0.0, -0.995711235633),
        ];

        let mut result = [Vec3A::zero(), Vec3A::zero(), Vec3A::zero(), Vec3A::zero()];

        geometric_slerp_multiple(p1, p2, &[0, 1, 2, 3], &mut result);

        for (&expected, &result) in expected.iter().zip(result.iter()) {
            assert!((expected - result).length() <= EPSILON);
        }
    }

    #[test]
    fn new() {
        let x = Hexasphere::new(0, |_| ());
        x.get_indices(0, &mut Vec::new());
    }

    #[test]
    fn one() {
        let x = Hexasphere::new(1, |_| ());
        x.get_indices(0, &mut Vec::new());
    }

    #[test]
    fn second_layer_inner() {
        let x = Hexasphere::new(2, |_| ());
        x.get_indices(0, &mut Vec::new());
        let x = Hexasphere::new(3, |_| ());
        x.get_indices(0, &mut Vec::new());
        let x = Hexasphere::new(5, |_| ());
        x.get_indices(0, &mut Vec::new());
        let x = Hexasphere::new(6, |_| ());
        x.get_indices(0, &mut Vec::new());
    }

    #[test]
    fn indices_zero() {
        use super::add_indices_triangular;
        use super::TriangleContents;

        let mut buffer = Vec::new();

        add_indices_triangular(
            0,
            1,
            2,
            Forward(&[]),
            Forward(&[]),
            Forward(&[]),
            &TriangleContents::none(),
            &mut buffer,
        );

        assert_eq!(buffer, &[0, 1, 2]);
    }

    #[test]
    fn indices_one() {
        use super::add_indices_triangular;
        use super::TriangleContents;

        let mut buffer = Vec::new();

        add_indices_triangular(
            0,
            1,
            2,
            Forward(&[3]),
            Forward(&[4]),
            Forward(&[5]),
            &TriangleContents::none(),
            &mut buffer,
        );

        assert_eq!(buffer, &[0, 3, 5, 1, 4, 3, 2, 5, 4, 3, 4, 5,]);
    }

    #[test]
    fn indices_two() {
        use super::add_indices_triangular;
        use super::TriangleContents;

        let mut buffer = Vec::new();

        add_indices_triangular(
            0,
            3,
            6,
            Forward(&[1, 2]),
            Forward(&[4, 5]),
            Forward(&[7, 8]),
            &TriangleContents::One(9),
            &mut buffer,
        );

        assert_eq!(
            buffer,
            &[0, 1, 8, 3, 4, 2, 6, 7, 5, 2, 9, 1, 5, 9, 4, 8, 9, 7, 1, 9, 8, 4, 9, 2, 7, 9, 5,]
        );
    }

    // Really, we're testing for the rest.
    #[test]
    fn indices_three() {
        use super::add_indices_triangular;
        use super::TriangleContents;

        let mut buffer = Vec::new();

        add_indices_triangular(
            0,
            4,
            8,
            Forward(&[1, 2, 3]),
            Forward(&[5, 6, 7]),
            Forward(&[9, 10, 11]),
            &TriangleContents::Three {
                a: 12,
                b: 13,
                c: 14,
            },
            &mut buffer,
        );

        assert_eq!(
            buffer,
            &[
                0, 1, 11, 4, 5, 3, 8, 9, 7, 1, 12, 11, 5, 13, 3, 9, 14, 7, 1, 2, 12, 2, 13, 12, 5,
                6, 13, 6, 14, 13, 9, 10, 14, 10, 12, 14, 3, 13, 2, 7, 14, 6, 11, 12, 10,
            ][..]
        );
    }

    #[test]
    fn precision() {
        let sphere = Hexasphere::new(10, |_| ());

        for i in sphere.raw_points() {
            assert!(i.length() - 1.0 <= EPSILON);
        }
    }
}
