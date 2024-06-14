//!
//! Library for subdividing shapes made of triangles.
//!
//! This library defines `Subdivided<T, S>`. This struct
//! allows one to define a base shape using `S` and the
//! `BaseShape` trait, and to subdivide it using the
//! interpolation functions defined as part of `S`.
//!
//! This includes a few base shapes:
//!
//! - Icosahedron
//! - Tetrahedron
//! - Square
//! - Triangle
//! - Cube
//!
//! ## Example usage
//!
//! ```rust
//! use hexasphere::shapes::IcoSphere;
//!
//! fn main() {
//!     // Create a new sphere with 20 subdivisions
//!     // an no data associated with the vertices.
//!     let sphere = IcoSphere::new(20, |_| ());
//!
//!     let points = sphere.raw_points();
//!     for p in points {
//!         println!("{:?} is a point on the sphere!", p);
//!     }
//!     let indices = sphere.get_all_indices();
//!     for triangle in indices.chunks(3) {
//!         println!(
//!             "[{}, {}, {}] is a triangle on the resulting shape",
//!             triangle[0],
//!             triangle[1],
//!             triangle[2],
//!         );
//!     }
//! }
//! ```
//!
//! # Features
//! - `adjacency` allows the user to create neighbour maps from
//! the indices provided by the `Subdivided` struct.
//!

use glam::Vec3A;

use slice::Slice::*;
use slice::*;

#[cfg(feature = "adjacency")]
pub use adjacency::*;

pub mod interpolation;
pub mod shapes;
mod slice;

///
/// Defines the setup for a base shape, and the functions
/// used in interpolation.
///
/// If you want to use a different interpolation function,
/// implement this trait for a new type, and carry over only
/// the properties you'd like to keep:
///
/// ```rust
/// # use hexasphere::BaseShape;
/// use hexasphere::{Triangle, shapes::IcoSphereBase};
/// use glam::Vec3A;
/// // Uses linear interpolation instead of spherical.
/// struct FlatIcosahedron;
///
/// impl BaseShape for FlatIcosahedron {
///     // Keep the initial parameters.
///     fn initial_points(&self) -> Vec<Vec3A> {
///         IcoSphereBase.initial_points()
///     }
///
///     fn triangles(&self) -> Box<[Triangle]> {
///         IcoSphereBase.triangles()
///     }
///     const EDGES: usize = IcoSphereBase::EDGES;
///
///     // Swap out what you'd like to change.
///     fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
///         hexasphere::interpolation::lerp(a, b, p)
///     }
///
///     fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
///         hexasphere::interpolation::lerp_half(a, b)
///     }
///
///     fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
///         hexasphere::interpolation::lerp_multiple(a, b, indices, points);
///     }
/// }
/// ```
///
/// Or, create your own shape, by changing the values for
/// [`initial_points`], [`triangles`], and [`EDGES`]. Check
/// the documentation for these members individually on how
/// they should be implemented.
///
/// [`initial_points`]: #tymethod.initial_points
/// [`triangles`]: #tymethod.triangles
/// [`EDGES`]: #associatedconstant.EDGES
///
pub trait BaseShape {
    ///
    /// The initial vertices for the triangle. Note that
    /// `Vec3A::new` is not a `const fn()`, hence I recommend
    /// you use `lazy_static`. Check the source file for this
    /// crate and look for the constants module at the bottom
    /// for an example.
    ///
    /// Constraints on the points depend on the interpolation
    /// function used:
    /// - `slerp` requires normalized (magnitude 1) data.
    /// - `lerp` doesn't care.
    /// - `normalized_lerp` requires normalized (magnitude 1)
    /// data.
    ///
    fn initial_points(&self) -> Vec<Vec3A>;

    ///
    /// Base triangles for the shape.
    ///
    /// - The fields `a`, `b`, and `c` define the indices for
    /// the points of the triangles given the data present
    /// in `initial_points`. Note that this crate assumes
    /// points are in a counter clockwise ordering.
    /// - The fields `ab_edge`, `bc_edge`, `ca_edge` mark the
    /// index of the edge which `a`/`b`, `b`/`c`, and `c`/`a`
    /// border respectively. While theoretically you could give
    /// each triangle its own edge, minimizing edges saves on
    /// memory footprint and performance.
    /// - Triangles should be created through [`Triangle::new`].
    ///
    fn triangles(&self) -> Box<[Triangle]>;

    ///
    /// Number of unique edges defined in the contents of
    /// `triangles()`. This number is 5 for a square for
    /// example:
    /// ```text
    /// a - 0 - b
    /// |     / |
    /// 3   4   1
    /// | /     |
    /// d - 2 - c
    /// ```
    ///
    const EDGES: usize;

    ///
    /// Basic function used for interpolation. When `p` is
    /// `0.0`, `a` is expected. When `p` is `1.0`, `b` is
    /// expected. There are three options already implemented
    /// in this crate:
    /// - [`lerp`] implements linear interpolation.
    /// - [`geometric_slerp`] implements spherical
    /// interpolation. (Interpolates along an arc on a sphere)
    /// - [`normalized_lerp`] implements cheaper
    /// yet less accurate spherical interpolation. The accuracy
    /// decreases as the angle between the two points on the unit
    /// sphere increases.
    ///
    /// [`lerp`]: ../fn.lerp.html
    /// [`geometric_slerp`]: ../fn.geometric_slerp.html
    /// [`normalized_lerp`]: ../fn.normalized_lerp.html
    ///
    fn interpolate(&self, a: Vec3A, b: Vec3A, p: f32) -> Vec3A;

    ///
    /// If an optimization is available for the case where `p`
    /// is `0.5`, this function should implement it. This defaults
    /// to calling `interpolate(a, b, 0.5)` however.
    ///
    fn interpolate_half(&self, a: Vec3A, b: Vec3A) -> Vec3A {
        self.interpolate(a, b, 0.5)
    }

    ///
    /// If an optimization is available for the case where `p`
    /// varies but `a` and `b` don't this function should implement
    /// it.
    ///
    /// ### Parameters
    /// - `a`: start.
    /// - `b`: end.
    /// - `indices`: list of indices to index into `points`. `points`
    /// at the index should contain the result. The index (n) of an
    /// index should correspond with the nth point in a distribution
    /// of `indices.len()` points between (exclusive) `a` and `b`.
    /// - `points`: list of points where the results of the calculation
    /// should end up. To be indexed by values in `indices`.
    ///
    fn interpolate_multiple(&self, a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        for (percent, index) in indices.iter().enumerate() {
            let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

            points[*index as usize] = self.interpolate(a, b, percent);
        }
    }
}

///
/// Implemented in the case where the triangles on the shape
/// are both equilateral and identifiable from their normal.
///
/// This is only used in the cases of spherical shapes.
///
#[cfg(feature = "shape-extras")]
pub trait EquilateralBaseShape: BaseShape {
    ///
    /// Normals for each of the triangles provided by
    /// [`BaseShape::triangles`].
    ///
    fn triangle_normals() -> &'static [Vec3A];
    ///
    /// Minimum value for the dot product which one could use
    /// to determine that triangle being the closest.
    ///
    /// If the dot product between a vector and a triangle
    /// normal is greater than this, then that vector is
    /// guaranteed to be within that triangle.
    ///
    /// This is also the cosine of the angle of the cone
    /// whose circular edge lies on a triangle.
    ///
    fn triangle_min_dot() -> f32;
}

///
/// The edge between two main triangles.
///
#[derive(Debug)]
struct Edge {
    ///
    /// Indices of the points between the endpoints.
    ///
    points: Vec<u32>,
    ///
    /// Whether this edge has already been processed
    /// or not.
    ///
    done: bool,
}

impl Default for Edge {
    fn default() -> Self {
        Self {
            points: Vec::new(),
            done: true,
        }
    }
}

impl Edge {
    pub fn subdivide_n_times(&mut self, n: usize, points: &mut usize) {
        for _ in 0..n {
            self.points.push(*points as _);
            *points += 1;
        }
    }
}

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
        // Implementing this as a `Vec<TriangleContents>` would
        // probably be a perf. improvement someday, however not
        // something worth implementing right now.
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
    fn one(points: &mut usize) -> Self {
        let index = *points as u32;
        *points += 1;
        TriangleContents::One(index)
    }

    fn calculate_one(
        &self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        points: &mut [Vec3A],
        shape: &impl BaseShape,
    ) {
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), 2);
        match self {
            TriangleContents::One(idx) => {
                let p1 = points[ab[0] as usize];
                let p2 = points[bc[1] as usize];

                points[*idx as usize] = shape.interpolate_half(p1, p2);
            }
            _ => panic!("Did not find One variant."),
        }
    }

    ///
    /// Creates a `Three` variant from a `One` variant.
    ///
    fn three(&mut self, points: &mut usize) {
        use TriangleContents::*;

        match self {
            &mut One(x) => {
                *points += 2;

                *self = Three {
                    a: x,
                    b: *points as u32 - 2,
                    c: *points as u32 - 1,
                };
            }
            _ => panic!("Self is {:?} while it should be One", self),
        }
    }

    fn calculate_three(
        &self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut [Vec3A],
        shape: &impl BaseShape,
    ) {
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert_eq!(ab.len(), 3);

        match self {
            &TriangleContents::Three { a, b, c } => {
                let ab = points[ab[1] as usize];
                let bc = points[bc[1] as usize];
                let ca = points[ca[1] as usize];

                let a_val = shape.interpolate_half(ab, ca);
                let b_val = shape.interpolate_half(bc, ab);
                let c_val = shape.interpolate_half(ca, bc);

                points[a as usize] = a_val;
                points[b as usize] = b_val;
                points[c as usize] = c_val;
            }
            _ => panic!("Did not find Three variant."),
        }
    }

    ///
    /// Creates a `Six` variant from a `Three` variant.
    ///
    fn six(&mut self, points: &mut usize) {
        use TriangleContents::*;

        match self {
            &mut Three {
                a: a_index,
                b: b_index,
                c: c_index,
            } => {
                *points += 3;

                *self = Six {
                    a: a_index,
                    b: b_index,
                    c: c_index,
                    ab: *points as u32 - 3,
                    bc: *points as u32 - 2,
                    ca: *points as u32 - 1,
                };
            }
            _ => panic!("Found {:?} whereas a Three was expected", self),
        }
    }

    fn calculate_six(
        &self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut [Vec3A],
        shape: &impl BaseShape,
    ) {
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert_eq!(ab.len(), 4);

        use TriangleContents::*;

        match self {
            &Six {
                a: a_index,
                b: b_index,
                c: c_index,
                ab: ab_index,
                bc: bc_index,
                ca: ca_index,
            } => {
                let aba = points[ab[1] as usize];
                let abb = points[ab[2] as usize];
                let bcb = points[bc[1] as usize];
                let bcc = points[bc[2] as usize];
                let cac = points[ca[1] as usize];
                let caa = points[ca[2] as usize];

                let a = shape.interpolate_half(aba, caa);
                let b = shape.interpolate_half(abb, bcb);
                let c = shape.interpolate_half(bcc, cac);

                let ab = shape.interpolate_half(a, b);
                let bc = shape.interpolate_half(b, c);
                let ca = shape.interpolate_half(c, a);

                points[a_index as usize] = a;
                points[b_index as usize] = b;
                points[c_index as usize] = c;
                points[ab_index as usize] = ab;
                points[bc_index as usize] = bc;
                points[ca_index as usize] = ca;
            }
            _ => panic!("Found {:?} whereas a Three was expected", self),
        }
    }

    ///
    /// Subdivides this given the surrounding points.
    ///
    pub fn subdivide(&mut self, points: &mut usize) {
        use TriangleContents::*;

        match self {
            None => *self = Self::one(points),
            One(_) => self.three(points),
            Three { .. } => self.six(points),
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
                self.subdivide(points);
            }
            More {
                sides,
                contents,
                my_side_length,
                ..
            } => {
                *points += 3;
                let len = *points as u32;
                sides.extend_from_slice(&[len - 3, len - 2, len - 1]);
                *my_side_length += 1;

                contents.subdivide(points);
            }
        }
    }

    pub fn calculate(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut [Vec3A],
        shape: &impl BaseShape,
    ) {
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert!(ab.len() >= 2);

        use TriangleContents::*;

        match self {
            None => panic!(),
            One(_) => self.calculate_one(ab, bc, points, shape),
            Three { .. } => self.calculate_three(ab, bc, ca, points, shape),
            Six { .. } => self.calculate_six(ab, bc, ca, points, shape),
            &mut More {
                a: a_idx,
                b: b_idx,
                c: c_idx,
                ref mut sides,
                ref mut contents,
                ref mut my_side_length,
            } => {
                let side_length = *my_side_length as usize;

                let outer_len = ab.len();

                let aba = points[ab[1] as usize];
                let abb = points[ab[outer_len - 2] as usize];
                let bcb = points[bc[1] as usize];
                let bcc = points[bc[outer_len - 2] as usize];
                let cac = points[ca[1] as usize];
                let caa = points[ca[outer_len - 2] as usize];

                points[a_idx as usize] = shape.interpolate_half(aba, caa);
                points[b_idx as usize] = shape.interpolate_half(abb, bcb);
                points[c_idx as usize] = shape.interpolate_half(bcc, cac);

                let ab = &sides[0..side_length];
                let bc = &sides[side_length..side_length * 2];
                let ca = &sides[side_length * 2..];

                shape.interpolate_multiple(
                    points[a_idx as usize],
                    points[b_idx as usize],
                    ab,
                    points,
                );
                shape.interpolate_multiple(
                    points[b_idx as usize],
                    points[c_idx as usize],
                    bc,
                    points,
                );
                shape.interpolate_multiple(
                    points[c_idx as usize],
                    points[a_idx as usize],
                    ca,
                    points,
                );

                contents.calculate(Forward(ab), Forward(bc), Forward(ca), points, shape);
            }
        }
    }

    ///
    /// Indexes the AB edge.
    ///
    /// This is inclusive of A and B.
    ///
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

    ///
    /// Indexes the BC edge.
    ///
    /// This is inclusive of B and C.
    ///
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

    ///
    /// Indexes the CA edge.
    ///
    /// This is inclusive of C and A.
    ///
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

    ///
    /// Adds the indices in this portion of the triangle
    /// to the specified buffer in order.
    ///
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

    ///
    /// Adds the indices for a wireframe of the triangles
    /// in this portion of the triangle to the specified
    /// buffer in order.
    ///
    pub fn add_line_indices(
        &self,
        buffer: &mut Vec<u32>,
        delta: u32,
        mut breaks: impl FnMut(&mut Vec<u32>),
    ) {
        use TriangleContents::*;
        match self {
            None | One(_) | Three { .. } => {}
            &Six { ab, bc, ca, .. } => {
                buffer.extend_from_slice(&[ab + delta, bc + delta, ca + delta]);
                breaks(buffer);
            }
            &More {
                ref sides,
                my_side_length,
                ref contents,
                ..
            } => {
                let my_side_length = my_side_length as usize;
                let ab = &sides[0..my_side_length];
                let bc = &sides[my_side_length..my_side_length * 2];
                let ca = &sides[my_side_length * 2..];

                // Contents are always stored forward.
                add_line_indices_triangular(
                    delta,
                    Forward(ab),
                    Forward(bc),
                    Forward(ca),
                    &**contents,
                    buffer,
                );
                breaks(buffer);
                contents.add_line_indices(buffer, delta, breaks);
            }
        }
    }
}

#[derive(Clone, Debug)]
///
/// A main triangle on the base shape of a subdivided shape.
///
/// The specification of the library expects `a`, `b`, and `c`
/// to be in a counter-clockwise winding.
///
pub struct Triangle {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub ab_edge: usize,
    pub bc_edge: usize,
    pub ca_edge: usize,
    pub(crate) ab_forward: bool,
    pub(crate) bc_forward: bool,
    pub(crate) ca_forward: bool,

    pub(crate) contents: TriangleContents,
}

impl Default for Triangle {
    fn default() -> Self {
        Self {
            a: 0,
            b: 0,
            c: 0,
            ab_edge: 0,
            bc_edge: 0,
            ca_edge: 0,
            ab_forward: false,
            bc_forward: false,
            ca_forward: false,
            contents: TriangleContents::None,
        }
    }
}

impl Triangle {
    ///
    /// Creates a new `Triangle` given the data. This is done
    /// to avoid boilerplate.
    ///
    pub const fn new(
        a: u32,
        b: u32,
        c: u32,
        ab_edge: usize,
        bc_edge: usize,
        ca_edge: usize,
    ) -> Self {
        Self {
            a,
            b,
            c,
            ab_edge,
            bc_edge,
            ca_edge,

            ab_forward: false,
            bc_forward: false,
            ca_forward: false,

            contents: TriangleContents::None,
        }
    }

    fn calculate_edges(
        &mut self,
        edges: &mut [Edge],
        points: &mut [Vec3A],
        shape: &impl BaseShape,
    ) -> usize {
        let mut divide = |p1: u32, p2: u32, edge_idx: usize, forward: &mut bool| {
            if !edges[edge_idx].done {
                shape.interpolate_multiple(
                    points[p1 as usize],
                    points[p2 as usize],
                    &edges[edge_idx].points,
                    points,
                );

                edges[edge_idx].done = true;
                *forward = true;
            } else {
                *forward = false;
            }
        };

        divide(self.a, self.b, self.ab_edge, &mut self.ab_forward);
        divide(self.b, self.c, self.bc_edge, &mut self.bc_forward);
        divide(self.c, self.a, self.ca_edge, &mut self.ca_forward);

        edges[self.ab_edge].points.len()
    }

    ///
    /// Subdivides the edges and contents of this triangle.
    ///
    /// If `calculate` is set to `false`, then the points are
    /// simply added to the buffer and the indices recorded,
    /// but no calculations are performed.
    ///
    fn subdivide(&mut self, points: &mut usize, subdivision_level: usize) {
        if subdivision_level >= 1 {
            self.contents.subdivide(points);
        }
    }

    fn calculate(&mut self, edges: &mut [Edge], points: &mut [Vec3A], shape: &impl BaseShape) {
        let side_length = self.calculate_edges(edges, points, shape) + 1;

        if side_length > 2 {
            let ab = if self.ab_forward {
                Forward(&edges[self.ab_edge].points)
            } else {
                Backward(&edges[self.ab_edge].points)
            };
            let bc = if self.bc_forward {
                Forward(&edges[self.bc_edge].points)
            } else {
                Backward(&edges[self.bc_edge].points)
            };
            let ca = if self.ca_forward {
                Forward(&edges[self.ca_edge].points)
            } else {
                Backward(&edges[self.ca_edge].points)
            };
            self.contents.calculate(ab, bc, ca, points, shape);
        }
    }

    ///
    /// Appends the indices of all the subtriangles onto the
    /// specified buffer.
    ///
    fn add_indices(&self, buffer: &mut Vec<u32>, edges: &[Edge]) {
        let ab = if self.ab_forward {
            Forward(&edges[self.ab_edge].points)
        } else {
            Backward(&edges[self.ab_edge].points)
        };
        let bc = if self.bc_forward {
            Forward(&edges[self.bc_edge].points)
        } else {
            Backward(&edges[self.bc_edge].points)
        };
        let ca = if self.ca_forward {
            Forward(&edges[self.ca_edge].points)
        } else {
            Backward(&edges[self.ca_edge].points)
        };

        add_indices_triangular(self.a, self.b, self.c, ab, bc, ca, &self.contents, buffer);

        self.contents.add_indices(buffer);
    }

    ///
    /// Appends the indices of all the subtriangles' wireframes
    /// onto the specified buffer.
    ///
    fn add_line_indices(
        &self,
        buffer: &mut Vec<u32>,
        edges: &[Edge],
        delta: u32,
        mut breaks: impl FnMut(&mut Vec<u32>),
    ) {
        let ab = if self.ab_forward {
            Forward(&edges[self.ab_edge].points)
        } else {
            Backward(&edges[self.ab_edge].points)
        };
        let bc = if self.bc_forward {
            Forward(&edges[self.bc_edge].points)
        } else {
            Backward(&edges[self.bc_edge].points)
        };
        let ca = if self.ca_forward {
            Forward(&edges[self.ca_edge].points)
        } else {
            Backward(&edges[self.ca_edge].points)
        };

        add_line_indices_triangular(delta, ab, bc, ca, &self.contents, buffer);

        breaks(buffer);

        self.contents.add_line_indices(buffer, delta, breaks);
    }
}

///
/// A progressively subdivided shape which can record
/// the indices of the points and list out the individual
/// triangles of the resulting shape.
///
/// All base triangles specified by `S` in [`BaseShape`]
/// are expected to be in counter clockwise winding.
///
/// Points are preferably stored with coordinates less
/// than or equal to `1.0`. This is why all default shapes
/// lie on the unit sphere.
///
pub struct Subdivided<T, S: BaseShape> {
    points: Vec<Vec3A>,
    data: Vec<T>,
    triangles: Box<[Triangle]>,
    shared_edges: Box<[Edge]>,
    subdivisions: usize,
    shape: S,
}

impl<T, S: BaseShape + Default> Subdivided<T, S> {
    pub fn new(subdivisions: usize, generator: impl FnMut(Vec3A) -> T) -> Self {
        Self::new_custom_shape(subdivisions, generator, Default::default())
    }
}

impl<T, S: BaseShape> Subdivided<T, S> {
    ///
    /// Creates the base shape from `S` and subdivides it.
    ///
    /// - `subdivisions` specifies the number of times a subdivision
    /// will be created. In other terms, this is the number of auxiliary
    /// points between the vertices on the original shape.
    ///
    /// - `generator` is a function run once all the subdivisions are
    /// applied and its values are stored in an internal `Vec`.
    ///
    pub fn new_custom_shape(
        subdivisions: usize,
        generator: impl FnMut(Vec3A) -> T,
        shape: S,
    ) -> Self {
        let mut this = Self {
            points: shape.initial_points(),
            shared_edges: {
                let mut edges = Vec::new();
                edges.resize_with(S::EDGES, Edge::default);
                edges.into_boxed_slice()
            },
            triangles: shape.triangles(),
            subdivisions: 1,
            data: vec![],
            shape,
        };

        let mut new_points = this.points.len();

        for edge in &mut *this.shared_edges {
            edge.subdivide_n_times(subdivisions, &mut new_points);
            edge.done = false;
        }

        for triangle in &mut *this.triangles {
            for i in 0..subdivisions {
                triangle.subdivide(&mut new_points, i);
            }
        }

        let diff = new_points - this.points.len();
        this.points
            .extend(std::iter::repeat(Vec3A::ZERO).take(diff));

        for triangle in &mut *this.triangles {
            triangle.calculate(&mut *this.shared_edges, &mut this.points, &this.shape);
        }

        this.data = this.points.iter().copied().map(generator).collect();

        this
    }

    ///
    /// Subdivides all triangles. `calculate` signals whether or not
    /// to recalculate vertices (To not calculate vertices between many
    /// subdivisions).
    ///
    pub fn subdivide(&mut self) {
        for Edge { done, .. } in &mut *self.shared_edges {
            *done = false;
        }

        let mut new_points = self.points.len();

        let subdivision_level = self.shared_edges[0].points.len();

        for edge in &mut *self.shared_edges {
            edge.subdivide_n_times(1, &mut new_points);
            edge.done = false;
        }

        for triangle in &mut *self.triangles {
            triangle.subdivide(&mut new_points, subdivision_level);
        }

        let diff = new_points - self.points.len();
        self.points
            .extend(std::iter::repeat(Vec3A::ZERO).take(diff));

        for triangle in &mut *self.triangles {
            triangle.calculate(&mut *self.shared_edges, &mut self.points, &self.shape);
        }
    }

    ///
    /// The raw points created by the subdivision process.
    ///
    pub fn raw_points(&self) -> &[Vec3A] {
        &self.points
    }

    ///
    /// Appends the indices for the triangle into `buffer`.
    ///
    /// The specified triangle is a main triangle on the base
    /// shape. The range of this should be limited to the number
    /// of triangles in the base shape.
    ///
    /// Alternatively, use [`get_all_indices`] to get all the
    /// indices.
    ///
    /// [`get_all_indices`]: #method.get_all_indices
    ///
    pub fn get_indices(&self, triangle: usize, buffer: &mut Vec<u32>) {
        self.triangles[triangle].add_indices(buffer, &self.shared_edges);
    }

    ///
    /// Gets the indices for all main triangles in the shape.
    ///
    pub fn get_all_indices(&self) -> Vec<u32> {
        let mut buffer = Vec::new();

        for i in 0..self.triangles.len() {
            self.get_indices(i, &mut buffer);
        }

        buffer
    }

    ///
    /// Gets the wireframe indices for the contents of a specified triangle.
    ///
    /// `delta` is added to all of the indices pushed into the buffer, and
    /// is generally intended to be used to have a NaN vertex at zero. Set
    /// to zero to produce the indices as if there was no NaN vertex.
    ///
    /// `breaks` is run every time there is a necessary break in the line
    /// strip. Use this to, for example, swap out the buffer using
    /// [`std::mem::swap`], or push a NaN index into the buffer.
    ///
    pub fn get_line_indices(
        &self,
        buffer: &mut Vec<u32>,
        triangle: usize,
        delta: usize,
        breaks: impl FnMut(&mut Vec<u32>),
    ) {
        self.triangles[triangle].add_line_indices(buffer, &self.shared_edges, delta as u32, breaks);
    }

    ///
    /// Gets the wireframe indices for the specified edge of the base shape.
    ///
    /// See [`Self::get_line_indices`] for more on `delta`.
    ///
    pub fn get_major_edge_line_indices(&self, edge: usize, buffer: &mut Vec<u32>, delta: usize) {
        buffer.extend(
            self.shared_edges[edge]
                .points
                .iter()
                .map(|x| x + delta as u32),
        );
    }

    ///
    /// Gets the wireframe indices for all main triangles in the shape,
    /// as well as all edges.
    ///
    /// See [`Self::get_line_indices`] for more on `delta`, and `breaks`.
    ///
    pub fn get_all_line_indices(
        &self,
        delta: usize,
        mut breaks: impl FnMut(&mut Vec<u32>),
    ) -> Vec<u32> {
        let mut buffer = Vec::new();

        for i in 0..self.triangles.len() {
            self.get_line_indices(&mut buffer, i, delta, &mut breaks);
        }

        for i in 0..self.shared_edges.len() {
            self.get_major_edge_line_indices(i, &mut buffer, delta);
            breaks(&mut buffer);
        }

        buffer
    }

    ///
    /// Returns the number of subdivisions applied when this shape
    /// was created.
    ///
    pub fn subdivisions(&self) -> usize {
        self.subdivisions
    }

    ///
    /// Returns the custom data created by the generator function.
    ///
    pub fn raw_data(&self) -> &[T] {
        &self.data
    }

    ///
    /// Returns the mutable custom data created by the generator function.
    ///
    pub fn raw_data_mut(&mut self) -> &mut [T] {
        &mut self.data
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
    /// subdivisions * EDGES + INITIAL_POINTS
    /// ```
    ///
    pub fn shared_vertices(&self) -> usize {
        self.subdivisions * S::EDGES + self.shape.initial_points().len()
    }

    ///
    /// Linear distance between two points on this shape.
    ///
    pub fn linear_distance(&self, p1: u32, p2: u32, radius: f32) -> f32 {
        (self.points[p1 as usize] - self.points[p2 as usize]).length() * radius
    }
}

#[cfg(feature = "shape-extras")]
impl<T, S: BaseShape + EquilateralBaseShape> Subdivided<T, S> {
    ///
    /// Closest "main" triangle.
    ///
    /// Undefined results if the point is one of the vertices
    /// on the original base shape.
    ///
    pub fn main_triangle_intersect(point: Vec3A) -> usize {
        let point = point.normalize();
        let mut nearest = 0;
        let mut near_factor = point.dot(S::triangle_normals()[0]);

        if near_factor > S::triangle_min_dot() {
            return 0;
        }

        for (index, normal) in S::triangle_normals().iter().enumerate().skip(1) {
            let factor = normal.dot(point);
            if factor > near_factor {
                if factor > S::triangle_min_dot() {
                    return index;
                }
                nearest = index;
                near_factor = factor;
            }
        }

        nearest
    }

    ///
    /// Distance between two points on this sphere (assuming this
    /// is a sphere).
    ///
    pub fn spherical_distance(&self, p1: u32, p2: u32, radius: f32) -> f32 {
        self.points[p1 as usize]
            .dot(self.points[p2 as usize])
            .acos()
            * radius
    }
}

///
/// Adds the indices of the triangles in this "layer" of the triangle to
/// the buffer.
///
// The logic in this function has been worked out mostly on pen and paper
// and therefore it is difficult to read.
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

///
/// Adds the indices of the triangles in this "layer" of the triangle to
/// the buffer in line strip format.
///
/// This is used to create a wireframe look.
///
// Like the previous function, this logic has been worked out mostly on
// pen and paper and it is therefore difficult to read.
fn add_line_indices_triangular(
    delta: u32,
    ab: Slice<u32>,
    bc: Slice<u32>,
    ca: Slice<u32>,
    contents: &TriangleContents,
    buffer: &mut Vec<u32>,
) {
    if ab.len() == 1 {
        buffer.extend_from_slice(&[ab[0] + delta, bc[0] + delta, ca[0] + delta]);
        return;
    }

    buffer.reserve((ab.len() - 1) * 9);

    if ab.len() != 2 {
        for i in 0..ab.len() - 2 {
            buffer.push(contents.idx_ab(i) + delta);
        }

        for i in 0..bc.len() - 2 {
            buffer.push(contents.idx_bc(i) + delta);
        }

        for i in 0..ca.len() - 2 {
            buffer.push(contents.idx_ca(i) + delta);
        }
    } else {
        buffer.push(contents.idx_ca(0) + delta);
    }

    buffer.push(ab[0] + delta);

    for i in (1..ca.len()).rev() {
        buffer.push(ca[i] + delta);
        buffer.push(contents.idx_ca(i - 1) + delta);
    }

    buffer.push(ca[0] + delta);

    for i in (1..bc.len()).rev() {
        buffer.push(bc[i] + delta);
        buffer.push(contents.idx_bc(i - 1) + delta);
    }

    buffer.push(bc[0] + delta);

    for i in (1..ab.len()).rev() {
        buffer.push(ab[i] + delta);
        buffer.push(contents.idx_ab(i - 1) + delta);
    }
}

#[cfg(feature = "adjacency")]
///
/// Implements neighbour tracking.
///
mod adjacency {
    use tinyvec::ArrayVec;

    #[derive(Copy, Clone, Eq, PartialEq, Debug)]
    pub(crate) enum RehexState {
        Empty,
        Clear,
        TwoTwo,
        ThreeTwo,
        TwoTwoTwo,
        Complete,
    }

    /// Tracks the neighbours adjacent to each vertex by only using index data.
    ///
    /// The result preserves winding: the resulting array is wound around the
    /// center vertex in the same way that the source triangles were wound.
    pub struct AdjacencyBuilder {
        pub(crate) state: Vec<RehexState>,
        pub result: Vec<ArrayVec<[usize; 6]>>,
    }

    impl AdjacencyBuilder {
        pub fn new(points_len: usize) -> Self {
            let state = std::iter::repeat(RehexState::Empty)
                .take(points_len)
                .collect::<Vec<_>>();
            let result = std::iter::repeat(ArrayVec::new())
                .take(points_len)
                .collect::<Vec<_>>();
            Self { state, result }
        }

        pub fn add_indices(&mut self, indices: &[u32]) {
            for chunk in indices.chunks_exact(3) {
                let &[a, b, c] = chunk else { unreachable!() };

                self.add_triangle(a, b, c);
                self.add_triangle(c, a, b);
                self.add_triangle(b, c, a);
            }
        }

        pub fn finish(self) -> Vec<ArrayVec<[usize; 6]>> {
            self.result
        }

        fn add_triangle(&mut self, a: u32, b: u32, c: u32) {
            let (a, b, c) = (a as usize, b as usize, c as usize);
            let state = &mut self.state[a];
            if let RehexState::Complete = state {
                return;
            }

            let result = &mut self.result[a];

            match state {
                RehexState::Empty => {
                    result.extend([b, c]);
                    *state = RehexState::Clear;
                }
                RehexState::Clear => {
                    if result[result.len() - 1] == b {
                        if result[0] == c {
                            *state = RehexState::Complete;
                        } else {
                            result.push(c);
                            if result.len() == 6 {
                                *state = RehexState::Complete;
                            }
                        }
                    } else if result[0] == c {
                        result.insert(0, b);
                        if result.len() == 6 {
                            *state = RehexState::Complete;
                        }
                    } else {
                        *state = match result.len() {
                            2 => RehexState::TwoTwo,
                            3 => RehexState::ThreeTwo,
                            4 => RehexState::Complete,
                            _ => unreachable!(),
                        };
                        result.extend([b, c]);
                    }
                }
                RehexState::TwoTwo => {
                    if result[1] == b {
                        if result[2] == c {
                            *state = RehexState::Clear;
                        } else {
                            result.insert(2, c);
                            *state = RehexState::ThreeTwo;
                        }
                    } else if result[0] == c {
                        if result[3] == b {
                            let temp = result[2];
                            result.pop();
                            result.pop();
                            result.insert(0, temp);
                            result.insert(1, b);
                            *state = RehexState::Clear;
                        } else {
                            result.insert(0, b);
                            *state = RehexState::ThreeTwo;
                        }
                    } else if result[2] == c {
                        result.insert(0, b);
                        let t2 = result.swap_remove(2);
                        let t1 = result.swap_remove(1);
                        result.push(t1);
                        result.push(t2);
                        *state = RehexState::ThreeTwo;
                    } else {
                        result.extend([b, c]);
                        *state = RehexState::TwoTwoTwo;
                    }
                }
                RehexState::ThreeTwo => {
                    if result[2] == b {
                        if result[3] == c {
                            *state = RehexState::Clear;
                        } else {
                            result.insert(3, c);
                            *state = RehexState::Complete;
                        }
                    } else {
                        if result[4] == b {
                            result.pop();
                            let temp = result.pop().unwrap();
                            result.insert(0, b);
                            result.insert(0, temp);
                            *state = RehexState::Clear;
                        } else {
                            result.insert(0, b);
                            *state = RehexState::Complete;
                        }
                    }
                }
                RehexState::TwoTwoTwo => {
                    if (result[1] != b || result[2] != c)
                        && (result[3] != b || result[4] != c)
                        && (result[5] != b || result[0] != c)
                    {
                        let t2 = result.swap_remove(3);
                        let t1 = result.swap_remove(2);
                        result.extend([t1, t2]);
                    }
                    *state = RehexState::Complete;
                }
                RehexState::Complete => unreachable!(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::shapes::IcoSphere;
    use crate::Slice::Forward;
    use glam::Vec3A;

    // Starting points aren't _quite_ precise enough to use `f32::EPSILON`.
    const EPSILON: f32 = 0.0000002;

    #[test]
    fn slerp_one() {
        use crate::interpolation::geometric_slerp_half;
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
        use crate::interpolation::geometric_slerp_multiple;

        let p1 = Vec3A::new(0.0, -0.885330189449, 0.464962854054);
        let p2 = Vec3A::new(0.0, 0.946042343528, 0.324043028395);

        let expected = Vec3A::new(0.0, 0.0767208624118, 0.997052611085);

        let mut result = Vec3A::ZERO;
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

        let mut result = [Vec3A::ZERO, Vec3A::ZERO, Vec3A::ZERO, Vec3A::ZERO];

        geometric_slerp_multiple(p1, p2, &[0, 1, 2, 3], &mut result);

        for (&expected, &result) in expected.iter().zip(result.iter()) {
            assert!((expected - result).length() <= EPSILON);
        }
    }

    #[test]
    fn new() {
        let x = IcoSphere::new(0, |_| ());
        x.get_indices(0, &mut Vec::new());
    }

    #[test]
    fn one() {
        let x = IcoSphere::new(1, |_| ());
        x.get_indices(0, &mut Vec::new());
    }

    #[test]
    fn second_layer_inner() {
        let x = IcoSphere::new(2, |_| ());
        x.get_indices(0, &mut Vec::new());
        let x = IcoSphere::new(3, |_| ());
        x.get_indices(0, &mut Vec::new());
        let x = IcoSphere::new(5, |_| ());
        x.get_indices(0, &mut Vec::new());
        let x = IcoSphere::new(6, |_| ());
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
        let sphere = IcoSphere::new(10, |_| ());

        for i in sphere.raw_points() {
            assert!(i.length() - 1.0 <= EPSILON);
        }
    }

    #[test]
    fn line_one() {
        use super::add_line_indices_triangular;
        use super::TriangleContents;

        let mut buffer = Vec::new();

        add_line_indices_triangular(
            0,
            Forward(&[0]),
            Forward(&[1]),
            Forward(&[2]),
            &TriangleContents::none(),
            &mut buffer,
        );

        assert_eq!(buffer, &[0, 1, 2]);
    }

    #[test]
    fn line_two() {
        use super::add_line_indices_triangular;
        use super::TriangleContents;

        let mut buffer = Vec::new();

        add_line_indices_triangular(
            0,
            Forward(&[0, 1]),
            Forward(&[2, 3]),
            Forward(&[4, 5]),
            &TriangleContents::One(6),
            &mut buffer,
        );

        assert_eq!(buffer, &[6, 0, 5, 6, 4, 3, 6, 2, 1, 6]);
    }

    #[test]
    fn line_three() {
        use super::add_line_indices_triangular;
        use super::TriangleContents;

        let mut buffer = Vec::new();

        add_line_indices_triangular(
            0,
            Forward(&[0, 1, 2]),
            Forward(&[3, 4, 5]),
            Forward(&[6, 7, 8]),
            &TriangleContents::Three { a: 9, b: 10, c: 11 },
            &mut buffer,
        );

        assert_eq!(
            buffer,
            &[9, 10, 11, 0, 8, 9, 7, 11, 6, 5, 11, 4, 10, 3, 2, 10, 1, 9]
        );
    }

    #[cfg(feature = "adjacency")]
    mod adjacency {
        use crate::adjacency::RehexState;
        use crate::{adjacency::AdjacencyBuilder, shapes::IcoSphere};

        #[test]
        fn creation() {
            let sphere = IcoSphere::new(5, |_| ());

            let mut indices = Vec::new();

            for i in 0..20 {
                sphere.get_indices(i, &mut indices);
            }

            let mut builder = AdjacencyBuilder::new(sphere.raw_points().len());
            builder.add_indices(&indices);
            builder
                .state
                .iter()
                .for_each(|&state| assert_eq!(state, RehexState::Complete));
        }
    }
}
