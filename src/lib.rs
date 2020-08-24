//!
//! Library for hexagonally tiling a sphere.
//!
//! Note this is *not* an explicit rewrite of `hexasphere.js`.
//!

use glam::Vec3A;
use core::ops::Index;
use core::marker::PhantomData;

use Slice::*;

#[cfg(feature = "adjacency")]
pub use adjacency::*;

pub trait BaseShape {
    fn initial_points() -> &'static [Vec3A];
    fn triangles() -> &'static [Triangle];
    const EDGES: usize;

    fn interpolate(a: Vec3A, b: Vec3A, p: f32) -> Vec3A;

    fn interpolate_half(a: Vec3A, b: Vec3A) -> Vec3A {
        Self::interpolate(a, b, 0.5)
    }

    fn interpolate_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        for (percent, index) in indices.iter().enumerate() {
            let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

            points[*index as usize] = Self::interpolate(a, b, percent);
        }
    }
}

pub trait EquilateralBaseShape: BaseShape {
    fn triangle_normals() -> &'static [Vec3A];
    fn triangle_min_dot() -> f32;
}

pub struct IcoSphereBase;

impl BaseShape for IcoSphereBase {
    #[inline]
    fn initial_points() -> &'static [Vec3A] {
        &*consts::icosphere::INITIAL_POINTS
    }

    #[inline]
    fn triangles() -> &'static [Triangle] {
        &consts::icosphere::TRIANGLES
    }
    const EDGES: usize = consts::icosphere::EDGES;

    #[inline]
    fn interpolate(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        geometric_slerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(a: Vec3A, b: Vec3A) -> Vec3A {
        geometric_slerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A])  {
        geometric_slerp_multiple(a, b, indices, points);
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

pub type Hexasphere<T> = Subdivided<T, IcoSphereBase>;

pub struct TetraSphereBase;

impl BaseShape for TetraSphereBase {
    #[inline]
    fn initial_points() -> &'static [Vec3A] {
        &*consts::tetrasphere::INITIAL_POINTS
    }

    #[inline]
    fn triangles() -> &'static [Triangle] {
        &consts::tetrasphere::TRIANGLES
    }
    const EDGES: usize = consts::tetrasphere::EDGES;

    #[inline]
    fn interpolate(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        geometric_slerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(a: Vec3A, b: Vec3A) -> Vec3A {
        geometric_slerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        geometric_slerp_multiple(a, b, indices, points);
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

pub type TetraSphere<T> = Subdivided<T, TetraSphereBase>;

pub struct TriangleBase;

impl BaseShape for TriangleBase {
    #[inline]
    fn initial_points() -> &'static [Vec3A] {
        &*consts::triangle::INITIAL_POINTS
    }

    #[inline]
    fn triangles() -> &'static [Triangle] {
        core::slice::from_ref(&consts::triangle::TRIANGLE)
    }
    const EDGES: usize = consts::triangle::EDGES;

    #[inline]
    fn interpolate(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        lerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(a: Vec3A, b: Vec3A) -> Vec3A {
        lerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        lerp_multiple(a, b, indices, points);
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

pub type TrianglePlane<T> = Subdivided<T, TriangleBase>;

pub struct SquareBase;

impl BaseShape for SquareBase {
    #[inline]
    fn initial_points() -> &'static [Vec3A] {
        &*consts::square::INITIAL_POINTS
    }

    #[inline]
    fn triangles() -> &'static [Triangle] {
        &consts::square::TRIANGLES
    }
    const EDGES: usize = consts::square::EDGES;

    #[inline]
    fn interpolate(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
        lerp(a, b, p)
    }

    #[inline]
    fn interpolate_half(a: Vec3A, b: Vec3A) -> Vec3A {
        lerp_half(a, b)
    }

    #[inline]
    fn interpolate_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
        lerp_multiple(a, b, indices, points);
    }
}

pub type SquarePlane<T> = Subdivided<T, SquareBase>;

struct Edge {
    points: Vec<u32>,
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
    fn one<S: BaseShape>(ab: Slice<u32>, bc: Slice<u32>, points: &mut Vec<Vec3A>, calculate: bool) -> Self {
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), 2);
        let p1 = points[ab[0] as usize];
        let p2 = points[bc[1] as usize];
        let index = points.len() as u32;
        if calculate {
            points.push(S::interpolate_half(p1, p2));
        } else {
            points.push(Vec3A::zero());
        }
        TriangleContents::One(index)
    }

    ///
    /// Creates a `Three` variant from a `One` variant.
    ///
    fn three<S: BaseShape>(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut Vec<Vec3A>,
        calculate: bool,
    ) {
        use TriangleContents::*;

        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert_eq!(ab.len(), 3);

        match self {
            &mut One(x) => {
                let ab = points[ab[1] as usize];
                let bc = points[bc[1] as usize];
                let ca = points[ca[1] as usize];

                if calculate {
                    let a = S::interpolate_half(ab, ca);
                    let b = S::interpolate_half(bc, ab);
                    let c = S::interpolate_half(ca, bc);

                    points.extend_from_slice(&[b, c]);
                    points[x as usize] = a;
                } else {
                    points.extend_from_slice(&[Vec3A::zero(), Vec3A::zero()])
                }

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
    fn six<S: BaseShape>(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut Vec<Vec3A>,
        calculate: bool,
    ) {
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

                if calculate {
                    let a = S::interpolate_half(aba, caa);
                    let b = S::interpolate_half(abb, bcb);
                    let c = S::interpolate_half(bcc, cac);

                    let ab = S::interpolate_half(a, b);
                    let bc = S::interpolate_half(b, c);
                    let ca = S::interpolate_half(c, a);

                    points[a_index as usize] = a;
                    points[b_index as usize] = b;
                    points[c_index as usize] = c;
                    points.extend_from_slice(&[ab, bc, ca]);
                } else {
                    points.extend_from_slice(&[Vec3A::zero(), Vec3A::zero(), Vec3A::zero()])
                }

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
    pub fn subdivide<S: BaseShape>(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut Vec<Vec3A>,
        calculate: bool,
    ) {
        use TriangleContents::*;
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert!(ab.len() >= 2);
        match self {
            None => *self = Self::one::<S>(ab, bc, points, calculate),
            One(_) => self.three::<S>(ab, bc, ca, points, calculate),
            Three { .. } => self.six::<S>(ab, bc, ca, points, calculate),
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
                self.subdivide::<S>(ab, bc, ca, points, calculate);
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

                if calculate {
                    points[a_idx as usize] = S::interpolate_half(aba, caa);
                    points[b_idx as usize] = S::interpolate_half(abb, bcb);
                    points[c_idx as usize] = S::interpolate_half(bcc, cac);
                }

                let ab = &sides[0..side_length];
                let bc = &sides[side_length..side_length * 2];
                let ca = &sides[side_length * 2..];

                if calculate {
                    S::interpolate_multiple(
                        points[a_idx as usize],
                        points[b_idx as usize],
                        ab,
                        points,
                    );
                    S::interpolate_multiple(
                        points[b_idx as usize],
                        points[c_idx as usize],
                        bc,
                        points,
                    );
                    S::interpolate_multiple(
                        points[c_idx as usize],
                        points[a_idx as usize],
                        ca,
                        points,
                    );
                }

                contents.subdivide::<S>(Forward(ab), Forward(bc), Forward(ca), points, calculate);
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

#[derive(Clone, Debug)]
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
    pub fn new(a: u32, b: u32, c: u32, ab_edge: usize, bc_edge: usize, ca_edge: usize) -> Self {
        Self {
            a,
            b,
            c,
            ab_edge,
            bc_edge,
            ca_edge,
            ..Default::default()
        }
    }

    fn subdivide_edges<'a, S: BaseShape>(
        &'a mut self,
        edges: &mut [Edge],
        points: &mut Vec<Vec3A>,
        calculate: bool,
    ) -> usize {
        let mut divide = |p1: u32, p2: u32, edge_idx: usize, forward: &mut bool| {
            if !edges[edge_idx].done {
                edges[edge_idx].points.push(points.len() as u32);
                points.push(Vec3A::zero());

                if calculate {
                    S::interpolate_multiple(
                        points[p1 as usize],
                        points[p2 as usize],
                        &edges[edge_idx].points,
                        points,
                    );
                }

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

    fn subdivide<S: BaseShape>(
        &mut self,
        edges: &mut [Edge],
        points: &mut Vec<Vec3A>,
        calculate: bool,
    ) {
        let side_length = self.subdivide_edges::<S>(edges, points, calculate) + 1;

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
            self.contents.subdivide::<S>(ab, bc, ca, points, calculate);
        }
    }

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
pub struct Subdivided<T, S: BaseShape> {
    points: Vec<Vec3A>,
    data: Vec<T>,
    triangles: Box<[Triangle]>,
    shared_edges: Box<[Edge]>,
    subdivisions: usize,
    _phantom: PhantomData<S>,
}

impl<T, S: BaseShape> Subdivided<T, S> {
    ///
    /// Creates a hexagonally tiled sphere.
    ///
    /// The number of hexagons between pentagons is `subdivisions`, and
    /// `generator` permits you to create a unique value for every point
    /// on the sphere.
    ///
    pub fn new(subdivisions: usize, generator: impl FnMut(Vec3A) -> T) -> Self {
        let mut this = Self {
            points: S::initial_points().into(),
            shared_edges: {
                let mut edges = Vec::new();
                edges.resize_with(S::EDGES, Edge::default);
                edges.into_boxed_slice()
            },
            triangles: S::triangles().to_vec().into_boxed_slice(),
            subdivisions: 1,
            data: vec![],
            _phantom: PhantomData
        };

        match subdivisions {
            0 => {}
            1 => this.subdivide(true),
            x => {
                for _ in 0..x - 1 {
                    this.subdivide(false);
                }

                this.subdivide(true);
            }
        }

        this.data = this.points.iter().copied().map(generator).collect();

        this
    }

    ///
    /// Subdivides all triangles. `calculate` signals whether or not
    /// to recalculate vertices (To not calculate vertices between many
    /// subdivisions)
    ///
    fn subdivide(&mut self, calculate: bool) {
        for Edge { done, .. } in &mut *self.shared_edges {
            *done = false;
        }

        for triangle in &mut *self.triangles {
            triangle.subdivide::<S>(&mut *self.shared_edges, &mut self.points, calculate);
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
        self.subdivisions * S::EDGES + S::initial_points().len()
    }

    ///
    /// Linear distance between two points on this sphere.
    ///
    pub fn linear_distance(&self, p1: u32, p2: u32, radius: f32) -> f32 {
        (self.points[p1 as usize] - self.points[p2 as usize]).length() * radius
    }
}

impl<T, S: BaseShape + EquilateralBaseShape> Subdivided<T, S> {
    ///
    /// Closest "main" triangle.
    ///
    /// Undefined results if the point is one of the vertices
    /// on the original icosahedron.
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
    /// Distance between two points on this sphere.
    ///
    pub fn spherical_distance(&self, p1: u32, p2: u32, radius: f32) -> f32 {
        self.points[p1 as usize]
            .dot(self.points[p2 as usize])
            .acos()
            * radius
    }
}

impl<T> Hexasphere<T> {
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

/// Note: `a` and `b` should both be normalized for normalized results.
pub fn geometric_slerp(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
    let angle = a.dot(b).acos();

    let sin = angle.sin().recip();
    a * (((1.0 - p) * angle).sin() * sin) + b * ((p * angle).sin() * sin)
}

/// Note: `a` and `b` should both be normalized for normalized results.
pub fn geometric_slerp_half(a: Vec3A, b: Vec3A) -> Vec3A {
    let angle = a.dot(b).acos();

    let sin_denom = angle.sin().recip();
    let sin_numer = (angle * 0.5).sin();

    (a + b) * sin_denom * sin_numer
}

/// Note: `a` and `b` should both be normalized for normalized results.
pub fn geometric_slerp_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
    let angle = a.dot(b).acos();
    let sin = angle.sin().recip();

    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] =
            a * (((1.0 - percent) * angle).sin() * sin) + b * ((percent * angle).sin() * sin);
    }
}

pub fn normalized_lerp(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
    ((1.0 - p) * a + p * b).normalize()
}

pub fn normalized_lerp_half(a: Vec3A, b: Vec3A) -> Vec3A {
    ((a + b) * 0.5).normalize()
}

pub fn normalized_lerp_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] = ((1.0 - percent) * a + percent * b).normalize();
    }
}

pub fn lerp(a: Vec3A, b: Vec3A, p: f32) -> Vec3A {
    (1.0 - p) * a + p * b
}

pub fn lerp_half(a: Vec3A, b: Vec3A) -> Vec3A {
    (a + b) * 0.5
}

pub fn lerp_multiple(a: Vec3A, b: Vec3A, indices: &[u32], points: &mut [Vec3A]) {
    for (percent, index) in indices.iter().enumerate() {
        let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

        points[*index as usize] = (1.0 - percent) * a + percent * b;
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

#[cfg(feature = "adjacency")]
mod adjacency {
    use smallvec::SmallVec;
    use std::collections::HashMap;

    #[derive(Default, Clone, Debug)]
    pub struct AdjacentStore {
        pub(crate) subdivisions: usize,
        pub(crate) map: HashMap<u32, SmallVec<[u32; 6]>>,
    }

    impl AdjacentStore {
        pub fn new() -> Self {
            Self::default()
        }

        pub fn neighbours(&self, id: u32) -> Option<&[u32]> {
            self.map.get(&id).map(|x| &**x)
        }

        pub fn from_indices(indices: &[u32]) -> Self {
            let mut this = Self::new();
            this.add_triangle_indices(indices);
            this
        }

        pub fn add_triangle_indices(&mut self, triangles: &[u32]) {
            assert_eq!(triangles.len() % 3, 0);

            for triangle in triangles.chunks(3) {
                self.add_triangle([triangle[0], triangle[1], triangle[2]]);
            }
        }

        fn add_triangle(&mut self, [a, b, c]: [u32; 3]) {
            let mut add_triangle = |a, b, c| {
                let vec = self.map.entry(a).or_insert_with(SmallVec::new);
                if !vec.contains(&b) {
                    vec.push(b);
                }
                if !vec.contains(&c) {
                    vec.push(c);
                }
            };

            add_triangle(a, b, c);
            add_triangle(b, c, a);
            add_triangle(c, a, b);
        }
    }
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

    #[cfg(feature = "adjacency")]
    mod adjacency {
        use crate::{AdjacentStore, Hexasphere};

        #[test]
        fn creation() {
            let sphere = Hexasphere::new(0, |_| ());

            let mut indices = Vec::new();

            for i in 0..20 {
                sphere.get_indices(i, &mut indices);
            }

            let _ = AdjacentStore::from_indices(&indices);
        }

        #[test]
        fn correct_indices() {
            let sphere = Hexasphere::new(0, |_| ());

            let mut indices = Vec::new();

            for i in 0..20 {
                sphere.get_indices(i, &mut indices);
            }

            let store = AdjacentStore::from_indices(&indices);

            const REFERENCE_DATA: [[u32; 5]; 12] = [
                [1, 2, 3, 4, 5],
                [0, 2, 7, 6, 5],
                [0, 1, 3, 8, 7],
                [0, 4, 2, 9, 8],
                [0, 5, 10, 9, 3],
                [0, 1, 6, 10, 4],
                [5, 1, 7, 11, 10],
                [1, 2, 8, 11, 6],
                [2, 3, 9, 11, 7],
                [3, 4, 10, 11, 8],
                [4, 5, 6, 11, 9],
                [6, 7, 8, 9, 10],
            ];

            for i in 0..12 {
                let expected = REFERENCE_DATA[i as usize];
                let actual = store.neighbours(i).unwrap();
                assert_eq!(actual.len(), 5);
                let mut values = [0; 5];
                for (x, i) in actual.iter().enumerate() {
                    assert!(expected.contains(i));
                    values[x] += 1;
                }
                assert_eq!(values, [1; 5]);
            }
        }
    }
}

mod consts {
    pub mod square {
        use crate::{Triangle, TriangleContents};
        use glam::Vec3A;

        lazy_static::lazy_static! {
            pub(crate) static ref INITIAL_POINTS: [Vec3A; 4] = [
                Vec3A::new(std::f32::consts::FRAC_1_SQRT_2, 0.0, std::f32::consts::FRAC_1_SQRT_2),
                Vec3A::new(std::f32::consts::FRAC_1_SQRT_2, 0.0, -std::f32::consts::FRAC_1_SQRT_2),
                Vec3A::new(-std::f32::consts::FRAC_1_SQRT_2, 0.0, -std::f32::consts::FRAC_1_SQRT_2),
                Vec3A::new(-std::f32::consts::FRAC_1_SQRT_2, 0.0, std::f32::consts::FRAC_1_SQRT_2),
            ];

            pub(crate) static ref TRIANGLE_NORMALS: [Vec3A; 2] = [
                Vec3A::new(0.0, 1.0, 0.0),
                Vec3A::new(0.0, 1.0, 0.0),
            ];
        }

        pub(crate) const TRIANGLES: [Triangle; 2] = [
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

        pub(crate) const EDGES: usize = 5;
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

        pub(crate) const TRIANGLE: Triangle = Triangle {
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

        pub(crate) const EDGES: usize = 3;
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

        pub(crate) const TRIANGLES: [Triangle; 4] = [
            Triangle {
                a: 0,
                b: 1,
                c: 3,

                ab_edge: 0,
                bc_edge: 4,
                ca_edge: 3,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 1,
                b: 2,
                c: 3,

                ab_edge: 1,
                bc_edge: 5,
                ca_edge: 4,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 2,
                b: 0,
                c: 3,

                ab_edge: 2,
                bc_edge: 3,
                ca_edge: 5,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,
            },
            Triangle {
                a: 0,
                b: 2,
                c: 1,

                ab_edge: 2,
                bc_edge: 1,
                ca_edge: 0,
                ab_forward: false,
                bc_forward: false,
                ca_forward: false,
                contents: TriangleContents::None,

            }
        ];
        pub(crate) const EDGES: usize = 6;
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

        pub(crate) const TRIANGLES: [Triangle; 20] = [
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

        pub(crate) const EDGES: usize = 30;
    }
}
