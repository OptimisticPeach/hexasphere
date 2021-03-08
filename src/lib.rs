//!
//! Library for subdividing shapes made of triangles.
//!
//! This is for internal bevy use only, and is not
//! well documented. See [`hexasphere`](crates.io/crates/hexasphere)
//! for more information.
//!

use slice::*;
use slice::Slice::*;

pub mod interpolation;
pub mod shapes;
mod slice;

pub trait Vec3:
    std::ops::Add<Self, Output = Self> +
    std::ops::Mul<f32, Output = Self> +
    Copy
{
    const ZERO: Self;

    fn dot(self, other: Self) -> f32;
    fn normalize(self) -> Self;
    fn from_arr3(data: [f32; 3]) -> Self;
}

pub trait BaseShape<V: Vec3> {
    fn initial_points(&self) -> Vec<V>;

    fn triangles(&self) -> Box<[Triangle]>;

    const EDGES: usize;

    fn interpolate(&self, a: V, b: V, p: f32) -> V;


    fn interpolate_half(&self, a: V, b: V) -> V {
        self.interpolate(a, b, 0.5)
    }

    fn interpolate_multiple(&self, a: V, b: V, indices: &[u32], points: &mut [V]) {
        for (percent, index) in indices.iter().enumerate() {
            let percent = (percent + 1) as f32 / (indices.len() + 1) as f32;

            points[*index as usize] = self.interpolate(a, b, percent);
        }
    }
}

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

#[derive(Clone, Debug)]
enum TriangleContents {
    None,
    One(u32),
    Three { a: u32, b: u32, c: u32 },
    Six {
        a: u32,
        b: u32,
        c: u32,
        ab: u32,
        bc: u32,
        ca: u32,
    },
    More {
        a: u32,
        b: u32,
        c: u32,
        // Separated into three `my_side_length` segments
        // to save on extra allocations.
        sides: Vec<u32>,
        my_side_length: u32,
        // Implementing this as a `Vec<TriangleContents>` would
        // probably be a perf. improvement someday, however not
        // something worth implementing right now.
        contents: Box<TriangleContents>,
    },
}

impl TriangleContents {
    pub fn none() -> Self {
        Self::None
    }

    fn one<V: Vec3>(ab: Slice<u32>, bc: Slice<u32>, points: &mut Vec<V>, calculate: bool, shape: &impl BaseShape<V>) -> Self {
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), 2);
        let p1 = points[ab[0] as usize];
        let p2 = points[bc[1] as usize];
        let index = points.len() as u32;
        if calculate {
            points.push(shape.interpolate_half(p1, p2));
        } else {
            points.push(V::ZERO);
        }
        TriangleContents::One(index)
    }

    fn three<V: Vec3>(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut Vec<V>,
        calculate: bool,
        shape: &impl BaseShape<V>,
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
                    let a = shape.interpolate_half(ab, ca);
                    let b = shape.interpolate_half(bc, ab);
                    let c = shape.interpolate_half(ca, bc);

                    points.extend_from_slice(&[b, c]);
                    points[x as usize] = a;
                } else {
                    points.extend_from_slice(&[V::ZERO, V::ZERO])
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

    fn six<V: Vec3>(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut Vec<V>,
        calculate: bool,
        shape: &impl BaseShape<V>,
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
                    let a = shape.interpolate_half(aba, caa);
                    let b = shape.interpolate_half(abb, bcb);
                    let c = shape.interpolate_half(bcc, cac);

                    let ab = shape.interpolate_half(a, b);
                    let bc = shape.interpolate_half(b, c);
                    let ca = shape.interpolate_half(c, a);

                    points[a_index as usize] = a;
                    points[b_index as usize] = b;
                    points[c_index as usize] = c;
                    points.extend_from_slice(&[ab, bc, ca]);
                } else {
                    points.extend_from_slice(&[V::ZERO, V::ZERO, V::ZERO])
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

    pub fn subdivide<V: Vec3>(
        &mut self,
        ab: Slice<u32>,
        bc: Slice<u32>,
        ca: Slice<u32>,
        points: &mut Vec<V>,
        calculate: bool,
        shape: &impl BaseShape<V>,
    ) {
        use TriangleContents::*;
        assert_eq!(ab.len(), bc.len());
        assert_eq!(ab.len(), ca.len());
        assert!(ab.len() >= 2);
        match self {
            None => *self = Self::one(ab, bc, points, calculate, shape),
            One(_) => self.three(ab, bc, ca, points, calculate, shape),
            Three { .. } => self.six(ab, bc, ca, points, calculate, shape),
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
                self.subdivide(ab, bc, ca, points, calculate, shape);
            }
            &mut More {
                a: a_idx,
                b: b_idx,
                c: c_idx,
                ref mut sides,
                ref mut contents,
                ref mut my_side_length,
            } => {
                points.extend_from_slice(&[V::ZERO, V::ZERO, V::ZERO]);
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
                    points[a_idx as usize] = shape.interpolate_half(aba, caa);
                    points[b_idx as usize] = shape.interpolate_half(abb, bcb);
                    points[c_idx as usize] = shape.interpolate_half(bcc, cac);
                }

                let ab = &sides[0..side_length];
                let bc = &sides[side_length..side_length * 2];
                let ca = &sides[side_length * 2..];

                if calculate {
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
                }

                contents.subdivide(Forward(ab), Forward(bc), Forward(ca), points, calculate, shape);
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
    pub const fn new(a: u32, b: u32, c: u32, ab_edge: usize, bc_edge: usize, ca_edge: usize) -> Self {
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

    fn subdivide_edges<V: Vec3>(
        &mut self,
        edges: &mut [Edge],
        points: &mut Vec<V>,
        calculate: bool,
        shape: &impl BaseShape<V>,
    ) -> usize {
        let mut divide = |p1: u32, p2: u32, edge_idx: usize, forward: &mut bool| {
            if !edges[edge_idx].done {
                edges[edge_idx].points.push(points.len() as u32);
                points.push(V::ZERO);

                if calculate {
                    shape.interpolate_multiple(
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

    fn subdivide<V: Vec3>(
        &mut self,
        edges: &mut [Edge],
        points: &mut Vec<V>,
        calculate: bool,
        shape: &impl BaseShape<V>,
    ) {
        let side_length = self.subdivide_edges(edges, points, calculate, shape) + 1;

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
            self.contents.subdivide(ab, bc, ca, points, calculate, shape);
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

pub struct Subdivided<T, V: Vec3, S: BaseShape<V>> {
    points: Vec<V>,
    data: Vec<T>,
    triangles: Box<[Triangle]>,
    shared_edges: Box<[Edge]>,
    subdivisions: usize,
    shape: S,
}

impl<T, V: Vec3, S: BaseShape<V> + Default> Subdivided<T, V, S> {
    pub fn new(subdivisions: usize, generator: impl FnMut(V) -> T) -> Self {
        Self::new_custom_shape(subdivisions, generator, S::default())
    }
}

impl<T, V: Vec3, S: BaseShape<V>> Subdivided<T, V, S> {
    pub fn new_custom_shape(subdivisions: usize, generator: impl FnMut(V) -> T, shape: S) -> Self {
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

    fn subdivide(&mut self, calculate: bool) {
        for Edge { done, .. } in &mut *self.shared_edges {
            *done = false;
        }

        for triangle in &mut *self.triangles {
            triangle.subdivide(&mut *self.shared_edges, &mut self.points, calculate, &self.shape);
        }
    }

    fn raw_user_data(&self) -> &T {
        &self.data
    }

    pub fn raw_points(&self) -> &[V] {
        &self.points
    }

    pub fn get_indices(&self, triangle: usize, buffer: &mut Vec<u32>) {
        self.triangles[triangle].add_indices(buffer, &self.shared_edges);
    }

    pub fn get_all_indices(&self) -> Vec<u32> {
        let mut buffer = Vec::new();

        for i in 0..self.triangles.len() {
            self.get_indices(i, &mut buffer);
        }

        buffer
    }

    pub fn subdivisions(&self) -> usize {
        self.subdivisions
    }
}

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

#[cfg(test)]
mod tests {
    use crate::shapes::IcoSphere;
    use crate::Slice::Forward;
    use glam::Vec3A;

    impl crate::Vec3 for glam::Vec3A {
        const ZERO: Self = Vec3A::ZERO;

        fn dot(self, other: Self) -> f32 {
            self.dot(other)
        }

        fn normalize(self) -> Self {
            self.normalize()
        }

        fn from_arr3([x, y, z]: [f32; 3]) -> Self {
            Self::new(x, y, z)
        }
    }

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
        let x = IcoSphere::new(0, |_: Vec3A| ());
        x.get_indices(0, &mut Vec::new());
    }

    #[test]
    fn one() {
        let x = IcoSphere::new(1, |_: Vec3A| ());
        x.get_indices(0, &mut Vec::new());
    }

    #[test]
    fn second_layer_inner() {
        let x = IcoSphere::new(2, |_: Vec3A| ());
        x.get_indices(0, &mut Vec::new());
        let x = IcoSphere::new(3, |_: Vec3A| ());
        x.get_indices(0, &mut Vec::new());
        let x = IcoSphere::new(5, |_: Vec3A| ());
        x.get_indices(0, &mut Vec::new());
        let x = IcoSphere::new(6, |_: Vec3A| ());
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
        let sphere = IcoSphere::new(10, |_: Vec3A| ());

        for i in sphere.raw_points() {
            assert!(i.length() - 1.0 <= EPSILON);
        }
    }
}
