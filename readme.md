# `hexasphere`

Library for subdividing shapes, such as an icosahedron. This provides
abstractions for almost anything to be possible, however there are many
predetermined shapes.

## In the case of the icosphere:

### Geometry:
This starts off with an icoashedron, and then proceeds to subdivide it.
The subdivided points form the dual to a [goldberg polyhedron](https://en.wikipedia.org/wiki/Goldberg_polyhedron). 
In essence, each point on this sphere is either a hexagon or a pentagon. 

### Interpolation:
Points are interpolated using the [geometric `slerp`](https://en.wikipedia.org/wiki/Slerp#Geometric_Slerp)
function to preserve their accuracy.

---

### Features

- Base shapes
  - Icosahedron
  - Tetrahedron
  - Triangle
  - Square
  - Cube
- Interpolation functions
  - Spherical interpolation
  - Linear interpolation
  - Normalized linear interpolation
- Basic optimizations for `p == 0.5` in interpolation.
- An adjacency (neighbour) map which can be generated from the indices
that this library creates. 

A few more optimizations are coming, along with more options for generation:

- Remove nested triangle layers by placing them into their own `Vec`.
- Allow interpolation functions to have access to state. (This would
permit certain things such as sampling noise as part of the generation).
