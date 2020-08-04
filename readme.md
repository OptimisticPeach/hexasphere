# `hexasphere`

Library for tiling a sphere with hexagons and pentagons.

### Geometry:
This starts off with an icoashedron, and then proceeds to subdivide it.
The subdivided points form the dual to a [goldberg polyhedron](https://en.wikipedia.org/wiki/Goldberg_polyhedron). 
In essence, each point on this sphere is either a hexagon or a pentagon. 

### Interpolation:
Points are interpolated using the [geometric `slerp`](https://en.wikipedia.org/wiki/Slerp#Geometric_Slerp)
function to preserve their accuracy.

### Features
(and optimizations) will be coming as I work on the game that will use this library.

Some things I plan on doing are:

- Add support for finding the adjacent points to a point.
- Move to `smallvec` and friends, and store data more efficiently
  in general. 
