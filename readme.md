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

### `no_std` support

`no_std` support can be enabled by compiling with `--no-default-features` to
disable `std` support and `--features libm` for math functions that are only
defined in `std`. For example:

```toml
[dependencies]
hexasphere = { version = "15", default-features = false, features = ["libm"] }
```

To support both `std` and `no_std` builds in project, you can use the following
in your `Cargo.toml`:

```toml
[features]
default = ["std"]

std = ["hexasphere/std"]
libm = ["hexasphere/libm"]

[dependencies]
hexasphere = { version = "15", default-features = false }
```

# License

Hexasphere is distributed under the terms of either the MIT license, or the Apache License (Version 2.0)
at the user's choice.

See the files named LICENSE-MIT and LICENSE-APACHE2 relative to the root directory of this
project for more details.
