# Libigl CoMISo Integration

This repository provides a standalone integration of [CoMISo](https://www.graphics.rwth-aachen.de/software/comiso/) (Constrained Mixed-Integer Solver) with [Libigl](https://libigl.github.io), preserving the functionality that was removed from Libigl in [v2.6.0](https://github.com/libigl/libigl/releases/tag/v2.6.0) ([PR #2384](https://github.com/libigl/libigl/pull/2384)).

## Overview

CoMISo is a powerful solver for optimizing discrete quadratic energies subject to linear and integer constraints. It's particularly useful for geometry processing applications such as:

- **Mixed-Integer Quadrangulation (MIQ)** - Generate high-quality quad meshes
- **Frame Field Generation** - Compute smooth frame fields on surfaces
- **N-RoSy Field Design** - Create N-directional fields with singularities

This repository is based on the CoMISo integration from [Libigl v2.5.0](https://github.com/libigl/libigl/releases/tag/v2.5.0) with:
- Updated CMake build system using CPM (CMake Package Manager)
- Updated CoMISo and GMM++ dependencies
- Preserved Libigl tutorials 505 (MIQ) and 506 (FrameField)
- Modern C++14 compatibility

*Created by Qingnan Zhou for research purposes.*

## Features

### Available CoMISo Functions
- `igl::copyleft::comiso::miq` - Mixed-Integer Quadrangulation
- `igl::copyleft::comiso::nrosy` - N-directional field generation
- `igl::copyleft::comiso::frame_field` - Frame field computation

### Included Tutorials
- **Tutorial 505**: Mixed-Integer Quadrangulation (MIQ) demonstration
- **Tutorial 506**: Frame field generation and visualization

## Dependencies

This project automatically fetches and builds the following dependencies via CMake:

- **Libigl** - Core geometry processing library
- **CoMISo** - Constrained Mixed-Integer Solver
- **GMM++** - Generic Matrix Template Library
- **Eigen** - Linear algebra library
- **GLFW** - OpenGL windowing (for tutorials)

## Build Instructions

### Prerequisites
- CMake 3.14 or higher
- C++14 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- OpenGL development libraries (for tutorials)

### Basic Build
```bash
# Clone the repository
git clone https://github.com/qnzhou/libigl-comiso.git
cd libigl-comiso

# Create build directory
mkdir build && cd build

# Configure with tutorials enabled
cmake .. -DLIBIGL_BUILD_TUTORIALS=ON -DLIBIGL_GLFW=ON

# Build
make -j$(nproc)
```

## Usage

### Running Tutorials

After building with tutorials enabled:

```bash
# Run Mixed-Integer Quadrangulation tutorial
./505_MIQ

# Run Frame Field tutorial  
./506_FrameField
```

## Python Bindings

Python bindings are available for all CoMISo functions: `nrosy()`, `miq()`, and `frame_field()`.

### Build

```bash
# From the libigl-comiso folder
pip install .
```

### Usage Example

```python
import numpy as np
import igl_comiso

# Load triangle mesh (V: vertices, F: faces)
# V, F = load_your_mesh()

# Generate 4-RoSy field with constraints
b = np.array([0, 10, 20])  # constrained face indices
bc = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # directions
R, S = igl_comiso.nrosy(V, F, b, bc, N=4)

# Create quad mesh parameterization
UV, FUV = igl_comiso.miq(V, F, R, np.array([]), gradientSize=30.0)

print(f"Generated {FUV.shape[0]} quads from {F.shape[0]} triangles")
```


## License

This project inherits the licensing from its components:
- Libigl components: [MPL2](https://github.com/libigl/libigl/blob/main/LICENSE.MPL2) and [GPL](https://github.com/libigl/libigl/blob/main/LICENSE.GPL)
- CoMISo: [GPL v3](https://www.gnu.org/licenses/gpl-3.0.html)

## References

- [CoMISo: Constrained Mixed-Integer Solver](https://www.graphics.rwth-aachen.de/software/comiso/)
- [Libigl: A simple C++ geometry processing library](https://libigl.github.io)
- [Mixed-Integer Quadrangulation Paper](https://www.graphics.rwth-aachen.de/publication/03234/)

## Acknowledgments

- Original CoMISo developers at RWTH Aachen University
- Libigl development team
- Contributors to the geometry processing community
