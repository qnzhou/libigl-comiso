### Libigl CoMISo

This repo contains the [Libigl](https://libigl.github.io)'s
[CoMISo](https://www.graphics.rwth-aachen.de/software/comiso/) integration, which was removed in
Libigl [v2.6.0 release](https://github.com/libigl/libigl/releases/tag/v2.6.0)
([PR2384](https://github.com/libigl/libigl/pull/2384)). This repo is based on the CoMISo integration
from [Libigl v2.5.0](https://github.com/libigl/libigl/releases/tag/v2.5.0) with an updated cmake
setup and updated CoMISo and GMM++ dependency. Libigl tutorial 505 is also included in this repo.

This repo is created by Qingnan Zhou for pure research purposes.

#### Build and run

```sh
mkdir build
cd build
cmake .. -DLIBIGL_BUILD_TUTORIALS=On -DLIBIGL_GLFW=On
make
./tutorial_505
```
