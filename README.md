# FD3D_ModelPrepare
  Functions for medium discretization

## Repository
  [https://github.com/jianglq6/FD3D_ModelPrepare](https://github.com/jianglq6/FD3D_ModelPrepare)

## License
  [GNU GENERAL PUBLIC LICENSE v3](https://gnu.org/licenses/gpl-3.0.html)

## Authors
  Luqian Jiang (jlq0608@gmail.com)  
  Wei Zhang    (zhangwei@sustech.edu.cn)

## Directory layout
  * media/media_layer2model.cpp:  
    &ensp;functions of medium parameterization for layer-based velocity model 

  - media/media_grid2model.cpp:  
    &ensp;functions of medium parameterization for grid-based velocity model 

  + media/media_utility.cpp:  
    &ensp;the functions required for computing equivalent medium parameterization 

  - media/media_geometry3d.cpp:  
    &ensp;3D geometry-related functions  

## Availability and use of the program package

  For instructions on how to invoke these functions in seismic waveform simulation codes, please refer to https://github.com/zw-sustech/CGFD3D.
  
  If you use this code for your own research, please cite:
  
  * Jiang, L. and Zhang, W., 2021. TTI equivalent medium parametrization method for the seismic waveform modelling of heterogeneous media with coarse grids. Geophysical Journal International, 227(3), pp.2016-2043. DOI:[10.1093/gji/ggab310](https://doi.org/10.1093/gji/ggab310)
  * Jiang, L. and Zhang, W., 2024. Efficient implementation of equivalent medium parametrization in finite-difference seismic wave simulation methods. Geophysical Journal International, 239(1), pp. 675–693. DOI:[10.1093/gji/ggae286](https://doi.org/10.1093/gji/ggae286)
