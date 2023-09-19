/***********************************************************************
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version. 
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 ***************************************************************************/


#ifndef _MEDIA_BIN2MODEL_
#define _MEDIA_BIN2MODEL_

#include "media_geometry3d.hpp"
#include "media_read_file.hpp"

#ifdef __cplusplus
extern "C" {
#endif
int media_bin2model_el_iso(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    float xmin, float xmax,
    float ymin, float ymax,
    int grid_type,
    int *bin_order,    // eg, [2, 0, 1]=[z, x, y] 0:x, 1:y, 2:z
    int *bin_size,     // [ndim1, ndim2, ndim3],
    float  *bin_spacing,  // [dh1, dh2, dh3],
    float  *bin_origin,   // [h0_1, h0_2, h0_3],
    const char *bin_file_rho,
    const char *bin_file_vp,
    const char *bin_file_vs  );

#ifdef __cplusplus
}
#endif

void parameterization_bin_el_iso_loc(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    int grid_type,
    std::vector<float> &xvec, 
    std::vector<float> &yvec, 
    std::vector<float> &zvec,
    float *bin_rho,
    float *bin_vp,
    float *bin_vs ); 


#endif  // MEDIA_BIN2MODEL 
