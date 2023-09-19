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



#ifndef _MEDIA_READ_INTERFACE_FILE_
#define _MEDIA_READ_INTERFACE_FILE_
#define MAX_BUF_LEN 1024

#include <iostream>
#include <vector>
#include "media_utility.hpp"

FILE *gfopen(const char *filename, const char *mode);

void read_interface_file(
    const char *interface_file,
    inter_t *interfaces);

/* 
 * Just read the grid data within the given
 *  [Xmin, Xmax]\times[Ymin, Ymax] domain.
 */
void read_grid_file(
    const char *grid_file,
    // the given grid
    float Xmin, float Xmax,
    float Ymin, float Ymax,
    int &NL,
    std::vector<int> &NGz, // how many z-grid in each layer
    inter_t *interfaces);

// check whether the elevation[ng[i]-1] == elevation[ng[i]] 
int checkGridData(int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces, 
    const char *grid_file); 

void read_bin_file(
    const char *bin_file,
    float *var,
    int dimx, 
    int dimy, 
    int dimz,
    int *bin_start, 
    int *bin_end, 
    int *bin_size, 
    size_t bin_line, 
    size_t bin_slice);

#endif
