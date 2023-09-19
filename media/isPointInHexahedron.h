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

#ifndef ISPOINTINHEXAHEDRON_H
#define ISPOINTINHEXAHEDRON_H

#include <stdbool.h>
// for C code call
bool isPointInHexahedron(float px, float py, float pz,
                         float *vx, float *vy, float *vz);

bool isPointInHexahedron_strict(float px, float py, float pz,
                               float *vx, float *vy, float *vz);

#endif
