/***************************************************************************
 *
 * This function is used for medium parameterization.
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    06/2021: Created by Luqian Jiang 
 *
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
 
#include <iostream>
#include <vector>
#include <string.h>
#include <cmath>
#include <set>
#include <numeric>
#include "media_geometry3d.hpp"
#include "media_layer2model.hpp"
#include "media_read_file.hpp"
#include "media_utility.hpp"

// for reporting error, media_type int -> string
std::map<int, std::string> md_map = create_md2str_map();

extern int edgeTable[256];

/*============================= for C call =================================*/
//---- 0. one component
int media_layer2model_onecmp(float *var3d,
                             const float *x3d, 
                             const float *y3d, 
                             const float *z3d, 
                             size_t nx,
                             size_t ny, 
                             size_t nz,
                             int grid_type, 
                             const char *in_var_file,
                             const char *average_method,
                             int myid) // loc, har, ari
{
    inter_t interfaces;

    /*Read interface file*/
    read_interface_file(in_var_file, &interfaces);
  
    if (interfaces.media_type != ONE_COMPONENT) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_layer2model_onecmp() only supports one_component,\n"\
                        "       please check the media_type in %s! \n", 
                md_map[interfaces.media_type].c_str(), in_var_file);   
        fflush(stderr);
        exit(1);
    }

    if (strcmp(average_method, "loc") == 0) {
        parametrization_layer_onecmp_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, interfaces, var3d, myid);
    } else if (strcmp(average_method, "har") == 0) {      // harmonic
        parametrization_layer_onecmp_har(nx, ny, nz, x3d, y3d, z3d, grid_type, interfaces, var3d, myid);
    } else if (strcmp(average_method, "ari") == 0) {      //arithemtic
        parametrization_layer_onecmp_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, interfaces, var3d, myid);
    } else {                                                        // default = loc
        fprintf(stderr, "Error: Wrong average method %s for one_component. \n", average_method);        
        fflush(stderr);
        exit(1);
    }


    return 0;
}

//---- 1. acoustic isotropic
int media_layer2model_ac_iso(
        float *rho3d,
        float *kappa3d,
        const float *x3d, 
        const float *y3d, 
        const float *z3d, 
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file,
        const char *equivalent_medium_method,
        int myid) 
{

    inter_t interfaces;

    /*Read interface file*/
    read_interface_file(in_3lay_file, &interfaces);

    if (interfaces.media_type != ACOUSTIC_ISOTROPIC) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_layer2model_ac_iso() only supports acoustic_isotropic,\n"\
                        "       please check the media_type in %s! \n", 
                md_map[interfaces.media_type].c_str(), in_3lay_file);         
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_layer_ac_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, kappa3d, rho3d, myid);
    } else if (strcmp(equivalent_medium_method, "ari") == 0 ) {
        parametrization_layer_ac_iso_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, kappa3d, rho3d, myid);
    } else if (strcmp(equivalent_medium_method, "har") == 0 ) {
        parametrization_layer_ac_iso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, kappa3d, rho3d, myid);
    } else { 
        fprintf(stderr, "Error: Wrong average method %s for acoustic_isotropic media. \n",
                 equivalent_medium_method);        
        fflush(stderr);
        exit(1);
    }


    return 0; 
}

//----  2. elastic isotropic
int media_layer2model_el_iso(
        float *lam3d,
        float *mu3d,
        float *rho3d,
        const float *x3d, 
        const float *y3d, 
        const float *z3d, 
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file,
        const char *equivalent_medium_method,
        int myid) // 
{

    inter_t interfaces;

    /* Read interface file */
    read_interface_file(in_3lay_file, &interfaces);   

    if (interfaces.media_type != ELASTIC_ISOTROPIC) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_layer2model_el_iso() only supports elastic_isotropic,\n"\
                        "       please check the media_type in %s! \n", 
                md_map[interfaces.media_type].c_str(), in_3lay_file);        
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_layer_el_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d, myid);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_layer_el_iso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d, myid);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_layer_el_iso_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d, myid);
    } else { //default
        fprintf(stderr,"Error: Wrong parametrization method %s in media_layer2model_el_iso(), " \
                       "       if you want to use tti equivalent medium method, " \
                       "       please call the media_layer2model_el_aniso().\n", 
                       equivalent_medium_method);        
        fflush(stderr);
        exit(1);
    }

    return 0; 
}

//--- 3. elastic vti
int media_layer2model_el_vti(
        float *rho,
        float *c11,
        float *c33,
        float *c55,
        float *c66,
        float *c13,
        const float *x3d,
        const float *y3d,
        const float *z3d,
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file, 
        const char *equivalent_medium_method,
        int myid) 
{

    inter_t interfaces;

    /* Read interface file */
    read_interface_file(in_3lay_file, &interfaces);   

    int md_type = interfaces.media_type;

    if (md_type != ELASTIC_VTI_PREM && md_type != ELASTIC_VTI_THOMSEN && md_type != ELASTIC_VTI_CIJ) {
        fprintf(stderr, "Error: media_type=%s is not supported in media_layer2model_el_vti(),\n"\
                        "       it only supports elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,\n"\
                        "       please check the media_type in %s! \n", 
                md_map[interfaces.media_type].c_str(), in_3lay_file);      
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_layer_el_vti_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c55, c66, c13, rho, myid);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_layer_el_vti_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c55, c66, c13, rho, myid);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_layer_el_vti_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c55, c66, c13, rho, myid);
    } else { //default
        fprintf(stderr,"Error: Wrong parametrization method %s in media_layer2model_el_vti(), " \
                       "       if you want to use tti equivalent medium method, " \
                       "       please call the media_layer2model_el_aniso().\n", 
                       equivalent_medium_method);        
        fflush(stderr);
        exit(1);
    }

    return 0; 
}

//--- 4. elastic anisotropic/TTI
int media_layer2model_el_aniso(
        float *rho,
        float *c11, float *c12, float *c13,
        float *c14, float *c15, float *c16,
        float *c22, float *c23, float *c24,
        float *c25, float *c26, float *c33,
        float *c34, float *c35, float *c36,
        float *c44, float *c45, float *c46,
        float *c55, float *c56, float *c66,
        const float *x3d,
        const float *y3d,
        const float *z3d,
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file,
        const char *equivalent_medium_method,
        int myid) 
{

    inter_t interfaces;
    size_t siz_volume = nx*ny*nz;

    /* Read interface file */
    read_interface_file(in_3lay_file, &interfaces);   

    int md_type = interfaces.media_type;

    // the function does not support one component and acoustic wave
    // if the 
    if (md_type == ONE_COMPONENT){
        fprintf(stderr, "Error: media_type=one_component is not supported in media_layer2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_3lay_file);        
        fflush(stderr);
        exit(1);
    } else if (md_type == ACOUSTIC_ISOTROPIC){
        fprintf(stderr, "Error: media_type=acoustic_isotropic is not supported in media_layer2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_3lay_file);        
        fflush(stderr);
        exit(1);
    }

    //- isotropic: loc, har, tti
    if (md_type == ELASTIC_ISOTROPIC) {
        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_layer_el_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c13, c44, rho, myid);
            for (size_t i = 0; i < siz_volume; ++i) {
                c11[i] = c13[i] + 2.0*c44[i]; 
                c22[i] = c11[i]; c33[i] = c11[i]; 
                c12[i] = c13[i]; c23[i] = c13[i];
                c55[i] = c44[i]; c66[i] = c44[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }
        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_layer_el_iso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c13, c44, rho, myid);
            for (size_t i = 0; i < siz_volume; ++i) {
                c11[i] = c13[i] + 2.0*c44[i]; 
                c22[i] = c11[i]; c33[i] = c11[i]; 
                c12[i] = c13[i]; c23[i] = c13[i];
                c55[i] = c44[i]; c66[i] = c44[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }
        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_layer_el_iso_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c13, c44, rho, myid);
            for (size_t i = 0; i < siz_volume; ++i) {
                c11[i] = c13[i] + 2.0*c44[i]; 
                c22[i] = c11[i]; c33[i] = c11[i]; 
                c12[i] = c13[i]; c23[i] = c13[i];
                c55[i] = c44[i]; c66[i] = c44[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }
        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
            parametrization_layer_el_iso_tti(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, rho, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, myid);
        } else {
            fprintf(stderr, "Error: no such equivalent_medium_method: %s!\n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    // vti: loc, har, ari    
    } else if (md_type == ELASTIC_VTI_PREM || md_type == ELASTIC_VTI_THOMSEN || md_type == ELASTIC_VTI_CIJ) {

        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_layer_el_vti_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c11, c33, c55, c66, c13, rho, myid);
            for (size_t i = 0; i < siz_volume; i++) {
                c12[i] = c11[i]-2.0*c66[i];
                c22[i] = c11[i];
                c23[i] = c13[i];
                c44[i] = c55[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }

        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_layer_el_vti_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c11, c33, c55, c66, c13, rho, myid);
            for (size_t i = 0; i < siz_volume; i++) {
                c12[i] = c11[i]-2.0*c66[i];
                c22[i] = c11[i];
                c23[i] = c13[i];
                c44[i] = c55[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }

        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_layer_el_vti_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c11, c33, c55, c66, c13, rho, myid);

            for (size_t i = 0; i < siz_volume; i++) {
                c12[i] = c11[i]-2.0*c66[i];
                c22[i] = c11[i];
                c23[i] = c13[i];
                c44[i] = c55[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }

        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
// TODO: tti equivalent medium for vti
        } else {
            fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    } else {
        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_layer_el_aniso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, myid);
        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_layer_el_aniso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, myid);
        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_layer_el_aniso_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, myid);
        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
// TODO: tti equivalent medium for tti
        } else { //default
            fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    }

    return 0; 
}

/*=================================================================================*/

int AssignLayerMediaPara2Point(
    size_t ix, size_t iy, size_t iz,         /* To print Error messages */ 
    Point3 A,  
    inter_t &interfaces,
    int media_type,                /* the type can be found in media_utility.hpp */ 
    std::vector<float> &var)
{
    size_t  NI = interfaces.NI;
    size_t  NX = interfaces.NX;
    size_t  NY = interfaces.NY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;
    int isPointOnInter = 0; 

    size_t interface_slice = NX*NY;

    /* If out of the INTERFACE MESH area, exit! */
    PrintIsPointOutOfInterfaceRange(A, ix, iy, iz, MINX, MAXX, MINY, MAXY);

    std::vector<float> XVEC(NX), YVEC(NY);
    std::vector<float> elevation(NI, -FLT_MAX);
    for (size_t i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (size_t i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }

    /* 
     * For each interface, interpolate to get the elevation of 
     *  point A at the projection position of the interface mesh.
     */
    for (int ni = 0; ni < NI; ni++) {
        elevation[ni] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*interface_slice, A.x, A.y);
    }

    /* Find which material_index to use */
    int mi = findLastGreaterEqualIndex(A.z, elevation);

    if (mi > -1 && isEqual(A.z, elevation[mi])) 
        isPointOnInter = 1;
    
    CalPointValue_layer(media_type, interfaces, interface_slice, XVEC, YVEC, A, elevation, mi, var);

    return isPointOnInter;
}

//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for layer2model
void CalPointValue_layer(int media_type, 
                   inter_t &interfaces,
                   size_t slice, 
                   std::vector<float> &xvec,  /* interface mesh */
                   std::vector<float> &yvec,
                   Point3 &A,
                   std::vector<float> &elevation, /*the elevation of point A at the projection position of the interface mesh. */
                   int mi,
                   std::vector<float> &var)
{
    float dz = elevation[mi] - A.z;

    switch(media_type)
    {
    case ONE_COMPONENT: /* 0. var */
        /* If grid_z > elevation of top_interface, it given by the medium of top non-zero thickness layer */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.var + mi*slice, A.x, A.y); 
        } else {
            float var0      = BilinearInterpolation(xvec, yvec, interfaces.var      + mi*slice , A.x, A.y);
            float var_grad  = BilinearInterpolation(xvec, yvec, interfaces.var_grad + mi*slice , A.x, A.y);
            float var_pow   = BilinearInterpolation(xvec, yvec, interfaces.var_pow  + mi*slice , A.x, A.y);
            var[0]  = var0  + pow(dz, var_pow)* var_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    var0      = BilinearInterpolation(xvec, yvec, interfaces.var      + mi*slice , A.x, A.y);
                    var_grad  = BilinearInterpolation(xvec, yvec, interfaces.var_grad + mi*slice , A.x, A.y);
                    var_pow   = BilinearInterpolation(xvec, yvec, interfaces.var_pow  + mi*slice , A.x, A.y);
                    dz = elevation[mi] - A.z;
                    float var_top = var0  + pow(dz, var_pow)* var_grad;
                    var[0] = (var[0] + var_top)/2.0;                    
                }
            }
        }
    break;

    case ELASTIC_ISOTROPIC: /*1. rho, vp, vs*/
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs  + mi*slice , A.x, A.y);
        } else {
            float vp       = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice, A.x, A.y);
            float vp_grad  = BilinearInterpolation(xvec, yvec, interfaces.vp_grad + mi*slice, A.x, A.y);
            float vp_pow   = BilinearInterpolation(xvec, yvec, interfaces.vp_pow  + mi*slice, A.x, A.y);
            float vs       = BilinearInterpolation(xvec, yvec, interfaces.vs      + mi*slice, A.x, A.y);
            float vs_grad  = BilinearInterpolation(xvec, yvec, interfaces.vs_grad + mi*slice, A.x, A.y);
            float vs_pow   = BilinearInterpolation(xvec, yvec, interfaces.vs_pow  + mi*slice, A.x, A.y);
            float rho      = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad+ mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice, A.x, A.y);
            var[1] = vp  + pow(dz, vp_pow)* vp_grad;
            var[2] = vs  + pow(dz, vs_pow)* vs_grad;
            var[0] = rho + pow(dz,rho_pow)*rho_grad;
            if (var[1] <= var[2] * sqrt(2)) {
              fprintf(stderr, "Error: vp must larger than sqrt(2)*vs!\n");
              fprintf(stderr, "       please check the input md3lay file!\n");
              fflush(stderr);  exit(1);
            }
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    vp       = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice, A.x, A.y);
                    vp_grad  = BilinearInterpolation(xvec, yvec, interfaces.vp_grad + mi*slice, A.x, A.y);
                    vp_pow   = BilinearInterpolation(xvec, yvec, interfaces.vp_pow  + mi*slice, A.x, A.y);
                    vs       = BilinearInterpolation(xvec, yvec, interfaces.vs      + mi*slice, A.x, A.y);
                    vs_grad  = BilinearInterpolation(xvec, yvec, interfaces.vs_grad + mi*slice, A.x, A.y);
                    vs_pow   = BilinearInterpolation(xvec, yvec, interfaces.vs_pow  + mi*slice, A.x, A.y);
                    rho      = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
                    rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad+ mi*slice, A.x, A.y);
                    rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice, A.x, A.y);
                    dz = elevation[mi] - A.z;
                    float vp_top  = vp  + pow(dz,  vp_pow)* vp_grad;
                    float vs_top  = vs  + pow(dz,  vs_pow)* vs_grad;
                    float rho_top = rho + pow(dz, rho_pow)* rho_grad;
                    var[0] = (var[0] + rho_top)/2.0;
                    var[1] = (var[1] +  vp_top)/2.0;
                    var[2] = (var[2] +  vs_top)/2.0;
                }
            }
        }
    break;

    case ELASTIC_VTI_PREM: /*2. rho, vph, vpv, vsh, vsv, eta */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);     
        } else {
            float vph = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            float vpv = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            float vsh = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            float vsv = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float eta = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);
            float vph_grad = BilinearInterpolation(xvec, yvec, interfaces.vph_grad + mi*slice, A.x, A.y);
            float vpv_grad = BilinearInterpolation(xvec, yvec, interfaces.vpv_grad + mi*slice, A.x, A.y);
            float vsh_grad = BilinearInterpolation(xvec, yvec, interfaces.vsh_grad + mi*slice, A.x, A.y);
            float vsv_grad = BilinearInterpolation(xvec, yvec, interfaces.vsv_grad + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
            float eta_grad = BilinearInterpolation(xvec, yvec, interfaces.eta_grad + mi*slice, A.x, A.y);
            float vph_pow  = BilinearInterpolation(xvec, yvec, interfaces.vph_pow  + mi*slice, A.x, A.y);
            float vpv_pow  = BilinearInterpolation(xvec, yvec, interfaces.vpv_pow  + mi*slice, A.x, A.y);
            float vsh_pow  = BilinearInterpolation(xvec, yvec, interfaces.vsh_pow  + mi*slice, A.x, A.y);
            float vsv_pow  = BilinearInterpolation(xvec, yvec, interfaces.vsv_pow  + mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
            float eta_pow  = BilinearInterpolation(xvec, yvec, interfaces.eta_pow  + mi*slice, A.x, A.y);
            var[0] = rho + pow(dz, rho_pow)* rho_grad;
            var[1] = vph + pow(dz, vph_pow)* vph_grad;
            var[2] = vpv + pow(dz, vpv_pow)* vpv_grad;
            var[3] = vsh + pow(dz, vsh_pow)* vsh_grad;
            var[4] = vsv + pow(dz, vsv_pow)* vsv_grad;
            var[5] = eta + pow(dz, eta_pow)* eta_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    vph = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
                    vpv = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
                    vsh = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
                    vsv = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
                    rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
                    eta = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);
                    vph_grad = BilinearInterpolation(xvec, yvec, interfaces.vph_grad + mi*slice, A.x, A.y);
                    vpv_grad = BilinearInterpolation(xvec, yvec, interfaces.vpv_grad + mi*slice, A.x, A.y);
                    vsh_grad = BilinearInterpolation(xvec, yvec, interfaces.vsh_grad + mi*slice, A.x, A.y);
                    vsv_grad = BilinearInterpolation(xvec, yvec, interfaces.vsv_grad + mi*slice, A.x, A.y);
                    rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
                    eta_grad = BilinearInterpolation(xvec, yvec, interfaces.eta_grad + mi*slice, A.x, A.y);
                    vph_pow  = BilinearInterpolation(xvec, yvec, interfaces.vph_pow  + mi*slice, A.x, A.y);
                    vpv_pow  = BilinearInterpolation(xvec, yvec, interfaces.vpv_pow  + mi*slice, A.x, A.y);
                    vsh_pow  = BilinearInterpolation(xvec, yvec, interfaces.vsh_pow  + mi*slice, A.x, A.y);
                    vsv_pow  = BilinearInterpolation(xvec, yvec, interfaces.vsv_pow  + mi*slice, A.x, A.y);
                    rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
                    eta_pow  = BilinearInterpolation(xvec, yvec, interfaces.eta_pow  + mi*slice, A.x, A.y);
           
                    dz = elevation[mi] - A.z;
                    float rho_top = rho + pow(dz, rho_pow)* rho_grad;
                    float vph_top = vph + pow(dz, vph_pow)* vph_grad;
                    float vpv_top = vpv + pow(dz, vpv_pow)* vpv_grad;
                    float vsv_top = vsv + pow(dz, vsv_pow)* vsv_grad;
                    float vsh_top = vsh + pow(dz, vsh_pow)* vsh_grad;
                    float eta_top = eta + pow(dz, eta_pow)* eta_grad;
                    var[0] = (rho_top + var[0])/2.0;
                    var[1] = (vph_top + var[1])/2.0;
                    var[2] = (vpv_top + var[2])/2.0;
                    var[3] = (vsh_top + var[3])/2.0;
                    var[4] = (vsv_top + var[4])/2.0;
                    var[5] = (eta_top + var[5])/2.0;
                }
            }
        }   
    break;

    case ELASTIC_VTI_THOMSEN: /*3. rho, vp0, vs0, epsilon, delta, gamma */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho    + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp0    + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs0    + mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.epsilon+ mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.delta  + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.gamma  + mi*slice , A.x, A.y);    
        } else {
            float rho     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float vp0     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
            float vs0     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
            float epsil   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);  // epsilon
            float delta   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
            float gamma   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
            float rho_grad     = BilinearInterpolation(xvec, yvec, interfaces.rho_grad     + mi*slice, A.x, A.y);
            float vp0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vp0_grad     + mi*slice, A.x, A.y);
            float vs0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vs0_grad     + mi*slice, A.x, A.y);
            float epsilon_grad = BilinearInterpolation(xvec, yvec, interfaces.epsilon_grad + mi*slice, A.x, A.y);
            float delta_grad   = BilinearInterpolation(xvec, yvec, interfaces.delta_grad   + mi*slice, A.x, A.y);
            float gamma_grad   = BilinearInterpolation(xvec, yvec, interfaces.gamma_grad   + mi*slice, A.x, A.y);
            float rho_pow      = BilinearInterpolation(xvec, yvec, interfaces.rho_pow      + mi*slice, A.x, A.y);
            float vp0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vp0_pow      + mi*slice, A.x, A.y);
            float vs0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vs0_pow      + mi*slice, A.x, A.y);
            float epsilon_pow  = BilinearInterpolation(xvec, yvec, interfaces.epsilon_pow  + mi*slice, A.x, A.y);
            float delta_pow    = BilinearInterpolation(xvec, yvec, interfaces.delta_pow    + mi*slice, A.x, A.y);
            float gamma_pow    = BilinearInterpolation(xvec, yvec, interfaces.gamma_pow    + mi*slice, A.x, A.y);
            var[0] = rho     + pow(dz, rho_pow)    * rho_grad;
            var[1] = vp0     + pow(dz, vp0_pow)    * vp0_grad;
            var[2] = vs0     + pow(dz, vs0_pow)    * vs0_grad;
            var[3] = epsil   + pow(dz, epsilon_pow)* epsilon_grad;
            var[4] = delta   + pow(dz, delta_pow)  * delta_grad;
            var[5] = gamma   + pow(dz, gamma_pow)  * gamma_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    rho     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
                    vp0     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
                    vs0     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
                    epsil   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);  // epsilon
                    delta   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
                    gamma   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
                    rho_grad     = BilinearInterpolation(xvec, yvec, interfaces.rho_grad     + mi*slice, A.x, A.y);
                    vp0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vp0_grad     + mi*slice, A.x, A.y);
                    vs0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vs0_grad     + mi*slice, A.x, A.y);
                    epsilon_grad = BilinearInterpolation(xvec, yvec, interfaces.epsilon_grad + mi*slice, A.x, A.y);
                    delta_grad   = BilinearInterpolation(xvec, yvec, interfaces.delta_grad   + mi*slice, A.x, A.y);
                    gamma_grad   = BilinearInterpolation(xvec, yvec, interfaces.gamma_grad   + mi*slice, A.x, A.y);
                    rho_pow      = BilinearInterpolation(xvec, yvec, interfaces.rho_pow      + mi*slice, A.x, A.y);
                    vp0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vp0_pow      + mi*slice, A.x, A.y);
                    vs0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vs0_pow      + mi*slice, A.x, A.y);
                    epsilon_pow  = BilinearInterpolation(xvec, yvec, interfaces.epsilon_pow  + mi*slice, A.x, A.y);
                    delta_pow    = BilinearInterpolation(xvec, yvec, interfaces.delta_pow    + mi*slice, A.x, A.y);
                    gamma_pow    = BilinearInterpolation(xvec, yvec, interfaces.gamma_pow    + mi*slice, A.x, A.y);
                    dz = elevation[mi] - A.z;
                    float var0_top = rho     + pow(dz, rho_pow)    * rho_grad;
                    float var1_top = vp0     + pow(dz, vp0_pow)    * vp0_grad;
                    float var2_top = vs0     + pow(dz, vs0_pow)    * vs0_grad;
                    float var3_top = epsil   + pow(dz, epsilon_pow)* epsilon_grad;
                    float var4_top = delta   + pow(dz, delta_pow)  * delta_grad;
                    float var5_top = gamma   + pow(dz, gamma_pow)  * gamma_grad;
                    var[0] = (var[0] + var0_top) / 2.0;
                    var[1] = (var[1] + var1_top) / 2.0;
                    var[2] = (var[2] + var2_top) / 2.0;
                    var[3] = (var[3] + var3_top) / 2.0;
                    var[4] = (var[4] + var4_top) / 2.0;
                    var[5] = (var[5] + var5_top) / 2.0;
                }
            }
        }   
    break;

    case ELASTIC_VTI_CIJ: /*4. rho c11 c33 c55 c66 c13 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);   
        } else {
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c55 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice, A.x, A.y);
            float c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice, A.x, A.y);
            float c55_grad = BilinearInterpolation(xvec, yvec, interfaces.c55_grad + mi*slice, A.x, A.y);
            float c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice, A.x, A.y);
            float c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
            float c11_pow  = BilinearInterpolation(xvec, yvec, interfaces.c11_pow  + mi*slice, A.x, A.y);
            float c33_pow  = BilinearInterpolation(xvec, yvec, interfaces.c33_pow  + mi*slice, A.x, A.y);
            float c55_pow  = BilinearInterpolation(xvec, yvec, interfaces.c55_pow  + mi*slice, A.x, A.y);
            float c66_pow  = BilinearInterpolation(xvec, yvec, interfaces.c66_pow  + mi*slice, A.x, A.y);
            float c13_pow  = BilinearInterpolation(xvec, yvec, interfaces.c13_pow  + mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
            var[0] = rho + pow(dz, rho_pow)* rho_grad;
            var[1] = c11 + pow(dz, c11_pow)* c11_grad;
            var[2] = c33 + pow(dz, c33_pow)* c33_grad;
            var[3] = c55 + pow(dz, c55_pow)* c55_grad;
            var[4] = c66 + pow(dz, c66_pow)* c66_grad;
            var[5] = c13 + pow(dz, c13_pow)* c13_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
                    c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
                    c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
                    c55 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
                    c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
                    c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
                    c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice, A.x, A.y);
                    c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice, A.x, A.y);
                    c55_grad = BilinearInterpolation(xvec, yvec, interfaces.c55_grad + mi*slice, A.x, A.y);
                    c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice, A.x, A.y);
                    c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice, A.x, A.y);
                    rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
                    c11_pow  = BilinearInterpolation(xvec, yvec, interfaces.c11_pow  + mi*slice, A.x, A.y);
                    c33_pow  = BilinearInterpolation(xvec, yvec, interfaces.c33_pow  + mi*slice, A.x, A.y);
                    c55_pow  = BilinearInterpolation(xvec, yvec, interfaces.c55_pow  + mi*slice, A.x, A.y);
                    c66_pow  = BilinearInterpolation(xvec, yvec, interfaces.c66_pow  + mi*slice, A.x, A.y);
                    c13_pow  = BilinearInterpolation(xvec, yvec, interfaces.c13_pow  + mi*slice, A.x, A.y);
                    rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
                    float var0_top = rho + pow(dz, rho_pow)* rho_grad;
                    float var1_top = c11 + pow(dz, c11_pow)* c11_grad;
                    float var2_top = c33 + pow(dz, c33_pow)* c33_grad;
                    float var3_top = c55 + pow(dz, c55_pow)* c55_grad;
                    float var4_top = c66 + pow(dz, c66_pow)* c66_grad;
                    float var5_top = c13 + pow(dz, c13_pow)* c13_grad;
                    var[0] = (var[0] + var0_top) / 2.0;
                    var[1] = (var[1] + var1_top) / 2.0;
                    var[2] = (var[2] + var2_top) / 2.0;
                    var[3] = (var[3] + var3_top) / 2.0;
                    var[4] = (var[4] + var4_top) / 2.0;
                    var[5] = (var[5] + var5_top) / 2.0;

                }
            }
        }   
    break;

    case ELASTIC_TTI_THOMSEN: /*5. rho, vp0, vs0, epsilon, delta, gamma, azimuth, dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice , A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice , A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice , A.x, A.y);
        } else {
            float vp0     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
            float vs0     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
            float epsil   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);
            float delta   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
            float gamma   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
            float rho     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float azimuth = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
            float dip     = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
            float vp0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vp0_grad     + mi*slice, A.x, A.y);
            float vs0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vs0_grad     + mi*slice, A.x, A.y);
            float epsilon_grad = BilinearInterpolation(xvec, yvec, interfaces.epsilon_grad + mi*slice, A.x, A.y);
            float delta_grad   = BilinearInterpolation(xvec, yvec, interfaces.delta_grad   + mi*slice, A.x, A.y);
            float gamma_grad   = BilinearInterpolation(xvec, yvec, interfaces.gamma_grad   + mi*slice, A.x, A.y);
            float rho_grad     = BilinearInterpolation(xvec, yvec, interfaces.rho_grad     + mi*slice, A.x, A.y);
            float azimuth_grad = BilinearInterpolation(xvec, yvec, interfaces.azimuth_grad + mi*slice, A.x, A.y);
            float dip_grad     = BilinearInterpolation(xvec, yvec, interfaces.dip_grad     + mi*slice, A.x, A.y);
            float vp0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vp0_pow      + mi*slice, A.x, A.y);
            float vs0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vs0_pow      + mi*slice, A.x, A.y);
            float epsilon_pow  = BilinearInterpolation(xvec, yvec, interfaces.epsilon_pow  + mi*slice, A.x, A.y);
            float delta_pow    = BilinearInterpolation(xvec, yvec, interfaces.delta_pow    + mi*slice, A.x, A.y);
            float gamma_pow    = BilinearInterpolation(xvec, yvec, interfaces.gamma_pow    + mi*slice, A.x, A.y);
            float rho_pow      = BilinearInterpolation(xvec, yvec, interfaces.rho_pow      + mi*slice, A.x, A.y);
            float azimuth_pow  = BilinearInterpolation(xvec, yvec, interfaces.azimuth_pow  + mi*slice, A.x, A.y);
            float dip_pow      = BilinearInterpolation(xvec, yvec, interfaces.dip_pow      + mi*slice, A.x, A.y);

            var[0] = rho     + pow(dz, rho_pow)    * rho_grad;
            var[1] = vp0     + pow(dz, vp0_pow)    * vp0_grad;
            var[2] = vs0     + pow(dz, vs0_pow)    * vs0_grad;
            var[3] = epsil   + pow(dz, epsilon_pow)* epsilon_grad;
            var[4] = delta   + pow(dz, delta_pow)  * delta_grad;
            var[5] = gamma   + pow(dz, gamma_pow)  * gamma_grad;
            var[6] = azimuth + pow(dz, azimuth_pow)* azimuth_grad;
            var[7] = dip     + pow(dz, dip_pow)    * dip_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    vp0     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
                    vs0     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
                    epsil   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);
                    delta   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
                    gamma   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
                    rho     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
                    azimuth = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
                    dip     = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
                    vp0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vp0_grad     + mi*slice, A.x, A.y);
                    vs0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vs0_grad     + mi*slice, A.x, A.y);
                    epsilon_grad = BilinearInterpolation(xvec, yvec, interfaces.epsilon_grad + mi*slice, A.x, A.y);
                    delta_grad   = BilinearInterpolation(xvec, yvec, interfaces.delta_grad   + mi*slice, A.x, A.y);
                    gamma_grad   = BilinearInterpolation(xvec, yvec, interfaces.gamma_grad   + mi*slice, A.x, A.y);
                    rho_grad     = BilinearInterpolation(xvec, yvec, interfaces.rho_grad     + mi*slice, A.x, A.y);
                    azimuth_grad = BilinearInterpolation(xvec, yvec, interfaces.azimuth_grad + mi*slice, A.x, A.y);
                    dip_grad     = BilinearInterpolation(xvec, yvec, interfaces.dip_grad     + mi*slice, A.x, A.y);
                    vp0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vp0_pow      + mi*slice, A.x, A.y);
                    vs0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vs0_pow      + mi*slice, A.x, A.y);
                    epsilon_pow  = BilinearInterpolation(xvec, yvec, interfaces.epsilon_pow  + mi*slice, A.x, A.y);
                    delta_pow    = BilinearInterpolation(xvec, yvec, interfaces.delta_pow    + mi*slice, A.x, A.y);
                    gamma_pow    = BilinearInterpolation(xvec, yvec, interfaces.gamma_pow    + mi*slice, A.x, A.y);
                    rho_pow      = BilinearInterpolation(xvec, yvec, interfaces.rho_pow      + mi*slice, A.x, A.y);
                    azimuth_pow  = BilinearInterpolation(xvec, yvec, interfaces.azimuth_pow  + mi*slice, A.x, A.y);
                    dip_pow      = BilinearInterpolation(xvec, yvec, interfaces.dip_pow      + mi*slice, A.x, A.y);
                    float var0_top = rho     + pow(dz, rho_pow)    * rho_grad;
                    float var1_top = vp0     + pow(dz, vp0_pow)    * vp0_grad;
                    float var2_top = vs0     + pow(dz, vs0_pow)    * vs0_grad;
                    float var3_top = epsil   + pow(dz, epsilon_pow)* epsilon_grad;
                    float var4_top = delta   + pow(dz, delta_pow)  * delta_grad;
                    float var5_top = gamma   + pow(dz, gamma_pow)  * gamma_grad;
                    float var6_top = azimuth + pow(dz, azimuth_pow)* azimuth_grad;
                    float var7_top = dip     + pow(dz, dip_pow)    * dip_grad;
                    var[0] = (var[0] + var0_top)/2.0; 
                    var[1] = (var[1] + var1_top)/2.0; 
                    var[2] = (var[2] + var2_top)/2.0; 
                    var[3] = (var[3] + var3_top)/2.0; 
                    var[4] = (var[4] + var4_top)/2.0; 
                    var[5] = (var[5] + var5_top)/2.0; 
                    var[6] = (var[6] + var6_top)/2.0; 
                    var[7] = (var[7] + var7_top)/2.0; 
                }
            }
        }   
    break;

    case ELASTIC_ANISO_CIJ: /* 7. rho c11 c12 c13 c14 c15 c16 c22 c23 c24 c25 c26 c33 c34 c35 c36 c44 c45 c46 c55 c56 c66 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0]  = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1]  = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2]  = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            var[3]  = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[4]  = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            var[5]  = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            var[6]  = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            var[7]  = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            var[8]  = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            var[9]  = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            var[10] = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            var[11] = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            var[12] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[13] = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            var[14] = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            var[15] = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            var[16] = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            var[17] = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            var[18] = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            var[19] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[20] = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            var[21] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
        } else {
            float c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c12 = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            float c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float c14 = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            float c15 = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            float c16 = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            float c22 = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            float c23 = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            float c24 = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            float c25 = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            float c26 = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            float c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c34 = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            float c35 = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            float c36 = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            float c44 = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            float c45 = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            float c46 = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            float c55 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c56 = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            float c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11_pow = BilinearInterpolation(xvec, yvec, interfaces.c11_pow + mi*slice , A.x, A.y);
            float c12_pow = BilinearInterpolation(xvec, yvec, interfaces.c12_pow + mi*slice , A.x, A.y);
            float c13_pow = BilinearInterpolation(xvec, yvec, interfaces.c13_pow + mi*slice , A.x, A.y);
            float c14_pow = BilinearInterpolation(xvec, yvec, interfaces.c14_pow + mi*slice , A.x, A.y);
            float c15_pow = BilinearInterpolation(xvec, yvec, interfaces.c15_pow + mi*slice , A.x, A.y);
            float c16_pow = BilinearInterpolation(xvec, yvec, interfaces.c16_pow + mi*slice , A.x, A.y);
            float c22_pow = BilinearInterpolation(xvec, yvec, interfaces.c22_pow + mi*slice , A.x, A.y);
            float c23_pow = BilinearInterpolation(xvec, yvec, interfaces.c23_pow + mi*slice , A.x, A.y);
            float c24_pow = BilinearInterpolation(xvec, yvec, interfaces.c24_pow + mi*slice , A.x, A.y);
            float c25_pow = BilinearInterpolation(xvec, yvec, interfaces.c25_pow + mi*slice , A.x, A.y);
            float c26_pow = BilinearInterpolation(xvec, yvec, interfaces.c26_pow + mi*slice , A.x, A.y);
            float c33_pow = BilinearInterpolation(xvec, yvec, interfaces.c33_pow + mi*slice , A.x, A.y);
            float c34_pow = BilinearInterpolation(xvec, yvec, interfaces.c34_pow + mi*slice , A.x, A.y);
            float c35_pow = BilinearInterpolation(xvec, yvec, interfaces.c35_pow + mi*slice , A.x, A.y);
            float c36_pow = BilinearInterpolation(xvec, yvec, interfaces.c36_pow + mi*slice , A.x, A.y);
            float c44_pow = BilinearInterpolation(xvec, yvec, interfaces.c44_pow + mi*slice , A.x, A.y);
            float c45_pow = BilinearInterpolation(xvec, yvec, interfaces.c45_pow + mi*slice , A.x, A.y);
            float c46_pow = BilinearInterpolation(xvec, yvec, interfaces.c46_pow + mi*slice , A.x, A.y);
            float c55_pow = BilinearInterpolation(xvec, yvec, interfaces.c55_pow + mi*slice , A.x, A.y);
            float c56_pow = BilinearInterpolation(xvec, yvec, interfaces.c56_pow + mi*slice , A.x, A.y);
            float c66_pow = BilinearInterpolation(xvec, yvec, interfaces.c66_pow + mi*slice , A.x, A.y);
            float rho_pow = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice , A.x, A.y);
            float c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice , A.x, A.y);
            float c12_grad = BilinearInterpolation(xvec, yvec, interfaces.c12_grad + mi*slice , A.x, A.y);
            float c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice , A.x, A.y);
            float c14_grad = BilinearInterpolation(xvec, yvec, interfaces.c14_grad + mi*slice , A.x, A.y);
            float c15_grad = BilinearInterpolation(xvec, yvec, interfaces.c15_grad + mi*slice , A.x, A.y);
            float c16_grad = BilinearInterpolation(xvec, yvec, interfaces.c16_grad + mi*slice , A.x, A.y);
            float c22_grad = BilinearInterpolation(xvec, yvec, interfaces.c22_grad + mi*slice , A.x, A.y);
            float c23_grad = BilinearInterpolation(xvec, yvec, interfaces.c23_grad + mi*slice , A.x, A.y);
            float c24_grad = BilinearInterpolation(xvec, yvec, interfaces.c24_grad + mi*slice , A.x, A.y);
            float c25_grad = BilinearInterpolation(xvec, yvec, interfaces.c25_grad + mi*slice , A.x, A.y);
            float c26_grad = BilinearInterpolation(xvec, yvec, interfaces.c26_grad + mi*slice , A.x, A.y);
            float c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice , A.x, A.y);
            float c34_grad = BilinearInterpolation(xvec, yvec, interfaces.c34_grad + mi*slice , A.x, A.y);
            float c35_grad = BilinearInterpolation(xvec, yvec, interfaces.c35_grad + mi*slice , A.x, A.y);
            float c36_grad = BilinearInterpolation(xvec, yvec, interfaces.c36_grad + mi*slice , A.x, A.y);
            float c44_grad = BilinearInterpolation(xvec, yvec, interfaces.c44_grad + mi*slice , A.x, A.y);
            float c45_grad = BilinearInterpolation(xvec, yvec, interfaces.c45_grad + mi*slice , A.x, A.y);
            float c46_grad = BilinearInterpolation(xvec, yvec, interfaces.c46_grad + mi*slice , A.x, A.y);
            float c55_grad = BilinearInterpolation(xvec, yvec, interfaces.c55_grad + mi*slice , A.x, A.y);
            float c56_grad = BilinearInterpolation(xvec, yvec, interfaces.c56_grad + mi*slice , A.x, A.y);
            float c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice , A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice , A.x, A.y);
            var[0]  = rho + pow(dz, rho_pow) * rho_grad;
            var[1]  = c11 + pow(dz, c11_pow) * c11_grad; 
            var[2]  = c12 + pow(dz, c12_pow) * c12_grad;
            var[3]  = c13 + pow(dz, c13_pow) * c13_grad;
            var[4]  = c14 + pow(dz, c14_pow) * c14_grad;
            var[5]  = c15 + pow(dz, c15_pow) * c15_grad;
            var[6]  = c16 + pow(dz, c16_pow) * c16_grad;
            var[7]  = c22 + pow(dz, c22_pow) * c22_grad;
            var[8]  = c23 + pow(dz, c23_pow) * c23_grad;
            var[9]  = c24 + pow(dz, c24_pow) * c24_grad;
            var[10] = c25 + pow(dz, c25_pow) * c25_grad;
            var[11] = c26 + pow(dz, c26_pow) * c26_grad;
            var[12] = c33 + pow(dz, c33_pow) * c33_grad;
            var[13] = c34 + pow(dz, c34_pow) * c34_grad;
            var[14] = c35 + pow(dz, c35_pow) * c35_grad;
            var[15] = c36 + pow(dz, c36_pow) * c36_grad;
            var[16] = c44 + pow(dz, c44_pow) * c44_grad;
            var[17] = c45 + pow(dz, c45_pow) * c45_grad;
            var[18] = c46 + pow(dz, c46_pow) * c46_grad;
            var[19] = c55 + pow(dz, c55_pow) * c55_grad;
            var[20] = c56 + pow(dz, c56_pow) * c56_grad;
            var[21] = c66 + pow(dz, c66_pow) * c66_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
                    c12 = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
                    c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
                    c14 = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
                    c15 = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
                    c16 = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
                    c22 = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
                    c23 = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
                    c24 = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
                    c25 = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
                    c26 = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
                    c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
                    c34 = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
                    c35 = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
                    c36 = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
                    c44 = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
                    c45 = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
                    c46 = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
                    c55 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
                    c56 = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
                    c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
                    rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
                    c11_pow = BilinearInterpolation(xvec, yvec, interfaces.c11_pow + mi*slice , A.x, A.y);
                    c12_pow = BilinearInterpolation(xvec, yvec, interfaces.c12_pow + mi*slice , A.x, A.y);
                    c13_pow = BilinearInterpolation(xvec, yvec, interfaces.c13_pow + mi*slice , A.x, A.y);
                    c14_pow = BilinearInterpolation(xvec, yvec, interfaces.c14_pow + mi*slice , A.x, A.y);
                    c15_pow = BilinearInterpolation(xvec, yvec, interfaces.c15_pow + mi*slice , A.x, A.y);
                    c16_pow = BilinearInterpolation(xvec, yvec, interfaces.c16_pow + mi*slice , A.x, A.y);
                    c22_pow = BilinearInterpolation(xvec, yvec, interfaces.c22_pow + mi*slice , A.x, A.y);
                    c23_pow = BilinearInterpolation(xvec, yvec, interfaces.c23_pow + mi*slice , A.x, A.y);
                    c24_pow = BilinearInterpolation(xvec, yvec, interfaces.c24_pow + mi*slice , A.x, A.y);
                    c25_pow = BilinearInterpolation(xvec, yvec, interfaces.c25_pow + mi*slice , A.x, A.y);
                    c26_pow = BilinearInterpolation(xvec, yvec, interfaces.c26_pow + mi*slice , A.x, A.y);
                    c33_pow = BilinearInterpolation(xvec, yvec, interfaces.c33_pow + mi*slice , A.x, A.y);
                    c34_pow = BilinearInterpolation(xvec, yvec, interfaces.c34_pow + mi*slice , A.x, A.y);
                    c35_pow = BilinearInterpolation(xvec, yvec, interfaces.c35_pow + mi*slice , A.x, A.y);
                    c36_pow = BilinearInterpolation(xvec, yvec, interfaces.c36_pow + mi*slice , A.x, A.y);
                    c44_pow = BilinearInterpolation(xvec, yvec, interfaces.c44_pow + mi*slice , A.x, A.y);
                    c45_pow = BilinearInterpolation(xvec, yvec, interfaces.c45_pow + mi*slice , A.x, A.y);
                    c46_pow = BilinearInterpolation(xvec, yvec, interfaces.c46_pow + mi*slice , A.x, A.y);
                    c55_pow = BilinearInterpolation(xvec, yvec, interfaces.c55_pow + mi*slice , A.x, A.y);
                    c56_pow = BilinearInterpolation(xvec, yvec, interfaces.c56_pow + mi*slice , A.x, A.y);
                    c66_pow = BilinearInterpolation(xvec, yvec, interfaces.c66_pow + mi*slice , A.x, A.y);
                    rho_pow = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice , A.x, A.y);
                    c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice , A.x, A.y);
                    c12_grad = BilinearInterpolation(xvec, yvec, interfaces.c12_grad + mi*slice , A.x, A.y);
                    c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice , A.x, A.y);
                    c14_grad = BilinearInterpolation(xvec, yvec, interfaces.c14_grad + mi*slice , A.x, A.y);
                    c15_grad = BilinearInterpolation(xvec, yvec, interfaces.c15_grad + mi*slice , A.x, A.y);
                    c16_grad = BilinearInterpolation(xvec, yvec, interfaces.c16_grad + mi*slice , A.x, A.y);
                    c22_grad = BilinearInterpolation(xvec, yvec, interfaces.c22_grad + mi*slice , A.x, A.y);
                    c23_grad = BilinearInterpolation(xvec, yvec, interfaces.c23_grad + mi*slice , A.x, A.y);
                    c24_grad = BilinearInterpolation(xvec, yvec, interfaces.c24_grad + mi*slice , A.x, A.y);
                    c25_grad = BilinearInterpolation(xvec, yvec, interfaces.c25_grad + mi*slice , A.x, A.y);
                    c26_grad = BilinearInterpolation(xvec, yvec, interfaces.c26_grad + mi*slice , A.x, A.y);
                    c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice , A.x, A.y);
                    c34_grad = BilinearInterpolation(xvec, yvec, interfaces.c34_grad + mi*slice , A.x, A.y);
                    c35_grad = BilinearInterpolation(xvec, yvec, interfaces.c35_grad + mi*slice , A.x, A.y);
                    c36_grad = BilinearInterpolation(xvec, yvec, interfaces.c36_grad + mi*slice , A.x, A.y);
                    c44_grad = BilinearInterpolation(xvec, yvec, interfaces.c44_grad + mi*slice , A.x, A.y);
                    c45_grad = BilinearInterpolation(xvec, yvec, interfaces.c45_grad + mi*slice , A.x, A.y);
                    c46_grad = BilinearInterpolation(xvec, yvec, interfaces.c46_grad + mi*slice , A.x, A.y);
                    c55_grad = BilinearInterpolation(xvec, yvec, interfaces.c55_grad + mi*slice , A.x, A.y);
                    c56_grad = BilinearInterpolation(xvec, yvec, interfaces.c56_grad + mi*slice , A.x, A.y);
                    c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice , A.x, A.y);
                    rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice , A.x, A.y);
                    float var0_top = rho + pow(dz, rho_pow) * rho_grad;
                    float var1_top = c11 + pow(dz, c11_pow) * c11_grad; 
                    float var2_top = c12 + pow(dz, c12_pow) * c12_grad;
                    float var3_top = c13 + pow(dz, c13_pow) * c13_grad;
                    float var4_top = c14 + pow(dz, c14_pow) * c14_grad;
                    float var5_top = c15 + pow(dz, c15_pow) * c15_grad;
                    float var6_top = c16 + pow(dz, c16_pow) * c16_grad;
                    float var7_top = c22 + pow(dz, c22_pow) * c22_grad;
                    float var8_top = c23 + pow(dz, c23_pow) * c23_grad;
                    float var9_top = c24 + pow(dz, c24_pow) * c24_grad;
                    float var10_top = c25 + pow(dz, c25_pow) * c25_grad;
                    float var11_top = c26 + pow(dz, c26_pow) * c26_grad;
                    float var12_top = c33 + pow(dz, c33_pow) * c33_grad;
                    float var13_top = c34 + pow(dz, c34_pow) * c34_grad;
                    float var14_top = c35 + pow(dz, c35_pow) * c35_grad;
                    float var15_top = c36 + pow(dz, c36_pow) * c36_grad;
                    float var16_top = c44 + pow(dz, c44_pow) * c44_grad;
                    float var17_top = c45 + pow(dz, c45_pow) * c45_grad;
                    float var18_top = c46 + pow(dz, c46_pow) * c46_grad;
                    float var19_top = c55 + pow(dz, c55_pow) * c55_grad;
                    float var20_top = c56 + pow(dz, c56_pow) * c56_grad;
                    float var21_top = c66 + pow(dz, c66_pow) * c66_grad;
                    var[0]  = (var[0]  + var0_top ) / 2.0;
                    var[1]  = (var[1]  + var1_top ) / 2.0;
                    var[2]  = (var[2]  + var2_top ) / 2.0;
                    var[3]  = (var[3]  + var3_top ) / 2.0;
                    var[4]  = (var[4]  + var4_top ) / 2.0;
                    var[5]  = (var[5]  + var5_top ) / 2.0;
                    var[6]  = (var[6]  + var6_top ) / 2.0;
                    var[7]  = (var[7]  + var7_top ) / 2.0;
                    var[8]  = (var[8]  + var8_top ) / 2.0;
                    var[9]  = (var[9]  + var9_top ) / 2.0;
                    var[10] = (var[10] + var10_top) / 2.0;
                    var[11] = (var[11] + var11_top) / 2.0;
                    var[12] = (var[12] + var12_top) / 2.0;
                    var[13] = (var[13] + var13_top) / 2.0;
                    var[14] = (var[14] + var14_top) / 2.0;
                    var[15] = (var[15] + var15_top) / 2.0;
                    var[16] = (var[16] + var16_top) / 2.0;
                    var[17] = (var[17] + var17_top) / 2.0;
                    var[18] = (var[18] + var18_top) / 2.0;
                    var[19] = (var[19] + var19_top) / 2.0;
                    var[20] = (var[20] + var20_top) / 2.0;
                    var[21] = (var[21] + var21_top) / 2.0;
                }
            }
        }     
    break;

    case ELASTIC_TTI_BOND: /* 6. rho c11 c33 c55 c66 c13 azimuth dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
        } else {
            float c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c55 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice, A.x, A.y);
            float c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice, A.x, A.y);
            float c55_grad = BilinearInterpolation(xvec, yvec, interfaces.c55_grad + mi*slice, A.x, A.y);
            float c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice, A.x, A.y);
            float c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
            float c11_pow  = BilinearInterpolation(xvec, yvec, interfaces.c11_pow  + mi*slice, A.x, A.y);
            float c33_pow  = BilinearInterpolation(xvec, yvec, interfaces.c33_pow  + mi*slice, A.x, A.y);
            float c55_pow  = BilinearInterpolation(xvec, yvec, interfaces.c55_pow  + mi*slice, A.x, A.y);
            float c66_pow  = BilinearInterpolation(xvec, yvec, interfaces.c66_pow  + mi*slice, A.x, A.y);
            float c13_pow  = BilinearInterpolation(xvec, yvec, interfaces.c13_pow  + mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
            float azimuth      = BilinearInterpolation(xvec, yvec, interfaces.azimuth      + mi*slice, A.x, A.y);
            float dip          = BilinearInterpolation(xvec, yvec, interfaces.dip          + mi*slice, A.x, A.y);
            float azimuth_grad = BilinearInterpolation(xvec, yvec, interfaces.azimuth_grad + mi*slice, A.x, A.y);
            float dip_grad     = BilinearInterpolation(xvec, yvec, interfaces.dip_grad     + mi*slice, A.x, A.y);
            float azimuth_pow  = BilinearInterpolation(xvec, yvec, interfaces.azimuth_pow  + mi*slice, A.x, A.y);
            float dip_pow      = BilinearInterpolation(xvec, yvec, interfaces.dip_pow      + mi*slice, A.x, A.y);
            var[0] = rho + pow(dz, rho_pow)* rho_grad;
            var[1] = c11 + pow(dz, c11_pow)* c11_grad;
            var[2] = c33 + pow(dz, c33_pow)* c33_grad;
            var[3] = c55 + pow(dz, c55_pow)* c55_grad;
            var[4] = c66 + pow(dz, c66_pow)* c66_grad;
            var[5] = c13 + pow(dz, c13_pow)* c13_grad;
            var[6] = azimuth + pow(dz, azimuth_pow)* azimuth_grad;
            var[7] = dip     + pow(dz, dip_pow)    * dip_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
                    c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
                    c55 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
                    c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
                    c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
                    rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
                    c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice, A.x, A.y);
                    c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice, A.x, A.y);
                    c55_grad = BilinearInterpolation(xvec, yvec, interfaces.c55_grad + mi*slice, A.x, A.y);
                    c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice, A.x, A.y);
                    c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice, A.x, A.y);
                    rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
                    c11_pow  = BilinearInterpolation(xvec, yvec, interfaces.c11_pow  + mi*slice, A.x, A.y);
                    c33_pow  = BilinearInterpolation(xvec, yvec, interfaces.c33_pow  + mi*slice, A.x, A.y);
                    c55_pow  = BilinearInterpolation(xvec, yvec, interfaces.c55_pow  + mi*slice, A.x, A.y);
                    c66_pow  = BilinearInterpolation(xvec, yvec, interfaces.c66_pow  + mi*slice, A.x, A.y);
                    c13_pow  = BilinearInterpolation(xvec, yvec, interfaces.c13_pow  + mi*slice, A.x, A.y);
                    rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
                    azimuth      = BilinearInterpolation(xvec, yvec, interfaces.azimuth      + mi*slice, A.x, A.y);
                    dip          = BilinearInterpolation(xvec, yvec, interfaces.dip          + mi*slice, A.x, A.y);
                    azimuth_grad = BilinearInterpolation(xvec, yvec, interfaces.azimuth_grad + mi*slice, A.x, A.y);
                    dip_grad     = BilinearInterpolation(xvec, yvec, interfaces.dip_grad     + mi*slice, A.x, A.y);
                    azimuth_pow  = BilinearInterpolation(xvec, yvec, interfaces.azimuth_pow  + mi*slice, A.x, A.y);
                    dip_pow      = BilinearInterpolation(xvec, yvec, interfaces.dip_pow      + mi*slice, A.x, A.y); 
                    float var0_top = rho + pow(dz, rho_pow)* rho_grad;
                    float var1_top = c11 + pow(dz, c11_pow)* c11_grad;
                    float var2_top = c33 + pow(dz, c33_pow)* c33_grad;
                    float var3_top = c55 + pow(dz, c55_pow)* c55_grad;
                    float var4_top = c66 + pow(dz, c66_pow)* c66_grad;
                    float var5_top = c13 + pow(dz, c13_pow)* c13_grad;
                    float var6_top = azimuth + pow(dz, azimuth_pow)* azimuth_grad;
                    float var7_top = dip     + pow(dz, dip_pow)    * dip_grad; 
                    var[0] = (var[0] + var0_top)/2.0;
                    var[1] = (var[1] + var1_top)/2.0;
                    var[2] = (var[2] + var2_top)/2.0;
                    var[3] = (var[3] + var3_top)/2.0;
                    var[4] = (var[4] + var4_top)/2.0;
                    var[5] = (var[5] + var5_top)/2.0;
                    var[6] = (var[6] + var6_top)/2.0;
                    var[7] = (var[7] + var7_top)/2.0;             
                }
            }
        }   
    break;

    case ACOUSTIC_ISOTROPIC: /* 7. rho vp */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
        } else {
            float vp_grad  = BilinearInterpolation(xvec, yvec, interfaces.vp_grad + mi*slice, A.x, A.y);
            float vp       = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice, A.x, A.y);
            float vp_pow   = BilinearInterpolation(xvec, yvec, interfaces.vp_pow  + mi*slice, A.x, A.y);
            float rho      = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad+ mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice, A.x, A.y);
            var[0] = rho + pow(dz,rho_pow)*rho_grad;
            var[1] = vp  + pow(dz, vp_pow)* vp_grad;
            if (isEqual(dz, 0.0)) {
                mi = findFirstGreaterEqualIndex(A.z, elevation)-1;
                if (mi >= 0) {
                    vp_grad  = BilinearInterpolation(xvec, yvec, interfaces.vp_grad + mi*slice, A.x, A.y);
                    vp       = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice, A.x, A.y);
                    vp_pow   = BilinearInterpolation(xvec, yvec, interfaces.vp_pow  + mi*slice, A.x, A.y);
                    rho      = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
                    rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad+ mi*slice, A.x, A.y);
                    rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice, A.x, A.y);
                    float var0_top = rho + pow(dz,rho_pow)*rho_grad;
                    float var1_top = vp  + pow(dz, vp_pow)* vp_grad;
                    var[0] = (var[0] + var0_top)/2.0;
                    var[1] = (var[1] + var1_top)/2.0;
                }
            }
        }
    break;

    default: // for self-check
        fprintf(stderr,"Error: Unknow media (for code check, please contact Luqian Jiang)\n");
        fflush(stderr);
        exit(1);

    }     
} 

//- 0. assign the parameter directly (use the local values): one component 
void parametrization_layer_onecmp_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *var3d,
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                
                std::vector<float> var(1, 0.0);

                size_t indx =  i + j * siz_line + k * siz_slice;
                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y

                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }
                
                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ONE_COMPONENT, var);
    
                var3d[indx] = var[0];     
            }
        }
    }

}

//- 1. assign the parameter directly (use the local values): isotropic, acoustic 
void parametrization_layer_ac_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *kappa,
    float *rho3d,
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                std::vector<float> var(2, 0.0); // rho, vp

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y

                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ACOUSTIC_ISOTROPIC, var);

                /*calculate output kappa */
                float vp = var[1], rho = var[0];
                kappa[indx] = vp*vp*rho;
                rho3d[indx] = rho; 
            }
        }
    }
}


//- 2. assign the parameter directly (use the local values): elastic isotropic 
void parametrization_layer_el_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d, 
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                std::vector<float> var(3, 0.0); // rho, vp, vs

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_ISOTROPIC, var);

                /*calculate output lambda and mu */
                float vp = var[1], vs = var[2], rho = var[0];
                mu3d[indx]  = vs*vs*rho;
                lam3d[indx] = vp*vp*rho - 2.0*mu3d[indx];
                rho3d[indx] = rho; 
                
            }
        }
    }
}


//- 3. Assign the parameter directly (use the local values): elastic vti
void parametrization_layer_el_vti_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    int media_type = interfaces.media_type;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {


                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                std::vector<float> var(6, 0.0); 
                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, media_type, var);
                
                para2vti(var, media_type, 
                    c11[indx], c33[indx], c55[indx], c66[indx], c13[indx], rho[indx]);
                
            }
        }
    }
}

//- 4. Assign the parameter directly (use the local values): elastic tti
void parametrization_layer_el_aniso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c12,
    float *c13,
    float *c14,
    float *c15,
    float *c16,
    float *c22,
    float *c23,
    float *c24,
    float *c25,
    float *c26,
    float *c33,
    float *c34,
    float *c35,
    float *c36,
    float *c44,
    float *c45,
    float *c46,
    float *c55,
    float *c56,
    float *c66,
    float *rho,
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;
    int media_type = interfaces.media_type;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y

                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                std::vector<float> var(22, 0.0); 

                AssignLayerMediaPara2Point(i, j, k,
                            Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                            interfaces, media_type, var);

                para2tti(var, media_type, // return cij
                    c11[indx], c12[indx], c13[indx], c14[indx], c15[indx], c16[indx],
                    c22[indx], c23[indx], c24[indx], c25[indx], c26[indx],
                    c33[indx], c34[indx], c35[indx], c36[indx],
                    c44[indx], c45[indx], c46[indx], 
                    c55[indx], c56[indx], c66[indx], rho[indx]);
            }
        }
    }
}

//======================== averaging/equivalent medium method ==================================

/* 
 * For half-grid point, marked the materials number.  
 */
void MarkInterfaceNumber(
        int grid_type,
        float *Hx, float *Hy, float *Hz,
        size_t nx, size_t ny, size_t nz,
        int *MaterNum, // nx*ny*nz
        inter_t &interfaces) 
{
    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    size_t  NX = interfaces.NX;
    size_t  NY = interfaces.NY;
    size_t  NI = interfaces.NI;
    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t inter_slice = NX*NY;

    std::vector<float> XVEC(NX), YVEC(NY);
   
    // x-vector and y-vector of interfaces mesh
    for (size_t i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (size_t i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }

    
    if (grid_type != GRID_CURV) {

        // for every half grid(Hx[i], Hy[j]), store the elevation of every ni.
        // elevation[indx_in_xy][ni]
        std::vector<std::vector<float>> elevation(siz_slice, std::vector<float>(NI, -FLT_MAX));
        
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                size_t indx_in_slice = i + j*nx;
                for (size_t ni = 0; ni < NI; ni++) {
                    elevation[indx_in_slice][ni] = BilinearInterpolation(
                        XVEC, YVEC, interfaces.elevation+ni*inter_slice, Hx[i], Hy[j]);
                }

                for (size_t k = 0; k < nz; k++) {
                    size_t indx = indx_in_slice + k*siz_slice;
                    /* use which material */
                    int mi = -1;
                    if (grid_type == GRID_VMAP) {
                        mi = findLastGreaterEqualIndex(Hz[indx], elevation[indx_in_slice]);
                    } else {
                        mi = findLastGreaterEqualIndex(Hz[k], elevation[indx_in_slice]);
                    }
                    if (mi == -1) {
                        mi = findLastGreaterEqualIndex(elevation[indx_in_slice][0], 
                                                       elevation[indx_in_slice]);
                    }

                    MaterNum[indx] = mi;
                }
            }
        }

    } else {
        std::vector<float> elevation(NI, -FLT_MAX);
        for (size_t indx = 0; indx < siz_volume; indx++) {
            for (size_t ni = 0; ni < NI; ni++) {
                /* Get the elevation for the corresponding location */
                elevation[ni] = BilinearInterpolation(XVEC, YVEC, 
                    interfaces.elevation + ni*inter_slice, Hx[indx], Hy[indx]);
            }
            /* use which material */
            int mi = findLastGreaterEqualIndex(Hz[indx], elevation);
            if (mi == -1) {
                mi = findLastGreaterEqualIndex(elevation[0], elevation);
            }
            MaterNum[indx] = mi;
        }
    } 
    //else {
   //     fprintf(stderr,"Error: Unknow grid_type, please check the code! (for code check, please contact Luqian Jiang)");
   //     fflush(stderr);
   //     exit(1);
   // }
}

//- 0.1 Assign the parameter by volume harmonic averaging
//- one component
void parametrization_layer_onecmp_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *var3d,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_onecmp_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, var3d, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float vol_var = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    size_t nsg = (NG+1)*(NG+1)*(NG+1);

                    for (size_t isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(1, 0.0);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, ONE_COMPONENT, var);

                        vol_var += (1.0/var[0]);
                    }

                    var3d[indx] = nsg*1.0/vol_var;

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;
}

//- 0.1 Assign the parameter by volume arithmetic averaging
//- one component
void parametrization_layer_onecmp_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *var3d,
    int myid) 
{
    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_onecmp_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, var3d, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    // for printing progress
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float vol_var = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    size_t nsg = (NG+1)*(NG+1)*(NG+1);

                    // numerical integration
                    for (size_t isg = 0; isg < nsg; isg++) {
                        std::vector<float> var(1, 0.0);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, ONE_COMPONENT, var);

                        vol_var += var[0];
                    }

                    var3d[indx] = vol_var/(nsg*1.0);

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;
}

//- 1.0 Assign the parameter by volume harmonic averaging (kappa)
//- acoustic isortopic
void parametrization_layer_ac_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *kappa, 
    float *rho3d,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_ac_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, kappa, rho3d, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic averaging):\n\n";

    for (size_t j = 1; j < ny-1; ++j) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; ++k) {
            for (size_t i = 1; i < nx-1; ++i) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    Point3 *SubGrid = MeshSubdivide(M);
                    size_t nsg = (NG+1)*(NG+1)*(NG+1);

                    // recalculate the material value of the point
                    float ari_rho = 0.0, har_kappa = 0.0;

                    for (size_t isg = 0; isg < nsg; isg++) {
                        std::vector<float> var(2, 0.0);

                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, ACOUSTIC_ISOTROPIC, var);
                        
                        float vp = var[1], rho = var[0];

                        har_kappa += (1.0/(vp*vp*rho));
                        ari_rho   += rho;
                    }

                    rho3d[indx] = ari_rho/(1.0*nsg);
                    kappa[indx] = 1.0*nsg/har_kappa;

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;
}

//- 1.1 Assign the parameter by volume arithmetic averaging (kappa)
//- acoustic isortopic
void parametrization_layer_ac_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *kappa, 
    float *rho3d,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_ac_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, kappa, rho3d, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    Point3 *SubGrid = MeshSubdivide(M);
                    size_t nsg = (NG+1)*(NG+1)*(NG+1);

                    // recalculate the material value of the point
                    float ari_rho = 0.0, ari_kappa = 0.0;

                    for (size_t isg = 0; isg < nsg; isg++) {
                        std::vector<float> var(2, 0.0);

                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, ACOUSTIC_ISOTROPIC, var);
                        
                        float vp = var[1], rho = var[0];

                        ari_kappa += (vp*vp*rho);
                        ari_rho   += rho;
                    }

                    rho3d[indx] = ari_rho/(1.0*nsg);
                    kappa[indx] = ari_kappa/(1.0*nsg);

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;
}



//- 2.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic isotropic
// ref: Moczo, 2002, 3D Heterogeneous Staggered-Grid Finite-Difference Modeling of Seismic
//      Motion with Volume Harmonic and Arithmetic Averaging of Elastic Moduli and Densities
void parametrization_layer_el_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    float slowk = 1.0/(ny-2); // for print progress


    // assign the local value first.
    parametrization_layer_el_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, lam3d, mu3d, rho3d, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- volume arithmetic and harmonic averaging method:\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho   = 0.0;
                    float har_kappa = 0.0;
                    float har_mu    = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var);

                        float vp = var[1], vs = var[2], rho = var[0];
                        float mu = vs*vs*rho;
                        float lambda = vp*vp*rho - 2.0*mu;

                        ari_rho   += rho;
                        har_kappa += (1.0/(lambda + 2.0/3.0*mu));
                        har_mu    += (1.0/mu);
                    }

                    har_mu    = nsg*1.0/har_mu;
                    har_kappa = nsg*1.0/har_kappa;
                    ari_rho   = ari_rho/(nsg*1.0); 

                    lam3d[indx] = har_kappa - 2.0/3.0*har_mu;
                    mu3d[indx]  = har_mu;
                    rho3d[indx] = ari_rho;

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;
}

//- 2.1 Assign the parameter by volume arithmetic averaging method 
//- elasic isotropic
void parametrization_layer_el_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_el_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, lam3d, mu3d, rho3d, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid==0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float ari_lam = 0.0;
                    float ari_mu  = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3, 0.0);

                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var);

                        float vp = var[1], vs = var[2], rho = var[0];
                        float mu = vs*vs*rho;
                        float lambda = vp*vp*rho - 2.0*mu;

                        ari_rho += rho;
                        ari_lam += lambda;
                        ari_mu  += mu;
                    }

                    lam3d[indx] = ari_lam/(nsg*1.0);
                    mu3d[indx]  = ari_mu /(nsg*1.0);
                    rho3d[indx] = ari_rho/(nsg*1.0);

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;
}

//- 3.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic vti
void parametrization_layer_el_vti_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;    
    int media_type = interfaces.media_type;
    float slowk = 1.0/(ny-2); // for print progress


    // assign the local value first.
    parametrization_layer_el_vti_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, c11, c33, c55, c66, c13, rho, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float har_c11 = 0.0;
                    float har_c33 = 0.0;
                    float har_c55 = 0.0;
                    float har_c66 = 0.0;
                    float har_c13 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(6, 0.0);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, media_type, var);

                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c33_p = 0.0;
                        float c55_p = 0.0, c66_p = 0.0;
                        float c13_p = 0.0, rho_p = 0.0;
                        para2vti(var, media_type, 
                            c11_p, c33_p, c55_p, c66_p, c13_p, rho_p);

                        har_c11 += (1.0/c11_p);
                        har_c13 += (1.0/c13_p);
                        har_c33 += (1.0/c33_p);
                        har_c55 += (1.0/c55_p);
                        har_c66 += (1.0/c66_p);
                        ari_rho += rho_p;
                    }

                    c11[indx] = 1.0*nsg/har_c11;
                    c13[indx] = 1.0*nsg/har_c13;
                    c33[indx] = 1.0*nsg/har_c33;
                    c55[indx] = 1.0*nsg/har_c55;
                    c66[indx] = 1.0*nsg/har_c66;
                    rho[indx] = ari_rho/nsg;

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;    
}

//- 3.1 Assign the parameter by volume arithmetic averaging method 
//- elasic vti
void parametrization_layer_el_vti_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    int media_type = interfaces.media_type;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_el_vti_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, c11, c33, c55, c66, c13, rho, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float ari_c11 = 0.0;
                    float ari_c33 = 0.0;
                    float ari_c55 = 0.0;
                    float ari_c66 = 0.0;
                    float ari_c13 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(6, 0.0);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, media_type, var);

                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c33_p = 0.0;
                        float c55_p = 0.0, c66_p = 0.0;
                        float c13_p = 0.0, rho_p = 0.0;
                        para2vti(var, media_type, 
                            c11_p, c33_p, c55_p, c66_p, c13_p, rho_p);

                        ari_c11 += c11_p;
                        ari_c13 += c13_p;
                        ari_c33 += c33_p;
                        ari_c55 += c55_p;
                        ari_c66 += c66_p;
                        ari_rho += rho_p;
                    }

                    c11[indx] = ari_c11/nsg;
                    c13[indx] = ari_c13/nsg;
                    c33[indx] = ari_c33/nsg;
                    c55[indx] = ari_c55/nsg;
                    c66[indx] = ari_c66/nsg;
                    rho[indx] = ari_rho/nsg;

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;    
}


//- 4.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic tti
void parametrization_layer_el_aniso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c12,
    float *c13,
    float *c14,
    float *c15,
    float *c16,
    float *c22,
    float *c23,
    float *c24,
    float *c25,
    float *c26,
    float *c33,
    float *c34,
    float *c35,
    float *c36,
    float *c44,
    float *c45,
    float *c46,
    float *c55,
    float *c56,
    float *c66,
    float *rho,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    int media_type = interfaces.media_type;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_el_aniso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, c11, c12, c13, c14, c15, c16,
        c22, c23, c24, c25, c26, c33, c34, c35, c36, 
        c44, c45, c46, c55, c56, c66, rho, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float har_c11 = 0.0;
                    float har_c12 = 0.0;
                    float har_c13 = 0.0;
                    float har_c14 = 0.0;
                    float har_c15 = 0.0;
                    float har_c16 = 0.0;
                    float har_c22 = 0.0;
                    float har_c23 = 0.0;
                    float har_c24 = 0.0;
                    float har_c25 = 0.0;
                    float har_c26 = 0.0;
                    float har_c33 = 0.0;
                    float har_c34 = 0.0;
                    float har_c35 = 0.0;
                    float har_c36 = 0.0;
                    float har_c44 = 0.0;
                    float har_c45 = 0.0;
                    float har_c46 = 0.0;
                    float har_c55 = 0.0;
                    float har_c56 = 0.0;
                    float har_c66 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(22, 0.0);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, media_type, var);

                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c12_p = 0.0, c13_p = 0.0;
                        float c14_p = 0.0, c15_p = 0.0, c16_p = 0.0;
                        float c22_p = 0.0, c23_p = 0.0, c24_p = 0.0;
                        float c25_p = 0.0, c26_p = 0.0, c33_p = 0.0;
                        float c34_p = 0.0, c35_p = 0.0, c36_p = 0.0;
                        float c44_p = 0.0, c45_p = 0.0, c46_p = 0.0;
                        float c55_p = 0.0, c56_p = 0.0, c66_p = 0.0;
                        float rho_p = 0.0;

                        para2tti(var, media_type, // return cij
                            c11_p, c12_p, c13_p, c14_p, c15_p, c16_p,
                            c22_p, c23_p, c24_p, c25_p, c26_p,
                            c33_p, c34_p, c35_p, c36_p,
                            c44_p, c45_p, c46_p, 
                            c55_p, c56_p, c66_p, rho_p);

                        har_c11 += (1.0/c11_p);
                        har_c12 += (1.0/c12_p);
                        har_c13 += (1.0/c13_p);
                        har_c14 += (1.0/c14_p);
                        har_c15 += (1.0/c15_p);
                        har_c16 += (1.0/c16_p);
                        har_c22 += (1.0/c22_p);
                        har_c23 += (1.0/c23_p);
                        har_c24 += (1.0/c24_p);
                        har_c25 += (1.0/c25_p);
                        har_c26 += (1.0/c26_p);
                        har_c33 += (1.0/c33_p);
                        har_c34 += (1.0/c34_p);
                        har_c35 += (1.0/c35_p);
                        har_c36 += (1.0/c36_p);
                        har_c44 += (1.0/c44_p);
                        har_c45 += (1.0/c45_p);
                        har_c46 += (1.0/c46_p);
                        har_c55 += (1.0/c55_p);
                        har_c56 += (1.0/c56_p);
                        har_c66 += (1.0/c66_p);
                        
                        ari_rho += rho_p;
                    }

                    c11[indx] = (1.0*nsg)/har_c11;
                    c12[indx] = (1.0*nsg)/har_c12;
                    c13[indx] = (1.0*nsg)/har_c13;
                    c14[indx] = (1.0*nsg)/har_c14;
                    c15[indx] = (1.0*nsg)/har_c15;
                    c16[indx] = (1.0*nsg)/har_c16;
                    c22[indx] = (1.0*nsg)/har_c22;
                    c23[indx] = (1.0*nsg)/har_c23;
                    c24[indx] = (1.0*nsg)/har_c24;
                    c25[indx] = (1.0*nsg)/har_c25;
                    c26[indx] = (1.0*nsg)/har_c26;
                    c33[indx] = (1.0*nsg)/har_c33;
                    c34[indx] = (1.0*nsg)/har_c34;
                    c35[indx] = (1.0*nsg)/har_c35;
                    c36[indx] = (1.0*nsg)/har_c36;
                    c44[indx] = (1.0*nsg)/har_c44;
                    c45[indx] = (1.0*nsg)/har_c45;
                    c46[indx] = (1.0*nsg)/har_c46;
                    c55[indx] = (1.0*nsg)/har_c55;
                    c56[indx] = (1.0*nsg)/har_c56;
                    c66[indx] = (1.0*nsg)/har_c66;

                    rho[indx] = ari_rho/nsg;

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;    
}

//- 4.1 Assign the parameter by volume arithmetic averaging method 
//- elasic tti
void parametrization_layer_el_aniso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c12,
    float *c13,
    float *c14,
    float *c15,
    float *c16,
    float *c22,
    float *c23,
    float *c24,
    float *c25,
    float *c26,
    float *c33,
    float *c34,
    float *c35,
    float *c36,
    float *c44,
    float *c45,
    float *c46,
    float *c55,
    float *c56,
    float *c66,
    float *rho,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;
    int media_type = interfaces.media_type;
    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_el_aniso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, c11, c12, c13, c14, c15, c16,
        c22, c23, c24, c25, c26, c33, c34, c35, c36, 
        c44, c45, c46, c55, c56, c66, rho, myid);

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slowk*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);

                // conter-clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    
                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float ari_c11 = 0.0;
                    float ari_c12 = 0.0;
                    float ari_c13 = 0.0;
                    float ari_c14 = 0.0;
                    float ari_c15 = 0.0;
                    float ari_c16 = 0.0;
                    float ari_c22 = 0.0;
                    float ari_c23 = 0.0;
                    float ari_c24 = 0.0;
                    float ari_c25 = 0.0;
                    float ari_c26 = 0.0;
                    float ari_c33 = 0.0;
                    float ari_c34 = 0.0;
                    float ari_c35 = 0.0;
                    float ari_c36 = 0.0;
                    float ari_c44 = 0.0;
                    float ari_c45 = 0.0;
                    float ari_c46 = 0.0;
                    float ari_c55 = 0.0;
                    float ari_c56 = 0.0;
                    float ari_c66 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(22, 0.0);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, media_type, var);

                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c12_p = 0.0, c13_p = 0.0;
                        float c14_p = 0.0, c15_p = 0.0, c16_p = 0.0;
                        float c22_p = 0.0, c23_p = 0.0, c24_p = 0.0;
                        float c25_p = 0.0, c26_p = 0.0, c33_p = 0.0;
                        float c34_p = 0.0, c35_p = 0.0, c36_p = 0.0;
                        float c44_p = 0.0, c45_p = 0.0, c46_p = 0.0;
                        float c55_p = 0.0, c56_p = 0.0, c66_p = 0.0;
                        float rho_p = 0.0;

                        para2tti(var, media_type, // return cij
                            c11_p, c12_p, c13_p, c14_p, c15_p, c16_p,
                            c22_p, c23_p, c24_p, c25_p, c26_p,
                            c33_p, c34_p, c35_p, c36_p,
                            c44_p, c45_p, c46_p, 
                            c55_p, c56_p, c66_p, rho_p);

                        ari_c11 += c11_p;
                        ari_c12 += c12_p;
                        ari_c13 += c13_p;
                        ari_c14 += c14_p;
                        ari_c15 += c15_p;
                        ari_c16 += c16_p;
                        ari_c22 += c22_p;
                        ari_c23 += c23_p;
                        ari_c24 += c24_p;
                        ari_c25 += c25_p;
                        ari_c26 += c26_p;
                        ari_c33 += c33_p;
                        ari_c34 += c34_p;
                        ari_c35 += c35_p;
                        ari_c36 += c36_p;
                        ari_c44 += c44_p;
                        ari_c45 += c45_p;
                        ari_c46 += c46_p;
                        ari_c55 += c55_p;
                        ari_c56 += c56_p;
                        ari_c66 += c66_p;
                        ari_rho += rho_p;
                    }

                    c11[indx] = ari_c11/nsg;
                    c12[indx] = ari_c12/nsg;
                    c13[indx] = ari_c13/nsg;
                    c14[indx] = ari_c14/nsg;
                    c15[indx] = ari_c15/nsg;
                    c16[indx] = ari_c16/nsg;
                    c22[indx] = ari_c22/nsg;
                    c23[indx] = ari_c23/nsg;
                    c24[indx] = ari_c24/nsg;
                    c25[indx] = ari_c25/nsg;
                    c26[indx] = ari_c26/nsg;
                    c33[indx] = ari_c33/nsg;
                    c34[indx] = ari_c34/nsg;
                    c35[indx] = ari_c35/nsg;
                    c36[indx] = ari_c36/nsg;
                    c44[indx] = ari_c44/nsg;
                    c45[indx] = ari_c45/nsg;
                    c46[indx] = ari_c46/nsg;
                    c55[indx] = ari_c55/nsg;
                    c56[indx] = ari_c56/nsg;
                    c66[indx] = ari_c66/nsg;
                    
                    rho[indx] = ari_rho/nsg;

                    if(SubGrid != nullptr) delete []SubGrid;
                }

            }
        }
    }

    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;    
}




void parametrization_layer_el_iso_tti(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *rho,
    float *c11,
    float *c12,
    float *c13,
    float *c14,
    float *c15,
    float *c16,
    float *c22,
    float *c23,
    float *c24,
    float *c25,
    float *c26,
    float *c33,
    float *c34,
    float *c35,
    float *c36,
    float *c44,
    float *c45,
    float *c46,
    float *c55,
    float *c56,
    float *c66,
    int myid)  
{
    size_t siz_line   = nx;
    size_t siz_slice  = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    size_t NI = interfaces.NI;

    float slowk = 1.0/(ny-2); // for print progress

    // assign the local value first.
    parametrization_layer_el_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz, grid_type, 
                interfaces, c13, c44, rho, myid);
    for (size_t i = 0; i < siz_volume; ++i) {
        c11[i] = c13[i] + 2.0*c44[i]; 
        c22[i] = c11[i]; c33[i] = c11[i]; 
        c12[i] = c13[i]; c23[i] = c13[i];
        c55[i] = c44[i]; c66[i] = c44[i];
        c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
        c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
        c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
        c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
    }

    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    /* mark the interface number at the half grid. */
    int *MaterNum = new int[siz_volume];

    MarkInterfaceNumber(grid_type, Hx, Hy, Hz, 
        nx, ny, nz, MaterNum, interfaces);

    int inter_nx = interfaces.NX, inter_ny = interfaces.NY;
    float inter_dx = interfaces.DX, inter_dy = interfaces.DY;
    float inter_x0 = interfaces.MINX, inter_y0 = interfaces.MINY;
    size_t inter_slice = inter_nx*inter_ny;

    if (myid == 0)
        std::cout << "- TTI equivalent medium parametrization:\n\n";
    for (size_t j = 1; j < ny-1; j++) {
      if (myid == 0) printProgress(slowk*j);
      
      for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
          size_t indx =  i + j * siz_line + k * siz_slice; 
          // check if the mesh have different values;
          std::vector<int> v(8);

          // clockwise, Consistent with the schematic diagram
          v[0] = MaterNum[indx-1-siz_line-siz_slice];
          v[1] = MaterNum[indx-1         -siz_slice];
          v[2] = MaterNum[indx           -siz_slice];
          v[3] = MaterNum[indx  -siz_line-siz_slice];
          v[4] = MaterNum[indx-1-siz_line];
          v[5] = MaterNum[indx-1         ];
          v[6] = MaterNum[indx           ];
          v[7] = MaterNum[indx  -siz_line];


          // There is more than one medium value in the half-grid mesh, 
          int medium_num = NumOfValues(v, NI);
          if ( medium_num > 1) {

            Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                    siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

            // There are 2 media in the cell.  
            // One interface, tti; two or more interface, ort
            bool is_ort = false;
            if (medium_num > 2) {
              is_ort = true;
            } 
            else  // medium_num==2
            {  
              int avg_number = std::accumulate(v.begin(), v.end(), 0);
              avg_number /= 8;

              // passed by which interface
              int interface_label = -1;
              for (auto &v_label:v) 
                if (v_label > interface_label)
                  interface_label = v_label;

              // compute cubeState (the index into the edgeTable)
              int cubeState = 0;
              if (v[0] <= avg_number) cubeState |= 1;
              if (v[1] <= avg_number) cubeState |= 2;
              if (v[2] <= avg_number) cubeState |= 4;
              if (v[3] <= avg_number) cubeState |= 8;
              if (v[4] <= avg_number) cubeState |= 16;
              if (v[5] <= avg_number) cubeState |= 32;
              if (v[6] <= avg_number) cubeState |= 64;
              if (v[7] <= avg_number) cubeState |= 128;

              // get the intersection Points, use set to avoid duplication
              std::set<Point3> intersectionList;
              // the interface pass through the cell two or more times, 
              //   can not compute the normal of the interface plane, apply ort
              int edgeState = edgeTable[cubeState];
              if (edgeState == 0) 
                is_ort = true;

              int listsize = 0;
              getInersectList(inter_nx, inter_ny, inter_dx, inter_dy, inter_x0, inter_y0, 
                interfaces.elevation+(avg_number+1)*inter_slice, M, edgeState, intersectionList );
              listsize = intersectionList.size();
              
              if(listsize < 3 || listsize > 6) {
                is_ort = true;
              } else {  // listsize in [3,6]
                Plane inter_plane(0.0, 0.0, 0.0, 0.0);
                bool isfitplane= FitHeightPlane(intersectionList, inter_plane);
                if (!isfitplane) {
                  is_ort = true;
                }

                Vector3 a(0.0, 0.0, 0.0), b(0.0, 0.0, 0.0), c(0.0, 0.0, 0.0);
                buildCoord_fromPlane(inter_plane, a, b, c); 

                // cal TTI equivalent medium parameters
                if (!is_ort) {
                  Point3 *SubGrid = MeshSubdivide(M);
                //  Vector3 d = Cross(c, b);
                  cal_iso_tti_parametrization_val(i, j, k, interfaces, SubGrid, a, b, c,
                    c11[indx], c12[indx], c13[indx], c14[indx], c15[indx], c16[indx], 
                    c22[indx], c23[indx], c24[indx], c25[indx], c26[indx], 
                    c33[indx], c34[indx], c35[indx], c36[indx], 
                    c44[indx], c45[indx], c46[indx],
                    c55[indx], c56[indx], c66[indx], rho[indx]);
                  if(SubGrid != nullptr) delete []SubGrid;
                } // tti

              } // listsize

            }  // medium num=2

            // numerical volume averaging
            if (is_ort) {
              Point3 *SubGrid = MeshSubdivide(M);
              cal_ort_represent_val(i,j,k, interfaces, SubGrid,
                c11[indx], c12[indx], c13[indx], 
                c22[indx], c23[indx], c33[indx],
                c44[indx], c55[indx], c66[indx], rho[indx]);
              if(SubGrid != nullptr) delete []SubGrid;
            } 
          } // medium_num > 1

        }
      }
    }
    if (Hx != nullptr) delete [] Hx;
    if (Hy != nullptr) delete [] Hy;
    if (Hz != nullptr) delete [] Hz;
    if (MaterNum != nullptr) delete [] MaterNum;
}

//
// calculate TTI parameterization parameters
// for elasic isotropic
// Jiang and Zhang, 2021, TTI equivalent medium parametrization method for the seismic 
//    waveform modelling of heterogeneous media with coarse grids
//
void cal_iso_tti_parametrization_val(
    int i, int j, int k,
    inter_t &interfaces, Point3 *SubGrid, 
    const Vector3 &a, const Vector3 &b, const Vector3 &c, 
    float &c11, float &c12, float &c13, float &c14, float &c15, float &c16, 
    float &c22, float &c23, float &c24, float &c25, float &c26,
    float &c33, float &c34, float &c35, float &c36,
    float &c44, float &c45, float &c46,
    float &c55, float &c56, float &c66, float &rho3d)
{
  int nsg = (NG+1)*(NG+1)*(NG+1);

  float ari_rho = 0.0;
  float mu_har = 0.0, mu_ari = 0.0;
  float lam2mu_ari = 0.0, lam_ari = 0.0;
  float lam2mu_har = 0.0, lam_lam2mu = 0.0;
  float lam2_lam2mu = 0.0;

  for (int isg = 0; isg < nsg; isg++) {
    std::vector<float> var(3);
    AssignLayerMediaPara2Point(i, j, k,
        SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var);

    float vp=var[1], vs=var[2], rho=var[0];
    float mu = vs*vs*rho;
    float lam2mu = vp*vp*rho;
    float lam = lam2mu - 2.0*mu;

    lam2mu_har  += (1.0/lam2mu);
    lam2mu_ari  += lam2mu;
    lam_ari     += lam;
    lam_lam2mu  += (lam/lam2mu);
    lam2_lam2mu += (lam*lam/lam2mu);
    mu_har      += (1.0/mu);
    mu_ari      += mu;
    ari_rho += rho;
  }
  rho3d = ari_rho/nsg;
  lam2mu_har = nsg*1.0/lam2mu_har;
  lam2mu_ari /= nsg;
  lam_ari /= nsg;
  lam2_lam2mu /= nsg;
  lam_lam2mu /= nsg;
  mu_har = nsg*1.0/mu_har;
  mu_ari /= nsg;

  float Cz = lam2mu_ari - lam2_lam2mu + lam_lam2mu*lam_lam2mu*lam2mu_har;
  float Dz = lam_ari - lam2_lam2mu + lam_lam2mu*lam_lam2mu*lam2mu_har;
  float Bz = lam_lam2mu * lam2mu_har;
  float Az = lam2mu_har;

  BondTransform_CoordSystem(Cz, Dz, Bz, 0.0, 0.0, 0.0, 
                            Cz, Bz, 0.0, 0.0, 0.0, Az,
                            0.0, 0.0, 0.0, mu_har, 0.0, 0.0,
                            mu_har, 0.0, mu_ari, 
                            a, b, c, 
                            c11, c12, c13, c14, c15, c16, 
                            c22, c23, c24, c25, c26, c33,
                            c34, c35, c36, c44, c45, c46,
                            c55, c56, c66);
}

bool getInterfaceIntersection(
  int nx, int ny, float dx, float dy, 
  float x0, float y0, const float *elevation,
  const Point3 &v1, const Point3 &v2,
  Point3 &interP)
{
  float edge_maxx = std::max(v1.x, v2.x);
  float edge_minx = std::min(v1.x, v2.x);
  float edge_maxy = std::max(v1.y, v2.y);
  float edge_miny = std::min(v1.y, v2.y);

  int ix_min = (edge_minx-x0)/dx;
  int iy_min = (edge_miny-y0)/dy;
  int ix_max = (edge_maxx-x0)/dx;
  int iy_max = (edge_maxy-y0)/dy;
  if (ix_max == nx-1) ix_max -= 1;
  if (iy_max == ny-1) iy_max -= 1;

  bool hasInter = false;
  // find the intersection points Within the region
  for (int ix = ix_min; ix <= ix_max; ix++) {
    for (int iy = iy_min; iy <= iy_max; iy++) {
      float x0_interp = x0+ix*dx;
      float y0_interp = x0+iy*dy;
      float x1_interp = x0_interp+dx;
      float y1_interp = y0_interp+dy;
      size_t indx = iy*nx + ix;

      interP = PlaneBilinearIntersectSegment(x0_interp, y0_interp,
              x1_interp, y1_interp, elevation[indx], elevation[indx+1],
              elevation[indx+nx], elevation[indx+1+nx], v1, v2, hasInter);
      if (hasInter) {
        return true;
      } 
    }
  }

  return hasInter;
}


void getInersectList(int nx, int ny, float dx, float dy, float x0, float y0, 
                const float *elev, Mesh3 M, int edgeState, 
                std::set<Point3> &intersectionList)
{
    Point3 P;
    if (edgeState & 1)  // edge: v[0]--v[1] 
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[0], M.v[1], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 2)  // edge: v[1] -- v[2]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[1], M.v[2], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 4)  // edge: v[2] -- v[3]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[2], M.v[3], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 8)  // edge: v[0] -- v[3]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[0], M.v[3], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 16)  // edge: v[4] -- v[5]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[4], M.v[5], P);
      if (has_inter) intersectionList.insert(P);
    }

    if (edgeState & 32)  // edge: v[5] -- v[6]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[5], M.v[6], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 64)  // edge: v[6] -- v[7]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[6], M.v[7], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 128)  // edge: v[4] -- v[7]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[4], M.v[7], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 256)  // edge: v[0] -- v[4]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[0], M.v[4], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 512)  // edge: v[1] -- v[5]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[1], M.v[5], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 1024)  // edge: v[2] -- v[6]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[2], M.v[6], P);
      if (has_inter) intersectionList.insert(P);
    }
    if (edgeState & 2048)  // edge: v[3] -- v[7]
    {
      bool has_inter = getInterfaceIntersection(nx, ny, dx, dy, x0, y0, elev, M.v[3], M.v[7], P);
      if (has_inter) intersectionList.insert(P);
    }
}


//
// calculate orthorhombic representation val
//- for elasic isotropic
//  Kristek et al., 2017, An orthorhombic representation
//
void cal_ort_represent_val(
    int indx_i, int indx_j, int indx_k,
    inter_t &interfaces, Point3 *SubGrid,
    float &c11, float &c12, float &c13, 
    float &c22, float &c23, float &c33,
    float &c44, float &c55, float &c66, float &rho3d)
{
    int nsg = (NG+1)*(NG+1)*(NG+1);
    float *lam = new float[nsg];
    float *mu  = new float[nsg];
    float *lam2mu = new float[nsg];

    float ari_rho = 0.0;
    for (int isg = 0; isg < nsg; isg++) {
      std::vector<float> var(3);
      AssignLayerMediaPara2Point(indx_i, indx_j, indx_k,
          SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var);

      float vp=var[1], vs=var[2], rho=var[0];
      mu[isg]     = vs*vs*rho;
      lam2mu[isg] = vp*vp*rho;
      lam[isg]    = lam2mu[isg] - 2.0*mu[isg];
      ari_rho += rho;
    }
    rho3d = ari_rho/nsg;

    int siz_line = NG+1;
    int siz_slice = (NG+1) * siz_line;

    // calculate \Pi_x (c11) 
    float Pi_x = 0.0;
    for (int i = 0; i < NG+1; i++) {
      float tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;

      for (int j = 0; j < NG+1; j++) {
        for (int k = 0; k < NG+1; k++) {
          int indx = i + j * siz_line + k * siz_slice;
          tmp1 += (lam2mu[indx]-lam[indx]*lam[indx]/lam2mu[indx]);
          tmp2 += (lam[indx]/lam2mu[indx]);
          tmp3 += (1.0/lam2mu[indx]);
        }
      }

      tmp1 /= siz_slice;
      tmp2 /= siz_slice;
      tmp3 = (siz_slice*1.0)/tmp3;

      Pi_x += (1.0 / (tmp1 + tmp2*tmp2 * tmp3));
    }
    c11 = (siz_line*1.0)/Pi_x; 

    // calculate \Pi_y (c22)
    float Pi_y = 0.0;
    for (int j = 0; j <= NG; j++) {
      float tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
      for (int i = 0; i <= NG; i++) {
        for (int k = 0; k <= NG; k++) {
          int indx = i + j * siz_line + k * siz_slice;
          tmp1 += (lam2mu[indx]-lam[indx]*lam[indx]/lam2mu[indx]);
          tmp2 += (lam[indx]/lam2mu[indx]);
          tmp3 += (1.0/lam2mu[indx]);          
        }
      }
      tmp1 /= siz_slice;
      tmp2 /= siz_slice;
      tmp3 = (siz_slice*1.0)/tmp3;

      Pi_y += (1.0 / (tmp1 + tmp2*tmp2 * tmp3));      
    }
    c22 = (siz_line*1.0)/Pi_y;

    // calculate \Pi_z (c33)
    float Pi_z = 0.0;
    for (int k = 0; k <= NG; k++) {
      float tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
      for (int i = 0; i <= NG; i++) {
        for (int j = 0; j <= NG; j++) {
          int indx = i + j * siz_line + k * siz_slice;
          tmp1 += (lam2mu[indx]-lam[indx]*lam[indx]/lam2mu[indx]);
          tmp2 += (lam[indx]/lam2mu[indx]);
          tmp3 += (1.0/lam2mu[indx]);          
        }
      }
      tmp1 /= siz_slice;
      tmp2 /= siz_slice;
      tmp3 = (siz_slice*1.0)/tmp3;

      Pi_z += (1.0 / (tmp1 + tmp2*tmp2 * tmp3));      
    }
    c33 = (siz_line*1.0)/Pi_z;

    // \lambda_{xy}, c12;  c66
    float Lxy1 = 0.0, Lxy2 = 0.0;
    float muxy = 0.0;
    for (int i = 0; i <= NG; i++) {
      for (int j = 0; j <= NG; j++) {

        float tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
        float mu_ = 0.0, M = 0.0;
        for (int k = 0; k <= NG; k++) {
          int indx = i + j * siz_line + k * siz_slice;
          tmp1 += (lam2mu[indx] - lam[indx]*lam[indx]/lam2mu[indx]);
          tmp2 += (lam[indx]/lam2mu[indx]);
          tmp3 += (lam[indx] - lam[indx]*lam[indx]/lam2mu[indx]);
          M   += (1.0/lam2mu[indx]);
          mu_ += mu[indx]; 
        }
        tmp1 /= siz_line; 
        tmp2 /= siz_line;
        tmp3 /= siz_line;
        M = (siz_line * 1.0)/M;
        mu_ /= siz_line;

        muxy += (1.0/mu_);
        Lxy1 += (1.0/ (tmp1 + tmp2*tmp2*M));
        Lxy2 += ((tmp3+tmp2*tmp2*M) / (tmp1+tmp2*tmp2*M));
      }
    }
    Lxy1 = (1.0*siz_slice)/Lxy1;
    Lxy2 /= siz_slice;
    c12 = Lxy1 * Lxy2;
    c66 = (1.0*siz_slice)/muxy;

    // \lambda_{xz}, c13;  c55
    float Lxz1 = 0.0, Lxz2 = 0.0;
    float muxz = 0.0;
    for (int i = 0; i <= NG; i++) {
      for (int k = 0; k <= NG; k++) {

        float tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
        float mu_ = 0.0, M = 0.0;
        for (int j = 0; j <= NG; j++) {
          int indx = i + j * siz_line + k * siz_slice;

          tmp1 += (lam2mu[indx] - lam[indx]*lam[indx]/lam2mu[indx]);
          tmp2 += (lam[indx]/lam2mu[indx]);
          tmp3 += (lam[indx] - lam[indx]*lam[indx]/lam2mu[indx]);
          M   += (1.0/lam2mu[indx]);
          mu_ += mu[indx]; 
        }
        tmp1 /= siz_line; 
        tmp2 /= siz_line;
        tmp3 /= siz_line;
        M = (siz_line * 1.0)/M;
        mu_ /= siz_line;

        muxz += (1.0/mu_);
        Lxz1 += (1.0/ (tmp1 + tmp2*tmp2*M));
        Lxz2 += ((tmp3+tmp2*tmp2*M) / (tmp1+tmp2*tmp2*M));
      }
    }
    Lxz1 = (1.0*siz_slice)/Lxz1;
    Lxz2 /= siz_slice;
    c12 = Lxz1 * Lxz2;
    c55 = (1.0*siz_slice)/muxz;

    // \lambda_{yz}, c23;  c44
    float Lyz1 = 0.0, Lyz2 = 0.0;
    float muyz = 0.0;
    for (int k = 0; k <= NG; k++) {
      for (int j = 0; j <= NG; j++) {

        float tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
        float mu_ = 0.0, M = 0.0;
        for (int i = 0; i <= NG; i++) {
          int indx = i + j * siz_line + k * siz_slice;
          tmp1 += (lam2mu[indx] - lam[indx]*lam[indx]/lam2mu[indx]);
          tmp2 += (lam[indx]/lam2mu[indx]);
          tmp3 += (lam[indx] - lam[indx]*lam[indx]/lam2mu[indx]);
          M   += (1.0/lam2mu[indx]);
          mu_ += mu[indx]; 
        }
        tmp1 /= siz_line; 
        tmp2 /= siz_line;
        tmp3 /= siz_line;
        M = (siz_line * 1.0)/M;
        mu_ /= siz_line;

        muyz += (1.0/mu_);
        Lyz1 += (1.0/ (tmp1 + tmp2*tmp2*M));
        Lyz2 += ((tmp3+tmp2*tmp2*M) / (tmp1+tmp2*tmp2*M));
      }
    }
    Lyz1 = (1.0*siz_slice)/Lyz1;
    Lyz2 /= siz_slice;
    c12 = Lyz1 * Lyz2;
    c44 = (1.0*siz_slice)/muyz;


    if (lam    != nullptr) delete[] lam; 
    if (mu     != nullptr) delete[] mu;
    if (lam2mu != nullptr) delete[] lam2mu;  
}

