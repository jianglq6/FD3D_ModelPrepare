/******************************************************************************
 *
 * 3D geometry functions.
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
 *******************************************************************************/


#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include "media_geometry3d.hpp"
//using namespace std;

int dcmp(float x) {
	if (fabs(x) < EPS) 
		return 0;
	else
		return x < 0?-1:1;
}

Vector3 operator + (const Vector3 &A, const Vector3 &B) {
	return Vector3(A.x+B.x, A.y+B.y, A.z+B.z);
}

Vector3 operator - (const Vector3 &A, const Vector3 &B) {
	return Vector3(A.x-B.x, A.y-B.y, A.z-B.z);
}

Vector3 operator * (const Vector3 &A, float p) {
	return Vector3(A.x*p, A.y*p, A.z*p);
}

Vector3 operator / (const Vector3 &A, float p) {
	return Vector3(A.x/p, A.y/p, A.z/p);
}

float Dot(const Vector3 &A, const Vector3 &B) {
	return A.x*B.x + A.y*B.y + A.z*B.z;
}


float Length(const Vector3 &A) {
	return sqrt(Dot(A, A));
}

float DistanceP2P(const Point3 &A, const Point3 &B) {
	return Length(B-A);	
}

float Angle(const Vector3 &A, const Vector3 &B) {
	return acos(Dot(A, B)/Length(A)/Length(B));
}

// directional cosine of two vectors
float CosinDir(const Vector3 &A, const Vector3 &B) {
  return Dot(A, B)/Length(A)/Length(B);
}

float DistanceToPlane(const Point3 &p, const Point3 &p0, const Vector3 &n) {
	return fabs(Dot(p-p0, n));
}

Point3 GetPlaneProjection(const Point3 &p, const Point3 &p0, const Vector3 &n) {
	return p - n*Dot(p-p0, n);
}

Vector3 Cross(const Vector3 &A, const Vector3 &B) {
	return Vector3(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x);
}

float Area2(const Point3 &A, const Point3 &B, const Point3 &C) {
	return Length(Cross(B-A, C-A));
}

//==== for plane
Plane::Plane(const Plane &P) {
  a = P.a;
  b = P.b;
  c = P.c;
  d = P.d;
}

Plane::Plane(float a_, float b_, float c_, float d_) {
  a = a_;
  b = b_;
  c = c_;
  d = d_;
}

Plane::Plane(const Vector3 &NormalizedNormal, float d_)
{
  a = NormalizedNormal.x;
  b = NormalizedNormal.y;
  c = NormalizedNormal.z;
  d = d_;
}

Plane Plane::Normalize()
{
  Plane Result;
  float Distance = sqrtf(a * a + b * b + c * c);
  Result.a = a / Distance;
  Result.b = b / Distance;
  Result.c = c / Distance;
  Result.d = d / Distance;
  return Result;
}

Point3 Plane::IntersectSegment(const Point3 &V1, const Point3 &V2, bool &isInter) const
{
  isInter = true;
  Vector3 Diff = V1 - V2;
  float Denominator = a * Diff.x + b * Diff.y + c * Diff.z;
  if(Denominator == 0.0f) {
    isInter = false;
    return V1;
  }
  float u = (a * V1.x + b * V1.y + c * V1.z + d) / Denominator;

  if (u > 1 || u < 0) isInter = false;

  return (V1 + (V2 - V1)*u);
}

float Plane::DotNormal(const Plane &P, const Vector3 &V)
{
  return P.a * V.x + P.b * V.y + P.c * V.z;
}


/*
bool PointInTri(Point3 P, Point3 P0, Point3 P1, Point3 P2) {
	float area1 = Area2(P, P0, P1);
	float area2 = Area2(P, P1, P2);
	float area3 = Area2(P, P2, P0);
	return dcmp(area1+area2+area3 - Area2(P0, P1, P2)) == 0;
}

float Volume6(Point3 A, Point3 B, Point3 C, Point3 D) {
	return Dot(D-A, Cross(B-A, C-A));
}
*/

/* 
 * Check whether the point is in the polyhedron.
 * Note: The hexahedron must be convex!
 */
bool isPointInPolyhedron(const Point3 &p, const std::vector<Face> &fs) {
	for (Face const &f:fs) {
		Vector3 p2f = f.v[0]-p;
		Vector3 A = f.normal();
		float sign = Dot(p2f, f.normal());
		sign /= Length(p2f); // normalization, for stability

		constexpr float bound = -1e-15;
		if(sign < 0.0) return false;
	}
	return true;
}

/*
 * Input: vx, vy, vz are the EIGHT vertexes of the hexahedron 
 *
 *    ↑ +z       4----6
 *    |         /|   /|
 *             / 0--/-2
 *            5----7 /
 *            |/   |/
 *            1----3
 *
 *
 */
bool isPointInHexahedron(float px, float py, float pz,
                         float *vx, float *vy, float *vz) 
{

	Point3 P(px, py, pz); 

	/* 
	 * Just for cgfd3D, in which the grid mesh maybe not a hexahedron,
	 */
	std::vector<Face> hexa{
		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[4], vy[4], vz[4]}, Point3{vx[6], vy[6], vz[6]}} }, 
//		Face{ {Point3{vx[6], vy[6], vz[6]}, Point3{vx[2], vy[2], vz[2]}, Point3{vx[0], vy[0], vz[0]}} }, // back

		Face{ {Point3{vx[7], vy[7], vz[7]}, Point3{vx[5], vy[5], vz[5]}, Point3{vx[1], vy[1], vz[1]}} },
//		Face{ {Point3{vx[1], vy[1], vz[1]}, Point3{vx[3], vy[3], vz[3]}, Point3{vx[7], vy[7], vz[7]}} }, // front

		Face{ {Point3{vx[5], vy[5], vz[5]}, Point3{vx[4], vy[4], vz[4]}, Point3{vx[0], vy[0], vz[0]}} }, 
//		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[1], vy[1], vz[1]}, Point3{vx[5], vy[5], vz[5]}} }, // left

		Face{ {Point3{vx[2], vy[2], vz[2]}, Point3{vx[6], vy[6], vz[6]}, Point3{vx[7], vy[7], vz[7]}} }, 
//		Face{ {Point3{vx[7], vy[7], vz[7]}, Point3{vx[3], vy[3], vz[3]}, Point3{vx[2], vy[2], vz[2]}} }, // right

		Face{ {Point3{vx[4], vy[4], vz[4]}, Point3{vx[5], vy[5], vz[5]}, Point3{vx[7], vy[7], vz[7]}} }, 
//		Face{ {Point3{vx[7], vy[7], vz[7]}, Point3{vx[6], vy[6], vz[6]}, Point3{vx[4], vy[4], vz[4]}} }, // top

//		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[2], vy[2], vz[2]}, Point3{vx[3], vy[3], vz[3]}} }, 
		Face{ {Point3{vx[3], vy[3], vz[3]}, Point3{vx[1], vy[1], vz[1]}, Point3{vx[0], vy[0], vz[0]}} }, // bottom
	};

    return isPointInPolyhedron(P, hexa);
}

/*
bool SegmentQuadIntersection(
  Point3 &p0, Point3 &p1, Point3 &p2, Point3 &p3,    // for quadrilateral in 3D
  Point3 &v0, Point3 &v1,    // for segment
  )
{



}

*/
