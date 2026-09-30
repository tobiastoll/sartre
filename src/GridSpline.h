//==============================================================================
//  GridSpline.h
//
//  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich 
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:35:41 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
//       
//  Implementation of grid splines (a type of spline interpolations)       
//       
//  This algorithm was written by Cristian Lalescu.       
//  See  http://arxiv.org/abs/0905.3564 for details.       
//  The code here is slightly modified for use in Sartre.       
//       
//==============================================================================        
#ifndef GridSpline_h        
#define GridSpline_h        
#include "GridSplineBetaPolynomials.h"       
#include "GridSplineInterpolateSums.h"       
       
template <unsigned int q> class local_scal_1D       
{       
	private:       
		double val[q];       
		int g;       
	public:       
		local_scal_1D(){g = (q-2)/2;}       
		~local_scal_1D(){}       
		inline double &operator()(int grid_point)       
		{       
			return val[(grid_point+g)%q];       
		}       
};       
       
template <unsigned int q> class local_scal_3D       
{       
	private:       
		double val[q][q][q];       
		int g;       
	public:       
		local_scal_3D(){g = (q-2)/2;}       
		~local_scal_3D(){}       
		inline double &operator()(int grid_pointx, int grid_pointy, int grid_pointz)       
		{       
			return val[(grid_pointx+g)%q][(grid_pointy+g)%q][(grid_pointz+g)%q];       
		}       
};       
       
inline void beta_n3(local_scal_1D<4> &bx, double x)       
{       
	bx(-1) = BETA_n3_q4_j2(1-x);       
	bx( 0) = BETA_n3_q4_j1(1-x);       
	bx( 1) = BETA_n3_q4_j1(x);       
	bx( 2) = BETA_n3_q4_j2(x);       
}       
       
inline void beta_n5(local_scal_1D<4> &bx, double x)       
{       
	bx(-1) = BETA_n5_q4_j2(1-x);       
	bx( 0) = BETA_n5_q4_j1(1-x);       
	bx( 1) = BETA_n5_q4_j1(x);       
	bx( 2) = BETA_n5_q4_j2(x);       
}       
       
inline void beta_n3(local_scal_1D<6> &bx, double x)       
{       
	bx(-2) = BETA_n3_q6_j3(1-x);       
	bx(-1) = BETA_n3_q6_j2(1-x);       
	bx( 0) = BETA_n3_q6_j1(1-x);       
	bx( 1) = BETA_n3_q6_j1(x);       
	bx( 2) = BETA_n3_q6_j2(x);       
	bx( 3) = BETA_n3_q6_j3(x);       
}       
       
inline void beta_n5(local_scal_1D<6> &bx, double x)       
{       
	bx(-2) = BETA_n5_q6_j3(1-x);       
	bx(-1) = BETA_n5_q6_j2(1-x);       
	bx( 0) = BETA_n5_q6_j1(1-x);       
	bx( 1) = BETA_n5_q6_j1(x);       
	bx( 2) = BETA_n5_q6_j2(x);       
	bx( 3) = BETA_n5_q6_j3(x);       
}       
       
inline void beta_n7(local_scal_1D<6> &bx, double x)       
{       
	bx(-2) = BETA_n7_q6_j3(1-x);       
	bx(-1) = BETA_n7_q6_j2(1-x);       
	bx( 0) = BETA_n7_q6_j1(1-x);       
	bx( 1) = BETA_n7_q6_j1(x);       
	bx( 2) = BETA_n7_q6_j2(x);       
	bx( 3) = BETA_n7_q6_j3(x);       
}       
       
inline void beta_n9(local_scal_1D<6> &bx, double x)       
{       
	bx(-2) = BETA_n9_q6_j3(1-x);       
	bx(-1) = BETA_n9_q6_j2(1-x);       
	bx( 0) = BETA_n9_q6_j1(1-x);       
	bx( 1) = BETA_n9_q6_j1(x);       
	bx( 2) = BETA_n9_q6_j2(x);       
	bx( 3) = BETA_n9_q6_j3(x);       
}       
       
inline void beta_n3(local_scal_1D<8> &bx, double x)       
{       
	bx(-3) = BETA_n3_q8_j4(1-x);       
	bx(-2) = BETA_n3_q8_j3(1-x);       
	bx(-1) = BETA_n3_q8_j2(1-x);       
	bx( 0) = BETA_n3_q8_j1(1-x);       
	bx( 1) = BETA_n3_q8_j1(x);       
	bx( 2) = BETA_n3_q8_j2(x);       
	bx( 3) = BETA_n3_q8_j3(x);       
	bx( 4) = BETA_n3_q8_j4(x);       
}       
       
inline void beta_n5(local_scal_1D<8> &bx, double x)       
{       
	bx(-3) = BETA_n5_q8_j4(1-x);       
	bx(-2) = BETA_n5_q8_j3(1-x);       
	bx(-1) = BETA_n5_q8_j2(1-x);       
	bx( 0) = BETA_n5_q8_j1(1-x);       
	bx( 1) = BETA_n5_q8_j1(x);       
	bx( 2) = BETA_n5_q8_j2(x);       
	bx( 3) = BETA_n5_q8_j3(x);       
	bx( 4) = BETA_n5_q8_j4(x);       
}       
       
inline void beta_n7(local_scal_1D<8> &bx, double x)       
{       
	bx(-3) = BETA_n7_q8_j4(1-x);       
	bx(-2) = BETA_n7_q8_j3(1-x);       
	bx(-1) = BETA_n7_q8_j2(1-x);       
	bx( 0) = BETA_n7_q8_j1(1-x);       
	bx( 1) = BETA_n7_q8_j1(x);       
	bx( 2) = BETA_n7_q8_j2(x);       
	bx( 3) = BETA_n7_q8_j3(x);       
	bx( 4) = BETA_n7_q8_j4(x);       
}       
       
inline void beta_n9(local_scal_1D<8> &bx, double x)       
{       
	bx(-3) = BETA_n9_q8_j4(1-x);       
	bx(-2) = BETA_n9_q8_j3(1-x);       
	bx(-1) = BETA_n9_q8_j2(1-x);       
	bx( 0) = BETA_n9_q8_j1(1-x);       
	bx( 1) = BETA_n9_q8_j1(x);       
	bx( 2) = BETA_n9_q8_j2(x);       
	bx( 3) = BETA_n9_q8_j3(x);       
	bx( 4) = BETA_n9_q8_j4(x);       
}       
       
inline void beta_n11(local_scal_1D<8> &bx, double x)       
{       
	bx(-3) = BETA_n11_q8_j4(1-x);       
	bx(-2) = BETA_n11_q8_j3(1-x);       
	bx(-1) = BETA_n11_q8_j2(1-x);       
	bx( 0) = BETA_n11_q8_j1(1-x);       
	bx( 1) = BETA_n11_q8_j1(x);       
	bx( 2) = BETA_n11_q8_j2(x);       
	bx( 3) = BETA_n11_q8_j3(x);       
	bx( 4) = BETA_n11_q8_j4(x);       
}       
       
inline void beta_n13(local_scal_1D<8> &bx, double x)       
{       
	bx(-3) = BETA_n13_q8_j4(1-x);       
	bx(-2) = BETA_n13_q8_j3(1-x);       
	bx(-1) = BETA_n13_q8_j2(1-x);       
	bx( 0) = BETA_n13_q8_j1(1-x);       
	bx( 1) = BETA_n13_q8_j1(x);       
	bx( 2) = BETA_n13_q8_j2(x);       
	bx( 3) = BETA_n13_q8_j3(x);       
	bx( 4) = BETA_n13_q8_j4(x);       
}       
       
double grid_spline(int spline_order, local_scal_3D<4> &f, double x, double y, double z)       
{       
	local_scal_1D<4> bx, by, bz;       
	switch(spline_order)       
	{       
		case 3:       
			beta_n3(bx,x);       
			beta_n3(by,y);       
			beta_n3(bz,z);       
			break;       
		case 5:       
			beta_n5(bx,x);       
			beta_n5(by,y);       
			beta_n5(bz,z);       
			break;       
	}       
	return SUM_3D_Q4_unsafe_macro(f,bx,by,bz);       
}       
       
double grid_spline(int spline_order, local_scal_3D<6> &f, double x, double y, double z)       
{       
	local_scal_1D<6> bx, by, bz;       
	switch(spline_order)       
	{       
		case 3:       
			beta_n3(bx,x);       
			beta_n3(by,y);       
			beta_n3(bz,z);       
			break;       
		case 5:       
			beta_n5(bx,x);       
			beta_n5(by,y);       
			beta_n5(bz,z);       
			break;       
		case 7:       
			beta_n7(bx,x);       
			beta_n7(by,y);       
			beta_n7(bz,z);       
			break;       
		case 9:       
			beta_n9(bx,x);       
			beta_n9(by,y);       
			beta_n9(bz,z);       
			break;       
	}       
	return SUM_3D_Q6_unsafe_macro(f,bx,by,bz);       
}       
       
double grid_spline(int spline_order, local_scal_3D<8> &f, double x, double y, double z)       
{       
	local_scal_1D<8> bx, by, bz;       
	switch(spline_order)       
	{       
		case 3:       
			beta_n3(bx,x);       
			beta_n3(by,y);       
			beta_n3(bz,z);       
			break;       
		case 5:       
			beta_n5(bx,x);       
			beta_n5(by,y);       
			beta_n5(bz,z);       
			break;       
		case 7:       
			beta_n7(bx,x);       
			beta_n7(by,y);       
			beta_n7(bz,z);       
			break;       
		case 9:       
			beta_n9(bx,x);       
			beta_n9(by,y);       
			beta_n9(bz,z);       
			break;       
		case 11:       
			beta_n11(bx,x);       
			beta_n11(by,y);       
			beta_n11(bz,z);       
			break;       
		case 13:       
			beta_n13(bx,x);       
			beta_n13(by,y);       
			beta_n13(bz,z);       
			break;       
	}       
	return SUM_3D_Q8_unsafe_macro(f,bx,by,bz);       
}       
       
#endif       
       
