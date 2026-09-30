//==============================================================================
//  GridSplineBetaPolynomials.h
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
//  Definitions for the spline polynomials of the second kind (beta-s)       
//  in addition, various wrappers (in case they are needed)       
//       
//       
//  This algorithm was written by Cristian Lalescu.       
//  See  http://arxiv.org/abs/0905.3564 for details.       
//  The code here is slightly modified for use in Sartre.       
//       
//==============================================================================        
#ifndef GridSplineBetaPolynomials_h       
#define GridSplineBetaPolynomials_h       
       
inline double BETA_n1_q2_j1    (double x)  { return ((x>=0 && x<=1) ? (1-x) : 0); }       
inline double OMEGA_n1_q2      (double x)  { return ((x>=0) ? BETA_n1_q2_j1(x) : BETA_n1_q2_j1(-x)); }       
       
inline double BETA_n3_q4_j1    (double x)  { return    ((x*((4-3*x)*x+1))/2); }       
inline double BETA_n3_q4_j2    (double x)  { return    (((x-1)*x*x)/2); }       
inline double OMEGA_n3_q4_half (double x)  { return ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n3_q4_j1(1-x) : BETA_n3_q4_j2(2-x)) : 0); }       
inline double OMEGA_n3_q4      (double x)  { return ((x>=0) ? OMEGA_n3_q4_half(x) : OMEGA_n3_q4_half(-x) ); }       
       
inline double BETA_n5_q4_j1    (double x)  { return    ((x*(x*(x*(x*(6*x-15)+9)+1)+1))/2); }       
inline double BETA_n5_q4_j2    (double x)  { return    ((x*x*x*((5-2*x)*x-3))/2); }       
inline double OMEGA_n5_q4_half (double x)  { return ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n5_q4_j1(1-x) : BETA_n5_q4_j2(2-x)) : 0); }       
inline double OMEGA_n5_q4      (double x)  { return ((x>=0) ? OMEGA_n5_q4_half(x) : OMEGA_n5_q4_half(-x) ); }       
       
inline double BETA_n3_q6_j3    (double x)  { return    (((1-x)*x*x)/12); }       
inline double BETA_n3_q6_j2    (double x)  { return    ((x*(x*(7*x-6)-1))/12); }       
inline double BETA_n3_q6_j1    (double x)  { return    ((x*((5-4*x)*x+2))/3); }       
inline double OMEGA_n3_q6_half (double x)  { return ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n3_q6_j1(1-x) : BETA_n3_q6_j2(2-x)) : BETA_n3_q6_j3(3-x)) : 0); }       
inline double OMEGA_n3_q6      (double x)  { return ((x>=0) ? OMEGA_n3_q6_half(x) : OMEGA_n3_q6_half(-x)); }       
       
inline double BETA_n5_q6_j3    (double x)  { return    ((x*x*x*(x*(5*x-12)+7))/24); }       
inline double BETA_n5_q6_j2    (double x)  { return    ((x*(x*(x*((61-25*x)*x-33)-1)-2))/24); }       
inline double BETA_n5_q6_j1    (double x)  { return    ((x*(x*(x*(x*(25*x-62)+33)+8)+8))/12); }       
inline double OMEGA_n5_q6_half (double x)  { return ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n5_q6_j1(1-x) : BETA_n5_q6_j2(2-x)) : BETA_n5_q6_j3(3-x)) : 0); }       
inline double OMEGA_n5_q6      (double x)  { return ((x>=0) ? OMEGA_n5_q6_half(x) : OMEGA_n5_q6_half(-x)); }       
       
inline double BETA_n7_q6_j3    (double x)  { return    ((x*x*x*x*(x*((49-14*x)*x-58)+23))/24); }       
inline double BETA_n7_q6_j2    (double x)  { return    ((x*(x*(x*(x*(x*(x*(70*x-245)+290)-114)+2)-1)-2))/24); }       
inline double BETA_n7_q6_j1    (double x)  { return    ((x*(x*(x*(x*(x*((245-70*x)*x-290)+113)-2)+8)+8))/12); }       
inline double OMEGA_n7_q6_half (double x)  { return ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n7_q6_j1(1-x) : BETA_n7_q6_j2(2-x)) : BETA_n7_q6_j3(3-x)) : 0); }       
inline double OMEGA_n7_q6      (double x)  { return ((x>=0) ? OMEGA_n7_q6_half(x) : OMEGA_n7_q6_half(-x)); }       
       
inline double BETA_n9_q6_j3    (double x)  { return    ((x*x*x*x*x*(x*(x*(x*(46*x-207)+354)-273)+80))/24); }       
inline double BETA_n9_q6_j2    (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*((1035-230*x)*x-1770)+1365)-400)+1)+2)-1)-2))/24); }       
inline double BETA_n9_q6_j1    (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(230*x-1035)+1770)-1365)+400)-2)-2)+8)+8))/12); }       
inline double OMEGA_n9_q6_half (double x)  { return ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n9_q6_j1(1-x) : BETA_n9_q6_j2(2-x)) : BETA_n9_q6_j3(3-x)) : 0); }       
inline double OMEGA_n9_q6      (double x)  { return ((x>=0) ? OMEGA_n9_q6_half(x) : OMEGA_n9_q6_half(-x)); }       
       
inline double BETA_n3_q8_j4    (double x)  { return    (((x-1)*x*x)/60); }       
inline double BETA_n3_q8_j3    (double x)  { return    ((x*((7-8*x)*x+1))/60); }       
inline double BETA_n3_q8_j2    (double x)  { return    ((x*(x*(12*x-9)-3))/20); }       
inline double BETA_n3_q8_j1    (double x)  { return    ((x*((6-5*x)*x+3))/4); }       
inline double OMEGA_n3_q8_half (double x)  { return ((x>=0 && x<=4) ? ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n3_q8_j1(1-x) : BETA_n3_q8_j2(2-x)) : BETA_n3_q8_j3(3-x)) : BETA_n3_q8_j4(4-x)):0); }       
inline double OMEGA_n3_q8      (double x)  { return ((x>=0) ? OMEGA_n3_q8_half(x) : OMEGA_n3_q8_half(-x)); }       
       
inline double BETA_n5_q8_j4    (double x)  { return    ((x*x*x*((19-8*x)*x-11))/180); }       
inline double BETA_n5_q8_j3    (double x)  { return    ((x*(x*(x*(x*(115*x-270)+147)+2)+6))/360); }       
inline double BETA_n5_q8_j2    (double x)  { return    ((x*(x*(x*((93-39*x)*x-45)-3)-6))/40); }       
inline double BETA_n5_q8_j1    (double x)  { return    ((x*(x*(x*(x*(59*x-145)+68)+27)+27))/36); }       
inline double OMEGA_n5_q8_half (double x)  { return ((x>=0 && x<=4) ? ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n5_q8_j1(1-x) : BETA_n5_q8_j2(2-x)) : BETA_n5_q8_j3(3-x)) : BETA_n5_q8_j4(4-x)):0); }       
inline double OMEGA_n5_q8      (double x)  { return ((x>=0) ? OMEGA_n5_q8_half(x) : OMEGA_n5_q8_half(-x)); }       
       
inline double BETA_n7_q8_j4    (double x)  { return    ((x*x*x*x*(x*(x*(89*x-311)+367)-145))/720); }       
inline double BETA_n7_q8_j3    (double x)  { return    ((x*(x*(x*(x*(x*((2178-623*x)*x-2566)+1010)-15)+4)+12))/720); }       
inline double BETA_n7_q8_j2    (double x)  { return    ((x*(x*(x*(x*(x*(x*(623*x-2179)+2565)-995)+40)-18)-36))/240); }       
inline double BETA_n7_q8_j1    (double x)  { return    ((x*(x*(x*(x*(x*((2180-623*x)*x-2566)+976)-39)+108)+108))/144); }       
inline double OMEGA_n7_q8_half (double x)  { return ((x>=0 && x<=4) ? ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n7_q8_j1(1-x) : BETA_n7_q8_j2(2-x)) : BETA_n7_q8_j3(3-x)) : BETA_n7_q8_j4(4-x)):0); }       
inline double OMEGA_n7_q8      (double x)  { return ((x>=0) ? OMEGA_n7_q8_half(x) : OMEGA_n7_q8_half(-x)); }       
       
inline double BETA_n9_q8_j4    (double x)  { return    ((x*x*x*x*x*(x*(x*((1305-290*x)*x-2231)+1719)-503))/720); }       
inline double BETA_n9_q8_j3    (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(2030*x-9135)+15617)-12032)+3524)-5)-15)+4)+12))/720); }       
inline double BETA_n9_q8_j2    (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*((9135-2030*x)*x-15617)+12031)-3525)+20)+40)-18)-36))/240); }       
inline double BETA_n9_q8_j1    (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(2030*x-9135)+15617)-12030)+3524)-39)-39)+108)+108))/144); }       
inline double OMEGA_n9_q8_half (double x)  { return ((x>=0 && x<=4) ? ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n9_q8_j1(1-x) : BETA_n9_q8_j2(2-x)) : BETA_n9_q8_j3(3-x)) : BETA_n9_q8_j4(4-x)):0); }       
inline double OMEGA_n9_q8      (double x)  { return ((x>=0) ? OMEGA_n9_q8_half(x) : OMEGA_n9_q8_half(-x)); }       
       
inline double BETA_n11_q8_j4   (double x)  { return    ((x*x*x*x*x*x*(x*(x*(x*(x*(1006*x-5533)+12285)-13785)+7829)-1802))/720); }       
inline double BETA_n11_q8_j3   (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(x*((38731-7042*x)*x-85995)+96495)-54803)+12615)+3)-5)-15)+4)+12))/720); }       
inline double BETA_n11_q8_j2   (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(7042*x-38731)+85995)-96495)+54803)-12616)-4)+20)+40)-18)-36))/240); }       
inline double BETA_n11_q8_j1   (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(x*((38731-7042*x)*x-85995)+96495)-54803)+12617)+3)-39)-39)+108)+108))/144); }       
inline double OMEGA_n11_q8_half(double x)  { return ((x>=0 && x<=4) ? ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n11_q8_j1(1-x) : BETA_n11_q8_j2(2-x)) : BETA_n11_q8_j3(3-x)) : BETA_n11_q8_j4(4-x)):0); }       
inline double OMEGA_n11_q8     (double x)  { return ((x>=0) ? OMEGA_n11_q8_half(x) : OMEGA_n11_q8_half(-x)); }       
       
inline double BETA_n13_q8_j4   (double x)  { return    ((x*x*x*x*x*x*x*(x*(x*(x*(x*((23426-3604*x)*x-63866)+93577)-77815)+34869)-6587))/720); }       
inline double BETA_n13_q8_j3   (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(25228*x-163982)+447062)-655039)+544705)-244083)+46109)+1)+3)-5)-15)+4)+12))/720); }       
inline double BETA_n13_q8_j2   (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(x*((163982-25228*x)*x-447062)+655039)-544705)+244083)-46109)-2)-4)+20)+40)-18)-36))/240); }       
inline double BETA_n13_q8_j1   (double x)  { return    ((x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(x*(25228*x-163982)+447062)-655039)+544705)-244083)+46109)+3)+3)-39)-39)+108)+108))/144); }       
inline double OMEGA_n13_q8_half(double x)  { return ((x>=0 && x<=4) ? ((x>=0 && x<=3) ? ((x>=0 && x<=2) ? ((x>=0 && x<=1) ? BETA_n13_q8_j1(1-x) : BETA_n13_q8_j2(2-x)) : BETA_n13_q8_j3(3-x)) : BETA_n13_q8_j4(4-x)):0); }       
inline double OMEGA_n13_q8     (double x)  { return ((x>=0) ? OMEGA_n13_q8_half(x) : OMEGA_n13_q8_half(-x)); }       
       
/********************************/       
// (probably inefficient) wrappers       
       
inline double BETA_n3_q4 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -1: return BETA_n3_q4_j2(1-x);       
		case  0: return BETA_n3_q4_j1(1-x);       
		case  1: return BETA_n3_q4_j1(x);       
		case  2: return BETA_n3_q4_j2(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n5_q4 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -1: return BETA_n5_q4_j2(1-x);       
		case  0: return BETA_n5_q4_j1(1-x);       
		case  1: return BETA_n5_q4_j1(x);       
		case  2: return BETA_n5_q4_j2(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n3_q6 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -2: return BETA_n3_q6_j3(1-x);       
		case -1: return BETA_n3_q6_j2(1-x);       
		case  0: return BETA_n3_q6_j1(1-x);       
		case  1: return BETA_n3_q6_j1(x);       
		case  2: return BETA_n3_q6_j2(x);       
		case  3: return BETA_n3_q6_j3(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n5_q6 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -2: return BETA_n5_q6_j3(1-x);       
		case -1: return BETA_n5_q6_j2(1-x);       
		case  0: return BETA_n5_q6_j1(1-x);       
		case  1: return BETA_n5_q6_j1(x);       
		case  2: return BETA_n5_q6_j2(x);       
		case  3: return BETA_n5_q6_j3(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n7_q6 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -2: return BETA_n7_q6_j3(1-x);       
		case -1: return BETA_n7_q6_j2(1-x);       
		case  0: return BETA_n7_q6_j1(1-x);       
		case  1: return BETA_n7_q6_j1(x);       
		case  2: return BETA_n7_q6_j2(x);       
		case  3: return BETA_n7_q6_j3(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n9_q6 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -2: return BETA_n9_q6_j3(1-x);       
		case -1: return BETA_n9_q6_j2(1-x);       
		case  0: return BETA_n9_q6_j1(1-x);       
		case  1: return BETA_n9_q6_j1(x);       
		case  2: return BETA_n9_q6_j2(x);       
		case  3: return BETA_n9_q6_j3(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n3_q8 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -3: return BETA_n3_q8_j4(1-x);       
		case -2: return BETA_n3_q8_j3(1-x);       
		case -1: return BETA_n3_q8_j2(1-x);       
		case  0: return BETA_n3_q8_j1(1-x);       
		case  1: return BETA_n3_q8_j1(x);       
		case  2: return BETA_n3_q8_j2(x);       
		case  3: return BETA_n3_q8_j3(x);       
		case  4: return BETA_n3_q8_j4(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n5_q8 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -3: return BETA_n5_q8_j4(1-x);       
		case -2: return BETA_n5_q8_j3(1-x);       
		case -1: return BETA_n5_q8_j2(1-x);       
		case  0: return BETA_n5_q8_j1(1-x);       
		case  1: return BETA_n5_q8_j1(x);       
		case  2: return BETA_n5_q8_j2(x);       
		case  3: return BETA_n5_q8_j3(x);       
		case  4: return BETA_n5_q8_j4(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n7_q8 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -3: return BETA_n7_q8_j4(1-x);       
		case -2: return BETA_n7_q8_j3(1-x);       
		case -1: return BETA_n7_q8_j2(1-x);       
		case  0: return BETA_n7_q8_j1(1-x);       
		case  1: return BETA_n7_q8_j1(x);       
		case  2: return BETA_n7_q8_j2(x);       
		case  3: return BETA_n7_q8_j3(x);       
		case  4: return BETA_n7_q8_j4(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n9_q8 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -3: return BETA_n9_q8_j4(1-x);       
		case -2: return BETA_n9_q8_j3(1-x);       
		case -1: return BETA_n9_q8_j2(1-x);       
		case  0: return BETA_n9_q8_j1(1-x);       
		case  1: return BETA_n9_q8_j1(x);       
		case  2: return BETA_n9_q8_j2(x);       
		case  3: return BETA_n9_q8_j3(x);       
		case  4: return BETA_n9_q8_j4(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n11_q8 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -3: return BETA_n11_q8_j4(1-x);       
		case -2: return BETA_n11_q8_j3(1-x);       
		case -1: return BETA_n11_q8_j2(1-x);       
		case  0: return BETA_n11_q8_j1(1-x);       
		case  1: return BETA_n11_q8_j1(x);       
		case  2: return BETA_n11_q8_j2(x);       
		case  3: return BETA_n11_q8_j3(x);       
		case  4: return BETA_n11_q8_j4(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_n13_q8 (int grid_point, double x)       
{       
	switch(grid_point)       
	{       
		case -3: return BETA_n13_q8_j4(1-x);       
		case -2: return BETA_n13_q8_j3(1-x);       
		case -1: return BETA_n13_q8_j2(1-x);       
		case  0: return BETA_n13_q8_j1(1-x);       
		case  1: return BETA_n13_q8_j1(x);       
		case  2: return BETA_n13_q8_j2(x);       
		case  3: return BETA_n13_q8_j3(x);       
		case  4: return BETA_n13_q8_j4(x);       
		default: return 0;       
	}       
}       
       
inline double BETA_q4 (int spline_order, int grid_point, double x)       
{       
	switch(spline_order)       
	{       
		case 3: return BETA_n3_q4(grid_point,x);       
		case 5: return BETA_n5_q4(grid_point,x);       
		default: return 0;       
	}       
}       
       
inline double BETA_q6 (int spline_order, int grid_point, double x)       
{       
	switch(spline_order)       
	{       
		case 3: return BETA_n3_q6(grid_point,x);       
		case 5: return BETA_n5_q6(grid_point,x);       
		case 7: return BETA_n7_q6(grid_point,x);       
		case 9: return BETA_n9_q6(grid_point,x);       
		default: return 0;       
	}       
}       
       
inline double BETA_q8 (int spline_order, int grid_point, double x)       
{       
	switch(spline_order)       
	{       
		case  3: return  BETA_n3_q8(grid_point,x);       
		case  5: return  BETA_n5_q8(grid_point,x);       
		case  7: return  BETA_n7_q8(grid_point,x);       
		case  9: return  BETA_n9_q8(grid_point,x);       
		case 11: return BETA_n11_q8(grid_point,x);       
		case 13: return BETA_n13_q8(grid_point,x);       
		default: return 0;       
	}       
}       
       
inline double BETA (int spline_order, int grid_size, int grid_point, double x)       
{       
	switch(grid_size)       
	{       
		case 4: return BETA_q4(spline_order, grid_point, x);       
		case 6: return BETA_q6(spline_order, grid_point, x);       
		case 8: return BETA_q8(spline_order, grid_point, x);       
		default: return 0;       
	}       
}       
       
#endif//__BETA_POLYNOMIALS_H__       
       
