/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#ifndef	__DISCRETIZE_H
#define	__DISCRETIZE_H

#include "jz_protein.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
For best performance, computation in VOTE uses fixed pointer arithmetics.
Distances/coordinates are discretized as:
	D=(int)d*128.0, or 1A = 128 unit
Angles are disctetized as:
	A=a*3.259493234522017e+02, or 1024 = PI

In a unit vector, the max value in each dimension is 32768, or,
the sum of squares equals to 1073741824. To calculate the angle between
two unit vectors, first calculate the dot of the two vector, right shift
19 bits, then look up table vote_res_acos_lut[].
-----------------------------------------------------------------------------*/

#define	VOTE_DISCRETIZE_COORD_FACTOR	128.0			/* 128  = 1A */
#define	VOTE_DISCRETIZE_ANGLE_FACTOR	3.259493234522017e+02	/* 1024 = PI */
#define	VOTE_DISCRETIZE_UNIT_VECTOR	32768			/* 32768= 1  */
#define	VOTE_ACOS_SHIFT			19
#define	VOTE_ACOS_LUT_SHIFT		2052

#define	VOTE_UNIT_VECTOR(V,d)						\
do{\
	JZ_VEC3_DOT((d),(V),(V));\
	(d)=(int)(sqrt(d+1));\
	JZ_VEC3_MUL((V),(V),VOTE_DISCRETIZE_UNIT_VECTOR);\
	JZ_VEC3_DIV((V),(V),(d));\
}while(0)	


/*-----------------------------------------------------------------------------
Convert residues in a protein into integer points to speed up calculation.
By default, the discretization uses 128 level per A. This is good if
square of distance difference is calculated.
-----------------------------------------------------------------------------*/

void	vote_discretize_residues(int *ri,jz_protein_residue *reses,int N);

#ifdef	__cplusplus
}
#endif

#endif

