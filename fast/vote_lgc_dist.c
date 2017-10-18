/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>

#include "misc.h"

/*-----------------------------------------------------------------------------
Calculate LGC_DIST: to get LGC. This is seperated to be included in the
pre-processing phase of the search program.

D[N+N] is an array of integer.
D[k+k] is the distance between P[i-2] and P[i+1]
D[k+k+1] is the distance between P[i-2] and P[i+2]
-----------------------------------------------------------------------------*/

void	vote_lgc_dist(int *D,int *pt,int N)
{
	int		i,d;

	JZ_VEC3_SUBSQ(d,pt,pt+9);
	D[4]=(int)(sqrt(d));
	JZ_VEC3_SUBSQ(d,pt,pt+12);
	D[5]=(int)(sqrt(d));
	D[0]=D[2]=D[4];
	D[1]=D[3]=D[5];

	for(pt+=3,N-=2,i=3;i<N;i++,pt+=3)
	{
		JZ_VEC3_SUBSQ(d,pt,pt+9);
		D[i+i]=(int)(sqrt(d));
		JZ_VEC3_SUBSQ(d,pt,pt+12);
		D[i+i+1]=(int)(sqrt(d));
	}
	JZ_VEC3_SUBSQ(d,pt,pt+9);
	D[i+i+2]=D[i+i]=(int)(sqrt(d));
	D[i+i+1]=D[i+i+3]=D[i+i-1];
}


