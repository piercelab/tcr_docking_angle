/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote_discretize.h"

/*-----------------------------------------------------------------------------
Convert residues in a protein into integer points to speed up calculation.
By default, the discretization uses 128 level per A. This is good if
square of distance difference is calculated.

ri should be initialized as int[N*3].
-----------------------------------------------------------------------------*/

void	vote_discretize_residues(int *ri,jz_protein_residue *reses,int N)
{
	int		k;
	float		*p;
	
	for(k=0;k<N;k++,ri+=3)
	{
		p=reses[k].center;
		ri[0]=(int)(p[0]*VOTE_DISCRETIZE_COORD_FACTOR);
		ri[1]=(int)(p[1]*VOTE_DISCRETIZE_COORD_FACTOR);
		ri[2]=(int)(p[2]*VOTE_DISCRETIZE_COORD_FACTOR);
	}
}

