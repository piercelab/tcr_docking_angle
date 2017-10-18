/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
Visualize a VOTE matrix
-----------------------------------------------------------------------------*/

void	vote_matrix(int *V,Voter *VP,int *ends,
	int M,int N,int L,int infinity)
{
	int		i,k,m,d;

	jz_intarray_set(V,infinity,M*N);
	for(k=i=0;i<M-2;i++,V+=N)
	for(m=ends[i];k<m;k++)
	{
		if(VP[k].nvote<0)continue;
		d=L-VP[k].nvote;
		if(d<0)d=0;
		V[VP[k].j]=VP[k].score+(int)(d*VC_PSEUDO_NEG);
	}
}

