/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>
#include <limits.h>

#include "vote.h"

/*-----------------------------------------------------------------------------
 * List the ends of columns of each allowed cell. This is a way to store
 * sparse matrix.
-----------------------------------------------------------------------------*/

void	vote_candidate(Voter *VC,int *ends,char *A,int M,int N)
{
	Voter		*pp;
	int		i,j;

	for(i=0,pp=VC;i<M;i++,A+=N)
	{
		for(j=0;j<N;j++)
		{
			if(A[j])continue;
			pp->i=i;
			pp->j=j;
			pp->link=NULL;
			pp->nvote=pp->score=0;
			pp++;
		}	
		ends[i]=pp-VC;
	}
}

