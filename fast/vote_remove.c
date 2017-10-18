/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
Remove an alignment
-----------------------------------------------------------------------------*/

void	vote_remove(Voter *VC,int nvoter,int *PAIRS,int len,int M,int N)
{
	int		i,j;
	int		i1,j1,i2,j2;
	
	for(i=j=0;i<nvoter&&j<len;)
	{
		if(VC[i].nvote<0)
		{
			i++;
			continue;
		}

		i1=VC[i].i;
		j1=VC[i].j;

		i2=PAIRS[j+j];
		j2=PAIRS[j+j+1];

		if(i1<i2)i++;
		else if(i1>i2)j++;
		else if(j1<j2)i++;
		else if(j1>j2)j++;
		else
		{
			VC[i].nvote=INT_MIN;
			i++;
			j++;
		}
	}
}


