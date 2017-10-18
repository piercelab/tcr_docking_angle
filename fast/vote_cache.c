/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
vote_cache.c	Cache system to stored calculated values of ES(i,j,x,y)
=============================================================================*/

#include "vote.h"
#include "basic.h"

/*-----------------------------------------------------------------------------
Clean-up a cache
-----------------------------------------------------------------------------*/

void	vote_cache_clean(int **cache,int M,int N)
{
	int	k,n;

	for(n=M*N,k=0;k<n;k++)
	{
		if(!cache[k])continue;
		JZ_ARRAY_FREE(cache[k]);
		cache[k]=NULL;
	}
}

/*-----------------------------------------------------------------------------
Initial cache
-----------------------------------------------------------------------------*/

void	vote_cache_initial(int ***cache,Relation *RR1,Relation *RR2,
	int *pairs,int len,int M,int N)
{
	int		i,j;
	int		i1,j1,i2,j2;
	int		s,a,b,d;
	int		*pc;
	Relation	*R1,*R2,*RRP1,*RRP2;

	if(len<=0)return;
	for(i=0;i<len;i++)
	{
		i1=pairs[i+i];
		j1=pairs[i+i+1];
		
		R1=RR1+i1*M;
		R2=RR2+j1*N;

		JZ_ARRAY_INIT(cache[i1][j1],len);
		pc=cache[i1][j1];

		for(j=0;j<i;j++)
		{
			i2=pairs[j+j];
			j2=pairs[j+j+1];
			pc[j]=cache[i2][j2][i];
		}
		pc[i]=0;
		for(j=i+1;j<len;j++)
		{
			i2=pairs[j+j];
			j2=pairs[j+j+1];

			RRP1=R1+i2;
			RRP2=R2+j2;
			VOTE_SCORE_PAIR(s,RRP1,RRP2,a,b,d);
			pc[j]=s;
		}
	}
}

/*-----------------------------------------------------------------------------
Get mapping between a new alignment and an old one
-----------------------------------------------------------------------------*/

void	vote_map(int *map,int *pairs,int len,int *pairs_old,int len_old)
{
	int		i,j,i1,j1;
	
	for(i=j=0;i<len;i++)
	{
		i1=pairs[i+i];
		j1=pairs[i+i+1];
		for(;j<len_old;j++)
		{
			if(pairs_old[j+j]<i1)continue;
			if(pairs_old[j+j+1]<j1)continue;
			break;
		}
		if(j>=len_old)map[i]=-1;
		else if(pairs_old[j+j]==i1&&pairs_old[j+j+1]==j1)map[i]=j++;
		else map[i]=-1;
	}
}


