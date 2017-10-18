/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "basic.h"
#include "vote.h"

/*-----------------------------------------------------------------------------
 * Get the contribution matrix, and contribution array, return the total score
-----------------------------------------------------------------------------*/

int	vote_CONTR_old(int **CTM,int *CTA,int *pairs,int len,int M,int N,
	int CTLD,Relation *RR1,Relation *RR2)
{
	int		i,j,i1,j1;
	int		score,s,a,b,d;
	Relation	*R1,*R2,*RRP1,*RRP2;
	
	for(i=0;i<len;i++)CTM[i][i]=0;
	
	for(i=0,--len;i<len;i++)
	{
		i1=pairs[i+i];
		j1=pairs[i+i+1];
		R1=RR1+i1*M;
		R2=RR2+j1*N;
		for(j=i+1;j<=len;j++)
		{
			if(pairs[j+j]==i1)continue;
			if(pairs[j+j+1]==j1)continue;
			RRP1=R1+pairs[j+j];
			RRP2=R2+pairs[j+j+1];
			VOTE_SCORE_PAIR(s,RRP1,RRP2,a,b,d);
			CTM[i][j]=CTM[j][i]=s;
		}
	}

	for(score=i=0;i<=len;i++)
	{
		for(s=j=0;j<=len;j++)s+=CTM[i][j];
		CTA[i]=s;
		score+=s;
	}

	return score/2;
}

/*-----------------------------------------------------------------------------
 * Get the contribution matrix, and contribution array, return the total score
-----------------------------------------------------------------------------*/

int	vote_CONTR(int **CTM,int *CTA,int *pairs,int len,int M,int N,
	int CTLD,Relation *RR1,Relation *RR2,int ***cache,int *map)
{
	int		i,j,i1,j1;
	int		score,s,a,b,d;
	int		*pc;
	Relation	*R1,*R2,*RRP1,*RRP2;
	
	for(i=0;i<len;i++)CTM[i][i]=0;
	
	for(i=0,--len;i<len;i++)
	{
		i1=pairs[i+i];
		j1=pairs[i+i+1];
		R1=RR1+i1*M;
		R2=RR2+j1*N;
		pc=cache[i1][j1];
		if(pc)
		{
			for(j=i+1;j<=len;j++)
			{
				if(map[j]>=0)s=pc[map[j]];
				else
				{
					RRP1=R1+pairs[j+j];
					RRP2=R2+pairs[j+j+1];
					VOTE_SCORE_PAIR(s,RRP1,RRP2,a,b,d);
				}
				CTM[i][j]=CTM[j][i]=s;
			}
		}
		else
		{
			for(j=i+1;j<=len;j++)
			{
				if(pairs[j+j]==i1)continue;
				if(pairs[j+j+1]==j1)continue;
				RRP1=R1+pairs[j+j];
				RRP2=R2+pairs[j+j+1];
				VOTE_SCORE_PAIR(s,RRP1,RRP2,a,b,d);
				CTM[i][j]=CTM[j][i]=s;
			}		
		}
	}

	for(score=i=0;i<=len;i++)
	{
		for(s=j=0;j<=len;j++)s+=CTM[i][j];
		CTA[i]=s;
		score+=s;
	}

	return score/2;
}

/*-----------------------------------------------------------------------------
 * Trim gaps in an alignment
-----------------------------------------------------------------------------*/

int	vote_trim_gaps(int *pairs,int len)
{
	int		i,j;

	if(!pairs)return -1;
	len+=len;
	for(i=j=0;j<len;j+=2)
	{
		if(pairs[j]<0||pairs[j+1]<0)continue;
		pairs[i]=pairs[j];
		pairs[i+1]=pairs[j+1];
		i+=2;
	}
	return i/2;
}

/*-----------------------------------------------------------------------------
 * Shrink initial alignment to get alignment core. Pairs that have low
 * contributions to the total score are dropped

 Count is the number of positive votes an pair received.
-----------------------------------------------------------------------------*/

int	vote_shrink(int *pairs,int len,int M,int N,
	Relation *RR1,Relation *RR2,int ***cache,int *map,int *pscore)
{
	int		**CTM,*CTA,*PM;
	int		i,k,min_id;
	int		contr_min;
	int		score;
	
	if(len<=4)
	{
		*pscore=0;
		return 0;
	}
	
	JZ_MATRIX_INIT(CTM,len+1,len);
	CTA=CTM[len];
	
	/* get contribution matrix and array */
	if(cache&&map)score=vote_CONTR(CTM,CTA,pairs,len,M,N,len,
		RR1,RR2,cache,map);
	else score=vote_CONTR_old(CTM,CTA,pairs,len,M,N,len,RR1,RR2);

	for(i=len;i>2;i--)
	{
		/* find the pair that contributes least */
		contr_min=INT_MAX;
		min_id=-1;
		for(k=0;k<len;k++)
		{
			if(pairs[k+k]<0)continue;
			if(CTA[k]>=contr_min)continue;
			min_id=k;
			contr_min=CTA[k];
		}
		if(min_id<0)break;
		
		/* Alignment core is defined as following:
		 * 1. Each pair has positive contribution
		 * 2. Contr_min and contr_max are close enough */
		if(contr_min>=0)break;
		
		/* CORE condition is not met, drop the lowest scoring pair */
		pairs[min_id+min_id]=-1;
		score-=contr_min;
		PM=CTM[min_id];
		for(k=0;k<len;k++)CTA[k]-=PM[k];
	}

	len=vote_trim_gaps(pairs,len);
	JZ_MATRIX_FREE(CTM);
	*pscore=score;
	return len;
}


