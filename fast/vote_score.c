/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
 * Get the contribution matrix, and contribution array, return the total score
-----------------------------------------------------------------------------*/

int	vote_score(int *pairs,int len,int M,int N,
	Relation *RR1,Relation *RR2)
{
	int		i,j,i1,j1;
	int		score,s,a,b,d;
	Relation	*R1,*R2,*RRP1,*RRP2;
	int		**CTM;
	int		*CTA;
	
	if(len<2)return 0;
	JZ_MATRIX_INIT(CTM,len+1,len);
	CTA=CTM[len];
	
	for(i=0;i<len;i++)
		for(j=0;j<len;j++)
			CTM[i][j]=0;
	
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

	JZ_MATRIX_FREE(CTM);
	score/=2;
	score+=VC_SCORE_LEN*len;
	return score;
}


