/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

void	vote_matrix(int *VOTE,Voter *CANDS,int *ENDS,
	int M,int N,int L,int infinity);

/*-----------------------------------------------------------------------------
Given a leader, find out all candidates that are alignable to the
leader. Get a new voting matrix out of the group.
-----------------------------------------------------------------------------*/

void	vote_follow(int *VOTE,Voter *CANDS,int nvoter,
	Relation *RR1,Relation *RR2,int M,int N,int leader)
{
	char		**group;
	int		i0,j0,i,j,k,K;
	Relation	*R1,*R2,*RP1,*RP2;
	float		degree;
	int		a,b,c,d;
	Voter		*NEW_CANDS;
	int		*NEW_ENDS;
	void		*chunk=NULL;
	int		n,L;

	JZ_CHUNK_INIT(chunk,1024*32);
	
	i0=CANDS[leader].i;
	j0=CANDS[leader].j;
	R1=RR1+i0*M;
	R2=RR2+j0*N;

	JZ_MATRIX_INIT(group,M,N);
	jz_charray_set(group[0],1,M*N);
	for(k=0;k<nvoter;k++)
	{
		if(CANDS[k].nvote<0)continue;
		i=CANDS[k].i;
		j=CANDS[k].j;
		group[i][j]=2;
	}

	/* test if each candidate is alignable to the leader. */
	for(k=0;k<nvoter;k++)
	{
		if(CANDS[k].nvote<0)continue;
		i=CANDS[k].i;
		j=CANDS[k].j;
		if(abs(i-i0)<VC_RR_NEAR_NEIGHBOR&&abs(j-j0)<VC_RR_NEAR_NEIGHBOR
			&&group[i][j]==2)
		{
			group[i][j]=0;
			continue;
		}

		RP1=R1+i;
		RP2=R2+j;
		if((RP1->mask[7]&RP2->mask[7])==0)continue;
		VOTE_SCORE_PAIR(d,RP1,RP2,a,b,c);
		if(d>VC_FOLLOW_THRESHOLD)group[i][j]=0;
	}

	/* seal breaches */
	for(i=0;i<M;i++)
	for(j=0;j<N;j++)
	{
		if(group[i][j]!=2)continue;
		if(i>0&&j>0&&group[i-1][j-1]==0)group[i][j]=0;
		if(i<M-1&&j<N-1&&group[i+1][j+1]==0)group[i][j]=0;
	}

	/* Get voting matrix with the group */
	JZ_ARRAY_INIT(NEW_CANDS,nvoter);
	JZ_ARRAY_INIT(NEW_ENDS,M);
	vote_candidate(NEW_CANDS,NEW_ENDS,group[0],M,N);
	K=NEW_ENDS[M-1];
	vote_tally(NEW_CANDS,NEW_ENDS,RR1,RR2,M,N,chunk);
	for(i=n=0;i<K;i++)
		if(NEW_CANDS[i].nvote>=0)
			n+=NEW_CANDS[i].nvote;
	degree=(float)n/(K*K);
	
	L=K*degree;
	k=(int)(L*VC_MATRIX_INFINITY);
	vote_matrix(VOTE,NEW_CANDS,NEW_ENDS,M,N,L,k);

	JZ_ARRAY_FREE(NEW_CANDS);
	JZ_ARRAY_FREE(NEW_ENDS);
	JZ_MATRIX_FREE(group);
	JZ_CHUNK_FREE(chunk);
}


