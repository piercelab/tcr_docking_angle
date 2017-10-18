/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
Vote in an band.
-----------------------------------------------------------------------------*/

#define	SHIFT		VC_VOTE_SHIFT

void	vote_tally_band(Voter **CANDS,int *ENDS,Relation *RR1,Relation *RR2,
	void *chunk,int M,int N,int di_min,int di_max,int dj_min,int dj_max)
{
	int		i,j,i1,j1,k,k1,i1_min,i1_max,m,n;
	int		a,b,d,bit,score;
	Relation	*R1,*R2,*RP1,*RP2;
	VoteEdge	*node;
	
	if(di_min<VC_VOTE_SHIFT)di_min=VC_VOTE_SHIFT;
	if(dj_min<VC_VOTE_SHIFT)dj_min=VC_VOTE_SHIFT;
	
	for(R1=RR1,k=i=0;i<M;i++,R1+=M)
	for(m=ENDS[i];k<m;k++)
	{
		j=CANDS[k]->j;
		R2=RR2+j*N;
		i1_min=i+di_min;
		i1_max=i+di_max;
		if(i1_max>M-1)i1_max=M-1;
		for(i1=i+SHIFT;i1<=i1_max;i1++)
		{
			k1=ENDS[i1-1];
			RP1=R1+i1;
			if(!RP1->mask[0])continue;
			
			n=ENDS[i1];

			if(i1<i1_min)j1=j+dj_min;
			else j1=j+SHIFT;
			
			for(;k1<n&&CANDS[k1]->j<j1;k1++);

			for(j1=j+dj_max;k1<n;k1++)
			{
				if(CANDS[k1]->j>j1)break;
				RP2=R2+CANDS[k1]->j;
				VOTE_SCORE(score,RP1,RP2,a,b,d,bit,continue);

				CANDS[k]->score+=score;
				CANDS[k]->nvote++;
				node=(VoteEdge*)JZ_CHUNK_ALLOC
					(chunk,sizeof(VoteEdge));
				node->source=CANDS[k1];
				node->score=score;
				node->next=CANDS[k]->link;
				CANDS[k]->link=node;

				CANDS[k1]->score+=score;
				CANDS[k1]->nvote++;
				node=(VoteEdge*)JZ_CHUNK_ALLOC
					(chunk,sizeof(VoteEdge));
				node->source=CANDS[k];
				node->score=score;
				node->next=CANDS[k1]->link;
				CANDS[k1]->link=node;
			}
		}
	}
}

/*-----------------------------------------------------------------------------
Stepwise voting
-----------------------------------------------------------------------------*/


void	vote_tally_fast(Voter *Voters,int nvoter,Relation *RR1,Relation *RR2,
	void *chunk,int M,int N,int L,float base,float fraction)
{
	int		i,j,k,m,n,level;
	double		alph=0,beta=0;
	double		di_min,di_max,dj_min,dj_max;
	Voter		**CANDS;
	int		*ENDS;
	int		*buf;
	int		*VOTE,*pv;
	int		nv;

	/* find the level of the voting graph */
	level=(int)(log((double)(M*N)/(L*L))/log(base)+0.5);
	if(level<0)level=0;

	/* calculate factor for sizes M and N */
	if(level>0)
	{
		alph=exp(log((double)M/L)/level);
		beta=exp(log((double)N/L)/level);
	}

	/* construct voter links */
	JZ_ARRAY_INIT(CANDS,nvoter);
	JZ_ARRAY_INIT(ENDS,M);
	JZ_ARRAY_INIT(buf,nvoter);
	for(i=0;i<nvoter;i++)
	{
		CANDS[i]=Voters+i;
		CANDS[i]->score=CANDS[i]->nvote=0;
	}
	
	for(i=0,j=0;j<nvoter;j++)
	{
		if(Voters[j].i<=i)continue;
		for(k=Voters[j].i;i<k;i++)ENDS[i]=j;
	}
	for(;i<M;i++)ENDS[i]=nvoter;

	di_min=0;
	di_max=L;
	dj_min=0;
	dj_max=L;
	n=nvoter;

	JZ_ARRAY_INIT(VOTE,M*N);
	jz_intarray_set(VOTE,INT_MIN,M*N);

	/* vote and screen for each level k */
	for(k=0;k<=level;k++)
	{
/*fprintf(stderr,"di: %d-%d, dj: %d-%d, n=%d\n",
	(int)di_min,(int)di_max,(int)dj_min,(int)dj_max,n);*/
	
		/* cast vote in band */
		vote_tally_band(CANDS,ENDS,RR1,RR2,chunk,M,N,
			(int)di_min,(int)di_max,(int)dj_min,(int)dj_max);

		if(k==level)break;
		
		/* map score to VOTE matrix */
		for(i=0;i<n;i++)
			VOTE[CANDS[i]->i*N+CANDS[i]->j]=CANDS[i]->nvote;
		
		/* eliminate a fraction of voters */
		for(i=0;i<n;i++)
			buf[i]=CANDS[i]->nvote;
		if(k==0)m=jz_select_i(buf,n,(int)(n*fraction*0.7));
		else m=jz_select_i(buf,n,(int)(n*fraction));
		if(m<2)m=2;

		for(i=0;i<n;i++)
		{
			nv=CANDS[i]->nvote;
			pv=VOTE+CANDS[i]->i*N+CANDS[i]->j;
			if(nv<m)
			{
				CANDS[i]->nvote=INT_MIN;
				pv[0]=INT_MIN;
			}
			nv=nv*VC_ELIM_NEIGHBOR_FACTOR;
			if(nv<pv[-1]||nv<pv[1]||nv<pv[-N]||nv<pv[N])
			{
				CANDS[i]->nvote=INT_MIN;
				pv[0]=INT_MIN;
			}	
		}		
		
		/* consolidate the change */
		vote_consolidate2(CANDS,n);

		for(i=j=0;j<n;j++)
		{
			if(CANDS[j]->nvote<0)continue;
			CANDS[i++]=CANDS[j];
		}
		for(n=i,i=0,j=0;j<n;j++)
		{
			if(CANDS[j]->i<=i)continue;
			for(m=CANDS[j]->i;i<m;i++)
				ENDS[i]=j;
		}
		for(;i<M;i++)ENDS[i]=n;

		/* update range for the next round of vote */
		di_min=di_max+1;
		dj_min=dj_max+1;
		di_max*=alph;
		dj_max*=beta;
	}

	JZ_ARRAY_FREE(CANDS);
	JZ_ARRAY_FREE(ENDS);
	JZ_ARRAY_FREE(buf);
	JZ_ARRAY_FREE(VOTE);
}


