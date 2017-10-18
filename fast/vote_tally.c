/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>
#include <limits.h>

#include "vote.h"

int	__count[128];

/*-----------------------------------------------------------------------------
 * Given allow matrix, get vote matrix, also return uncertainty
-----------------------------------------------------------------------------*/

#define	SHIFT		VC_VOTE_SHIFT

void	vote_tally(Voter *CANDS,int *ENDS,Relation *RR1,
	Relation *RR2,int M,int N,void *chunk)
{
	int		i,j,j1,k,i1,k1,m,n;
	int		a,b,d,score;
	Relation	*R1,*R2,*RP1,*RP2;
	VoteEdge	*node;
	int		bit;
	int		dj_max;
	int		dj_min;
	int		dj;
	
	n=ENDS[M-1];
	for(k=0;k<n;k++)
		CANDS[k].score=CANDS[k].nvote=0;

	/* the comparison between CANDS[k] and CANDS[k1] */
	for(R1=RR1,k=i=0;i<M-SHIFT;i++,R1+=M)
	for(m=ENDS[i];k<m;k++)
	{
		j=CANDS[k].j;
		R2=RR2+j*N;
		for(i1=i+SHIFT;i1<M;i1++)
		{
			k1=ENDS[i1-1];
			RP1=R1+i1;
			if(!RP1->mask[0])continue;
			
			n=ENDS[i1];
			for(j1=j+SHIFT;k1<n&&CANDS[k1].j<j1;k1++);

			dj_max=i1-i;
			dj_min=dj_max/8-10;
			dj_max=dj_max*8+10;
			if(dj_max>1024)dj_max=1024;

			for(;k1<n;k1++)
			{
				dj=CANDS[k1].j-j;
				if(dj>dj_max)continue;
				if(dj<dj_min)continue;
				RP2=R2+CANDS[k1].j;

				/* score calculation */
				VOTE_SCORE(score,RP1,RP2,a,b,d,bit,continue);

				/* credit the score */
				CANDS[k].score+=score;
				CANDS[k].nvote++;
				node=(VoteEdge*)JZ_CHUNK_ALLOC
					(chunk,sizeof(VoteEdge));
				node->source=CANDS+k1;
				node->score=score;
				node->next=CANDS[k].link;
				CANDS[k].link=node;
				
				CANDS[k1].score+=score;
				CANDS[k1].nvote++;
				node=(VoteEdge*)JZ_CHUNK_ALLOC
					(chunk,sizeof(VoteEdge));
				node->source=CANDS+k;
				node->score=score;
				node->next=CANDS[k1].link;
				CANDS[k1].link=node;
			}
		}
	}

	/*fprintf(stderr,"N=%d\n",n);*/
}



