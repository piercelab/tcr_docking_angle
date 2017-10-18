/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
Drop candidates with low DPP score or VOTE score. Gap penalty is
determined by the pseudo length of the alignment.
-----------------------------------------------------------------------------*/

void	vote_eliminate(int *VOTE,Voter *VC,int *ends,
	int K,int L,int M,int N)
{
	int		*DPP,*buf;
	int		gap_a,gap_b;
	int		i,j,k,m,n,d;
	int		*P,*pi;
	int		*pv;
	int		score;
	int		nvote_cutoff;
	
	/* Get DPP of VOTE */
	gap_a=VC_ELIM_GAPA*L;
	gap_b=VC_ELIM_GAPB*L;

	/* get DPP matrix */
	JZ_ARRAY_INIT(DPP,M*N);
	jz_dp_local_affine_potential(DPP,VOTE,M,N,gap_a,gap_b);
	
	/* save DPP score */
	JZ_ARRAY_INIT(buf,K+K);
	for(k=i=m=0,P=DPP,pi=buf;i<M;i++,P+=N)
	for(n=ends[i];k<n;k++)
	{
		if(VC[k].nvote<0)continue;
		d=P[VC[k].j];
		if(d>m)m=d;
		VC[k].dpp_score=d;
		pi[K]=VC[k].nvote;
		*pi++=d;
	}
	JZ_ARRAY_FREE(DPP);
	
	/* determine thresholds */
	d=jz_select_i(buf,K,(int)(K*VC_ELIM_DPP_FRACTION));
	m=(int)(m*VC_ELIM_DPP_LIMIT);
	if(d>m)d=m;
	nvote_cutoff=jz_select_i(buf+K,K,(int)(K*VC_ELIM_NVOTE_FRACTION));
	m=(int)(L*VC_ELIM_NVOTE_CUTOFF);
	if(nvote_cutoff<m)nvote_cutoff=m;
	JZ_ARRAY_FREE(buf);

	/* eliminate bad candidates */
	for(k=ends[0],i=1,pv=VOTE+N;i<M-1;i++,pv+=N)
	for(n=ends[i];k<n;k++)
	{
		/*fprintf(stderr,"k=%d\n",k);*/
		if(VC[k].nvote<0)continue;
		j=VC[k].j;
		if(VC[k].nvote<nvote_cutoff)
		{
			VC[k].nvote=INT_MIN;
			continue;
		}
		if(VC[k].dpp_score<d)
		{
			VC[k].nvote=INT_MIN;
			continue;
		}
		score=VC[k].score;
		score=pv[j]*VC_ELIM_NEIGHBOR_FACTOR;		
		if(score<pv[j-1]||score<pv[j+1]||score<pv[j-N]||score<pv[j+N])
			VC[k].nvote=INT_MIN;
	}
}


