/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>
#include <limits.h>

#include "basic.h"
#include "misc.h"
#include "vote.h"

/*-----------------------------------------------------------------------------
 * LCA flags
-----------------------------------------------------------------------------*/

#define	VOTE_LCA_ALLOW			0	/* Pair is allowed	*/
#define	VOTE_LCA_EDGE			1	/* Pair is on an edge	*/
#define	VOTE_LCA_LGC			2	/* LGC score too low	*/
#define	VOTE_LCA_DPPLGC			4	/* DPP-LGC two low	*/
#define	VOTE_LCA_ROW			8	/* ROW discrimination	*/
#define	VOTE_LCA_COL			16	/* COL discrimination	*/
#define	VOTE_LCA_ISOLATED		32	/* Isolated pair	*/
#define	VOTE_LCA_ISON			(~32)	/* Not isolated		*/
#define	VOTE_LCA_BREAK			64	/* small break		*/
#define	VOTE_LCA_VOTE			64	/* VOTE value		*/
#define	VOTE_LCA_DPPVOTE		64	/* DPP of VOTE		*/

/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------*/

int	vote_lca_length(int N)
{
	int		len;

	/*
	len=N/3;
	if(len<30)len=30;
	if(len>=N)len=N-1;
	return len;
	*/

#define	INV_LOG2			1.442695040888963e+00	

	if(N<=0)return 0;
	len=(int)(VC_LCA_LEN_MIN*pow(VC_LCA_LEN_BASE,
		log((double)N/VC_LCA_LEN_UNIT)*INV_LOG2));
	if(len<0)len=N;
	if(len>=N)len=N-1;

	return len;
}


void	vote_trim_lgc(char *LCA,int *LGC,int M,int N)
{
	int		*buf,*DPP,*pc,*p,d;
	char		*pa;
	int		i,j,k,m,n,s,L;
	int		*P;
	int		len_row,len_col;
	int		dpp_max;
	/*int		cutoff;*/
	
	/* initialize LCA */
	memset(LCA,0,M*N);
	len_row=vote_lca_length(N);
	len_col=vote_lca_length(M);
	
	/* edges are not allowed */
	L=N;
	M-=3;
	N-=3;
	pa=LCA+M*L;
	n=L+L+L;
	for(j=0;j<n;j++)
		LCA[j]=pa[j]=VOTE_LCA_EDGE;
	LCA+=n;
	LGC+=n;

	pa=LCA-3;
	for(i=M-2;i;i--,pa+=L)
		pa[0]=pa[1]=pa[2]=pa[3]=pa[4]=pa[5]=VOTE_LCA_EDGE;
	
	/* LGC threshold */
	if(VC_LCA_LGC_MIN>=0)
	{
		pa=LCA;
		pc=LGC;
		n=((M-3)*L)&~3;
		for(j=0;j<n;j+=4)
		{
			if(pc[j]<VC_LCA_LGC_MIN)pa[j]=VOTE_LCA_LGC;
			if(pc[j+1]<VC_LCA_LGC_MIN)pa[j+1]=VOTE_LCA_LGC;
			if(pc[j+2]<VC_LCA_LGC_MIN)pa[j+2]=VOTE_LCA_LGC;
			if(pc[j+3]<VC_LCA_LGC_MIN)pa[j+3]=VOTE_LCA_LGC;
		}
	}

	/* get the DPP matrix */
	JZ_ARRAY_INIT(DPP,M*L);
	jz_dp_local_affine_potential(DPP,LGC,M,L,
		VC_LCA_DPP_GAP_OPEN,VC_LCA_DPP_GAP_EXTEND);

	/* find the max in DPP-LGC */
	dpp_max=d=jz_intarray_max(DPP,(M-3)*L);

	d=(int)((float)d*VC_LCA_DPP_CUTOFF);
	
	/* by DPPLGC */
	n=L+L;
	pa=LCA;
	P=DPP;
	for(i=3;i<M;i++,pa+=L,P+=L)
	for(j=3;j<N;j++)
	{
		s=P[j];
		if(s<d)
		{
			pa[j]=VOTE_LCA_DPPLGC;
			continue;
		}
		k=0;
		s+=VC_LCA_DPP_ADD;
		if(s>P[j-1])k++;
		if(s>P[j-2])k++;
		if(s>P[j+1])k++;
		if(s>P[j+2])k++;
		if(s>P[j-n])k++;
		if(s>P[j-L])k++;
		if(s>P[j+L])k++;
		if(s>P[j+n])k++;
		if(k<VC_LCA_DPP_NN_MIN)pa[j]=VOTE_LCA_DPPLGC;
	}
	
	/* trim row by row */
	buf=DPP;
	n=len_row;
	pc=LGC;
	pa=LCA;
	for(i=3;i<M;i++,pc+=L,pa+=L)
	{
		for(p=buf,j=3;j<N;j++)
			if(!pa[j])*p++=pc[j];
		k=p-buf;
		/*if(k<n)continue;*/
		d=jz_select_i(buf,k,k/5);
		/*
		cutoff=(int)(k*0.9);
		if(k<cutoff)continue;
		d=jz_select_i(buf,k,k-cutoff);
		*/
		for(j=3;j<N;j++)
			if(!pa[j]&&(pc[j]<d))pa[j]=VOTE_LCA_ROW;
	}

	/* trim by column */
	n=len_col; /* vote_lca_length(M); */
	for(j=3;j<N;j++)
	{
		pa=LCA+j;
		pc=LGC+j;
		for(p=buf,i=M-3;i;i--,pa+=L,pc+=L)
			if(pa[0]==0||pa[0]==VOTE_LCA_ROW)*p++=pc[0];
		k=p-buf;
		/*if(k<n)continue;*/
		d=jz_select_i(buf,k,k/5);
		/*
		cutoff=(int)(k*0.9);
		if(k<cutoff)continue;
		d=jz_select_i(buf,k,k-cutoff);
		*/
		pa=LCA+j;
		pc=LGC+j;
		for(i=M-3;i;i--,pa+=L,pc+=L)
			if(pa[0]==0&&pc[0]<d)pa[0]=VOTE_LCA_COL;
	}
	
	/* add back pairs with best DPP */
	pa=LCA;
	pc=DPP;
	n=((M-3)*L)&~3;
	for(j=0;j<n;j+=4)
	{
		if(pc[j]>=dpp_max)pa[j]=0;
		if(pc[j+1]>=dpp_max)pa[j+1]=0;
		if(pc[j+2]>=dpp_max)pa[j+2]=0;
		if(pc[j+3]>=dpp_max)pa[j+3]=0;
	}
	
	/* trim isolated */
	pa=LCA;
	n=N+4;
	m=n+n;
	for(i=3;i<M;i++,pa+=L)
	for(j=3;j<N;j++)
	{
		k=0;
		if(!pa[j])k++;
		if(!(pa[j-n]))k++;
		if(!(pa[j-m]))k++;
		if(!(pa[j-m-n]))k++;
		if(!pa[j+n])k++;
		if(!pa[j+m])k++;
		if(!pa[j+m+n])k++;
		if(k<VC_LCA_ISOLATE_LIMIT)pa[j]=VOTE_LCA_ISOLATED;
		else if(pa[j])pa[j]=VOTE_LCA_BREAK;
	}
	
	/* seal small breaks */
	pa=LCA;
	for(i=3;i<M;i++,pa+=L)
	for(j=3;j<N;j++)
		if(pa[j]==VOTE_LCA_BREAK)pa[j]=0;
	
	/* free pointers */
	JZ_ARRAY_FREE(DPP);
}


