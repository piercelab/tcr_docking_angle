/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>
#include <limits.h>

#include "vote.h"

/*-----------------------------------------------------------------------------
 * Compare two proteins by local geometry. Return a score matrix.
 *
 * For each residue at position k, residues at positions k-2 to k+2 are
 * considered. The elastic score using 7 residues is assesed. This
 * version does not utilize SSE info at all.
 *
 * 	*----*----*----+----*----*----*
 * 
-----------------------------------------------------------------------------*/

void	vote_lgc(int *LGC,int *D1,int *D2,int M,int N)
{
	int		a,d,X;
	int		i,j,n,K;

	/* LGC calculation */
	n=N+1;
	K=(N-1);
	for(j=0;j<N;j++)LGC[j]=0;
	for(i=0;i<M-1;i++,LGC+=N)
	{
		X=D1[i+i];
		LGC[0]=0;
		for(j=0;j<K;j++)
		{
			d=D2[j+j];
			a=d+X+VC_LGC_RELAX;
			d-=X;
			if(d<0)d=-d;
			d=VC_LGC_W3*d/a;
			LGC[j]+=d;
			LGC[j+n]=d;
		}
		d=D2[j+j];
		a=d+X+VC_LGC_RELAX;
		d-=X;
		d=VC_LGC_W3*d/a*d/a;
		LGC[j]+=d;

		X=D1[i+i+1];
		for(j=0;j<N;j++)
		{
			d=D2[j+j+1];
			a=d+X+VC_LGC_RELAX;
			d-=X;
			if(d<0)d=-d;
			d=VC_LGC_W4*d/a;
			if(d<LGC[j])d=LGC[j];
			LGC[j]=VC_LGC_MAX-d;
			/*LGC[j]=VC_LGC_MAX-LGC[j]-d;*/
		}
	}

	/* i==M-1 */
	X=D1[i+i];
	LGC[0]=0;
	for(j=0;j<N-1;j++)
	{
		d=D2[j+j];
		a=d+X+VC_LGC_RELAX;
		d-=X;
		if(d<0)d=-d;
		d=VC_LGC_W3*d/a;
		LGC[j]+=d;
	}
	d=D2[j+j];
	a=d+X+256;
	d-=X;
	if(d<0)d=-d;
	d=VC_LGC_W3*d/a;
	LGC[j]+=d;

	X=D1[i+i+1];
	for(j=0;j<N;j++)
	{
		d=D2[j+j+1];
		a=d+X+VC_LGC_RELAX;
		d-=X;
		if(d<0)d=-d;
		d=VC_LGC_W4*d/a;
		if(d<LGC[j])d=LGC[j];
		LGC[j]=VC_LGC_MAX-d;
		/*LGC[j]=VC_LGC_MAX-LGC[j]-d;*/
	}

	LGC-=(M-1)*N;
	/*for(i=M*N-1;i;i--)LGC[i]+=32;*/
	/*vote_add(LGC,M,N);*/
}


