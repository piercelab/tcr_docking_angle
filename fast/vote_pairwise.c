/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>
#include <limits.h>

#include "vote.h"

/*-----------------------------------------------------------------------------
Align two proteins at residue level
-----------------------------------------------------------------------------*/

int	vote_pairwise(jz_protein *pt1,jz_protein *pt2,Alignment *align)
{
	int		M,N,len,score=0;
	int		*PTD1,*PTD2;
	int		*DIST1,*DIST2;
	int		*PAIRS;
	Relation	*RR1,*RR2;

	/* clear the alignment space */
	alignment_clear(align);

	/* size of the two protein structures */
	M=pt1->nres,N=pt2->nres;
	if(M<=6||N<=6)return -1;

	/* allocate space */
	JZ_ARRAY_INIT(PTD1,M*7+N*7);
	DIST1=PTD1+M*3;
	PTD2=DIST1+M+M;
	DIST2=PTD2+N*3;
	PAIRS=DIST2+N+N;
	JZ_ARRAY_INIT(RR1,M*M+N*N);
	RR2=RR1+M*M;
	
	/* pre-process protein pt1 */
	vote_discretize_residues(PTD1,pt1->residues,M);
	vote_lgc_dist(DIST1,PTD1,M);
	vote_relation(RR1,PTD1,M,0);

	/* pre-process protein pt2 */
	vote_discretize_residues(PTD2,pt2->residues,N);
	vote_lgc_dist(DIST2,PTD2,N);
	vote_relation(RR2,PTD2,N,1);

	/* call vote_align */
	len=vote_align(PTD1,PTD2,DIST1,DIST2,RR1,RR2,PAIRS,&score,
		M,N,pt1,pt2,-1);
	
	alignment_build(align,PAIRS,len);
	align->score=(float)score/1024;
	align->pvalue=align->score/sqrt(M*N+1);
	align->rmsd=vote_rmsd(PAIRS,len,pt1,pt2);

	/* free pointers and return */
	JZ_ARRAY_FREE(PTD1);
	JZ_ARRAY_FREE(RR1);

	return len;
}


