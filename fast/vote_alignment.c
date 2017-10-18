/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"
#include "basic.h"
#include "misc.h"

int	vote_alignment(int *PAIRS,int *VOTE,int *LGC,Relation *RR1,
	Relation *RR2,int *PTD1,int *PTD2,int M,int N,int L,int *score)
{
	int	gap_a,gap_b;
	int	len;

	/* dynamic programming to get initial alignment */
	gap_a=VC_ALIGN_GAPA*L;
	gap_b=VC_ALIGN_GAPB*L;
	len=jz_dp_local_affine(VOTE,M,N,gap_a,gap_b,
		(jz_pair_int*)PAIRS,0,score);

	/* Refinement of the alignment */
	len=vote_refine(PAIRS,len,VOTE,PTD1,PTD2,RR1,RR2,M,N,score);
	len=vote_extend((jz_pair_int*)PAIRS,len,M,N);
	len=vote_trim((jz_pair_int*)PAIRS,len,4);
	len=vote_extend((jz_pair_int*)PAIRS,len,M,N);
	len=vote_shrink(PAIRS,len,M,N,RR1,RR2,NULL,NULL,score);
	len=vote_trim((jz_pair_int*)PAIRS,len,4);
	*score=vote_score(PAIRS,len,M,N,RR1,RR2);
	/*len=vote_extra(PAIRS,len,M,N,RR1,RR2,&score);*/
	return len;
}


