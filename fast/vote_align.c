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
vote_align pairwise aligns two pre-processed protein structures. This is
the kernel of the VOTE algorithm. vote_pairwise() simply pre-process those
protein structures and call vote_align(). In vote_search(), database
proteins are already pre-processed. 
-----------------------------------------------------------------------------*/

int	vote_align(int *PTD1,int *PTD2,int *DIST1,int *DIST2,
	Relation *RR1,Relation *RR2,int *PAIRS,int *pscore,
	int M,int N,jz_protein *pt1,jz_protein *pt2,float abort_factor)
{
	int		*LGC;
	int		*VOTE;
	char		*ALLOW;
	Voter		*CANDS;
	int		nvoter,K,L;
	int		*ENDS;
	int		len=0,score=0;
	float		degree,last_degree=0;
	int		*TMP_PAIRS;
	int		tmp_len;
	int		tmp_score;
	int		round=0;
	int		i,k,m,n;
	double		factor;
	void		*chunk=NULL;

	/* calculate local geometric compatibility matrix */
	JZ_ARRAY_INIT(LGC,M*N);
	JZ_ARRAY_INIT(VOTE,M*N);
	vote_lgc(LGC,DIST1,DIST2,M,N);
	
	/* trim LGC to get initial ALLOW */
	JZ_ARRAY_INIT(ALLOW,M*N);
	vote_trim_lgc(ALLOW,LGC,M,N);
	
	/* get candidate list */
	JZ_ARRAY_INIT(ENDS,M);
	nvoter=jz_charray_count_zero(ALLOW,M*N);
	JZ_ARRAY_INIT(CANDS,nvoter);
	vote_candidate(CANDS,ENDS,ALLOW,M,N);
	JZ_CHUNK_INIT(chunk,1024*256);

	if(nvoter<3000)
		vote_tally(CANDS,ENDS,RR1,RR2,M,N,chunk);
	else
	{
		factor=0.145*(2-exp(-(M*N)/20000));	
		vote_tally_fast(CANDS,nvoter,RR1,RR2,chunk,M,N,64,2,factor);
	}	

	/* Iteration to get/update the voting matrix */
	for(;;)
	{
		/* count the number of remaining voters */
		for(i=K=n=0;i<nvoter;i++)
			if(CANDS[i].nvote>=0)K++,n+=CANDS[i].nvote;
		
		/* calculate degree of unamity */
		degree=(float)n/(K*K);

		/* pseudo length is the estimated alignment length */
		L=(int)(K*degree);

		/* abort if the alignment is hopeless */
		if(abort_factor>0)
		{
			if(L*L*L*L<abort_factor*M*N)
			{
				fprintf(stderr,"aborted.\n");
				goto __EXIT;
			}	
		}

		/* convert CANDS into VOTE matrix */
		vote_matrix(VOTE,CANDS,ENDS,M,N,L,
			(int)(L*VC_MATRIX_INFINITY));
		
		/*fprintf(stderr,"NV=%d, K=%d, degree=%.3e, n=%d\n",
			nvoter,K,degree,n);*/
		
		/* stop if degree acceptable or no progress */	
		if(degree>=VC_DEGREE_GREAT)break;
		if((degree-last_degree)<=(degree*VC_ELIM_DEGREE_LIMIT))break;
		if(K<(M+N))break;
		last_degree=degree;

		/* drop low scoring pairs */
		vote_eliminate(VOTE,CANDS,ENDS,K,degree,M,N);
		vote_consolidate(CANDS,nvoter);
	}


	/* the above iteration converged to global optimum */
	if(degree>=VC_DEGREE_GOOD)
	{
		len=vote_alignment(PAIRS,VOTE,LGC,RR1,RR2,PTD1,PTD2,
			M,N,L,&score);
		round++;	
	}

	/* possibly tangled with local minima */
	{
		JZ_ARRAY_INIT(TMP_PAIRS,2*(M+N));
		for(;;)
		{
			/* choose a leader by number of vote received */
			for(k=m=i=0;i<nvoter;i++)
			{
				if(CANDS[i].nvote<m)continue;
				m=CANDS[k=i].nvote;
			}

			/* choose followers and get an alignment */
			vote_follow(VOTE,CANDS,nvoter,RR1,RR2,M,N,k);

			tmp_len=vote_alignment(TMP_PAIRS,VOTE,LGC,RR1,RR2,
				PTD1,PTD2,M,N,L,&tmp_score);
			
			/* update the best alignment */
			if(tmp_score>score)
			{
				JZ_ARRAY_COPYF(PAIRS,TMP_PAIRS,tmp_len*2);
				len=tmp_len;
				score=tmp_score;
			}

			++round;
			if(round>=VC_ALIGN_ROUND)break;
	
			/* remove the alignment from the committee */
			vote_remove(CANDS,nvoter,TMP_PAIRS,tmp_len,M,N);
			vote_consolidate(CANDS,nvoter);
		}

		JZ_ARRAY_FREE(TMP_PAIRS);
	}

__EXIT:

	if(chunk)JZ_CHUNK_FREE(chunk);
	JZ_ARRAY_FREE(VOTE);
	JZ_ARRAY_FREE(ALLOW);
	JZ_ARRAY_FREE(ENDS);
	JZ_ARRAY_FREE(CANDS);
	*pscore=score;

	return len;
}




