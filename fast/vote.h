/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
VOTE_H		The VOTE protein structure alignment algorithm
=============================================================================*/

#ifndef	__VOTE_H
#define	__VOTE_H

#include <limits.h>
#include <math.h>

#include "basic.h"
#include "misc.h"
#include "jz_protein.h"
#include "Alignment.h"
#include "vote_discretize.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*=============================================================================
For best performance, computation in VOTE uses fixed point arithmetics.
Distances/coordinates are discretized as:
	D=(int)d*128.0, or 1A = 128 unit
Angles are disctetized as:
	A=a*3.259493234522017e+02, or 1024 = PI

In a unit vector, the max value in each dimension is 32768, or,
the sum of squares equals to 1073741824. To calculate the angle between
two unit vectors, first calculate the dot of the two vector, right shift
19 bits, then look up table vote_acos_lut[].
=============================================================================*/

#define	VOTE_DISCRETIZE_COORD_FACTOR	128.0			/* 128  = 1A */
#define	VOTE_DISCRETIZE_ANGLE_FACTOR	3.259493234522017e+02	/* 1024 = PI */
#define	VOTE_DISCRETIZE_UNIT_VECTOR	32768			/* 32768= 1  */
#define	VOTE_ACOS_SHIFT			19
#define	VOTE_ACOS_LUT_SHIFT		2052

/*-----------------------------------------------------------------------------
Make a vector a unit vector.
-----------------------------------------------------------------------------*/

#ifndef VOTE_UNIT_VECTOR
#define	VOTE_UNIT_VECTOR(V,d)						\
do{\
	JZ_VEC3_DOT((d),(V),(V));\
	(d)=(int32_t)(sqrt(d+1));\
	JZ_VEC3_MUL((V),(V),VOTE_DISCRETIZE_UNIT_VECTOR);\
	JZ_VEC3_DIV((V),(V),(d));\
}while(0)
#endif

/*-----------------------------------------------------------------------------
Convert residues in a protein into integer points to speed up calculation.
By default, the discretization uses 128 level per A. This is good if
square of distance difference is calculated.
-----------------------------------------------------------------------------*/

void	vote_discretize_residues(int *ri,jz_protein_residue *reses,int N);

/*=============================================================================
The modified elastic score with directionality terms. An integer score
will be placed in variable vs. Make sure R1, R2 are Relation* and a,b,d
are int variables.
=============================================================================*/

#define	VC_SCORE_DISTANCE_FACTOR	10240
#define	VC_SCORE_DISTANCE_FACTOR2	10240
#define	VC_SCORE_DIST_BUF		2560
#define	VC_SCORE_DIST_BUF2		2750
#define	VC_SCORE_ANGLE_BUF		2/3
#define	VC_SCORE_MAX			1024
#define	VC_SCORE_ANGLE_FACTOR		4
#define	VC_SCORE_ANGLE_FACTOR2		7/2
#define	VC_SCORE_LUTIDX_SHIFT		4
#define	VC_SCORE_SHIFT			13

/*-----------------------------------------------------------------------------
Score calculation for the vote_tally() routine. Negative votes are omitted.
-----------------------------------------------------------------------------*/

#define	VOTE_SCORE(vs,R1,R2,a,b,d,bit,skip)				\
{\
	bit=JZ_CASTR(R1->mask[0],int32_t)&JZ_CASTR(R2->mask[0],int32_t);\
	if((bit&0xFFFF)==0||(bit&0xFFFF0000)==0)skip;\
	bit=JZ_CASTR(R1->mask[2],int32_t)&JZ_CASTR(R2->mask[2],int32_t);\
	if((bit&0xFFFF)==0||(bit&0xFFFF0000)==0)skip;\
	bit=JZ_CASTR(R1->mask[4],int32_t)&JZ_CASTR(R2->mask[4],int32_t);\
	if((bit&0xFFFF)==0||(bit&0xFFFF0000)==0)skip;\
	if((R1->mask[6]&R2->mask[6])==0)skip;\
	d=R1->distance+R2->distance;\
	a=R1->distance-R2->distance;\
	b=d+VC_SCORE_DIST_BUF;\
	vs=VC_SCORE_DISTANCE_FACTOR*abs(a)/b;\
	if(vs>=VC_SCORE_MAX)skip;\
	a=R1->angles[4]-R2->angles[4];\
	if(a<0)a=-a;\
	b=R1->angles[5]-R2->angles[5];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	a=a*VC_SCORE_ANGLE_BUF;\
	b=R1->angles[0]-R2->angles[0];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	b=R1->angles[1]-R2->angles[1];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	b=R1->angles[2]-R2->angles[2];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	b=R1->angles[3]-R2->angles[3];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	a=a*VC_SCORE_ANGLE_FACTOR;\
	if(a>vs)vs=a;\
	vs=VC_SCORE_MAX-vs;\
	if(vs<=0)skip;\
	vs*=(vote_decay_lut[d>>VC_SCORE_LUTIDX_SHIFT]);\
	vs>>=VC_SCORE_SHIFT;\
}

/*-----------------------------------------------------------------------------
Score calculation for any given pair of pairs. 
-----------------------------------------------------------------------------*/

#define	VOTE_SCORE_PAIR(vs,R1,R2,a,b,d)				\
{\
	if(R1->mask[7]!=R2->mask[7])\
	{\
		vs=-1;\
	}\
	else if((R1->mask[0]==0)||(R2->mask[0]==0))\
	{\
		vs=-1;\
	}\
	else{\
	d=R1->distance+R2->distance;\
	a=R1->distance-R2->distance;\
	b=d+VC_SCORE_DIST_BUF2;\
	vs=VC_SCORE_DISTANCE_FACTOR2*abs(a)/b;\
	a=R1->angles[4]-R2->angles[4];\
	if(a<0)a=-a;\
	b=R1->angles[5]-R2->angles[5];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	a=a*VC_SCORE_ANGLE_BUF;\
	b=R1->angles[0]-R2->angles[0];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	b=R1->angles[1]-R2->angles[1];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	b=R1->angles[2]-R2->angles[2];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	b=R1->angles[3]-R2->angles[3];\
	if(b<0)b=-b;\
	if(b>a)a=b;\
	a=a*VC_SCORE_ANGLE_FACTOR2;\
	if(a>vs)vs=a;\
	vs=VC_SCORE_MAX-vs;\
	vs*=(vote_decay_lut[d>>VC_SCORE_LUTIDX_SHIFT]);\
	vs>>=VC_SCORE_SHIFT;\
	}\
}

/*-----------------------------------------------------------------------------
LUT's to calculate ES
-----------------------------------------------------------------------------*/

extern	int32_t		vote_decay_lut[];
extern	int16_t		vote_acos_lut[];
extern	uint16_t	vote_distance_bit_lut[];
extern	uint16_t	vote_distance_mask_lut[];
extern	uint16_t	vote_angle_bit_lut[];
extern	uint16_t	vote_angle_mask_lut[];
extern	uint16_t	vote_angle_mask_lut2[];


/*=============================================================================
Voting Constants (built-in parameters) governing the alignment algorithm.
LGC	Local Geometric Compatibility
LCA	Local Compatibility Allow
RR	Relationship
=============================================================================*/

#define	VC_LGC_W3			10000
#define	VC_LGC_W4			20000
#define	VC_LGC_MAX			2048
#define	VC_LGC_RELAX			128

#define	VC_LCA_LGC_MIN			0
#define	VC_LCA_DPP_GAP_OPEN		-16
#define	VC_LCA_DPP_GAP_EXTEND		-16
#define	VC_LCA_DPP_CUTOFF		0.33
#define	VC_LCA_DPP_ADD			800
#define	VC_LCA_DPP_NN_MIN		6
#define	VC_LCA_LEN_MIN			24
#define	VC_LCA_LEN_BASE			1.82
#define	VC_LCA_LEN_UNIT			64
#define	VC_LCA_ISOLATE_LIMIT		4

#define	VC_RR_NEAR_NEIGHBOR		1
#define	VC_RR_DISTANCE_MAX		8192
#define	VC_RR_ORDER_LT			1
#define	VC_RR_ORDER_GT			2
#define	VC_VOTE_SHIFT			2
#define	VC_PSEUDO_NEG			-8

#define	VC_ELIM_GAPA			-64
#define	VC_ELIM_GAPB			-8
#define	VC_ELIM_DPP_FRACTION		0.09
#define	VC_ELIM_NVOTE_FRACTION		0.12
#define	VC_ELIM_DPP_LIMIT		0.95
#define	VC_ELIM_DEGREE_LIMIT		0.03
#define	VC_DEGREE_GOOD			0.2
#define	VC_DEGREE_GREAT			0.4
#define	VC_ALIGN_ROUND			3
#define	VC_ELIM_NVOTE_CUTOFF		0.3
#define	VC_ELIM_NEIGHBOR_FACTOR		3
#define	VC_MATRIX_INFINITY		(-128)
#define	VC_ALIGN_GAPA			-64
#define	VC_ALIGN_GAPB			-8
#define	VC_FOLLOW_THRESHOLD		-128

#define	VC_REFINE_ALLOW_MAX(d)		((d)*11/9+420)
#define	VC_REFINE_ALLOW_MIN(d)		((d)*9/11-336)
#define	VC_REFINE_SEQ_TH(len)		((len)/4)
#define	VC_REFINE_BAD_SCORE		-128
#define	VC_REFINE_GAPA			-64
#define	VC_REFINE_GAPB			-8
#define	VC_SCORE_LEN			512

#define	VC_REFINE_ROUND			6


/*=============================================================================
Structured Types
=============================================================================*/

/*-----------------------------------------------------------------------------
Relation is the intra-molecular relationship between a pair of residues
within the same protein structure. We have distance and six angles.
Mask is used to quickly compare a pair of relations. 32 bytes in size.
-----------------------------------------------------------------------------*/

typedef	struct		Relation
{
	uint16_t	mask[8];
	uint16_t	angles[6];
	int32_t		distance;
}Relation;

/*-----------------------------------------------------------------------------
Voter is an vertex in the product graph in which each vertex represents an
possible pairing between residue i in protein A and residue j in protein B. 
-----------------------------------------------------------------------------*/

typedef	struct		Voter
{
	int32_t		j;
	int32_t		i;
	int32_t		score;
	int32_t		nvote;
	int32_t		dpp_score;
	struct		VoteEdge *link;
}Voter;

/*-----------------------------------------------------------------------------
VoteEdge is an positive vote between two vertices. All edges associated with
a vertex are organized into a linked list.
-----------------------------------------------------------------------------*/

typedef	struct		VoteEdge
{
	struct		VoteEdge	*next;
	struct		Voter		*source;
	int32_t		score;
}VoteEdge;

/*=============================================================================
Function prototypes
=============================================================================*/

int	vote_pairwise(jz_protein *pt1,jz_protein *pt2,Alignment *align);

int	vote_align(int *PTD1,int *PTD2,int *DIST1,int *DIST2,
	Relation *RR1,Relation *RR2,int *PAIRS,int *pscore,int M,int N,
	jz_protein *pt1,jz_protein *pt2,float abort_factor);

void	vote_lgc_dist(int *DIST,int *PTD,int N);

void	vote_lgc(int *LGC,int *DIST1,int *DIST2,int M,int N);

void	vote_trim_lgc(char *ALLOW,int *LGC,int M,int N);

void	vote_relation(Relation *RR,int *PTD,int N,int option);

void	vote_candidate(Voter *CANDS,int *ENDS,char *ALLOW,int M,int N);

void	vote_tally(Voter *CANDS,int *ENDS,Relation *RR1,
	Relation *RR2,int M,int N,void *chunk);

float	vote_degree(Voter *CANDS,int K);

int	vote_pseudo(Voter *CANDS,int K,float degree);

void	vote_matrix(int *VOTE,Voter *CANDS,int *ENDS,
	int M,int N,int L,int infinity);

void	vote_eliminate(int *VOTE,Voter *VC,int *ends,int K,int L,
	int M,int N);

int	vote_alignment(int *PAIRS,int *VOTE,int *LGC,Relation *RR1,
	Relation *RR2,int *PTD1,int *PTD2,int M,int N,int L,int *score);

int	vote_shrink(int *PAIRS,int len,int M,int N,
	Relation *RR1,Relation *RR2,int ***cache,int *map,int *pscore);

int	vote_refine(int *PAIRS,int len,int *VOTE,
	int *pt1,int *pt2,Relation *RR1,Relation *RR2,
	int M,int N,int *Score);

void	vote_follow(int *VOTE,Voter *CANDS,int nvoter,
	Relation *RR1,Relation *RR2,int M,int N,int leader);

void	vote_remove(Voter *VC,int nvoter,int *PAIRS,int len,int M,int N);

void	vote_print_matrix(int *MAT,int M,int N,int LD,char *fn);

void	vote_view_matrix(int *Mat,char *msg,int M,int N,int option);

void	vote_view_allow(char *Mat,char *msg,int M,int N);

float	**vote_RRD(jz_protein *pt);

float	vote_elastic(jz_pair_int *pairs,int len,float **RRD1,float **RRD2);

int	vote_score(int *pairs,int len,int M,int N,
	Relation *RR1,Relation *RR2);

int	vote_trim(jz_pair_int *pairs,int len,int limit);

int	vote_extend(jz_pair_int *pairs,int len,int M,int N);

void	vote_add(int *A,int M,int N);

int	vote_comittee(Voter *committee,int *ends,Voter *cands,
	int nvoter,int M);

void	vote_consolidate(Voter *VC,int nvoter);
void	vote_consolidate2(Voter **VC,int nvoter);

float	vote_zscore(float score,int M,int N);

void	vote_tally_fast(Voter *Voters,int nvoter,Relation *RR1,Relation *RR2,
	void *chunk,int M,int N,int L,float base,float fraction);
	
int	vote_db_save(jz_pair *protein_hash,char *fn);

jz_pair *vote_db_load(char *fn);

void	vote_cache_clean(int **cache,int M,int N);

void	vote_cache_initial(int ***cache,Relation *RR1,Relation *RR2,
	int *pairs,int len,int M,int N);

void	vote_map(int *map,int *pairs,int len,int *pairs_old,int len_old);

void	vote_cluster_pairs(FILE *fout,jz_pair *protein_hash,
	float len_ratio_max);

void	vote_cluster_identity(FILE *fout,FILE *fin,jz_pair *protein_hash);

int	vote_cluster_greedy(FILE *fout,FILE *fin,jz_pair *protein_hash,
	float len_limit,float score_min);

int	vote_cluster_transitive(FILE *fout,FILE *fin,jz_pair *protein_hash,
	float len_limit,float score_min);

int	vote_cluster_exact(FILE *fout,FILE *fin,jz_pair *protein_hash,
	float len_limit,float score_min,float time_limit);

int	vote_rms_fit(jz_pair_int *pairs,int N,jz_protein *pt1,jz_protein *pt2,
	float *C,float *R);

float	vote_rmsd(int *pairs,int len,jz_protein *pt1,jz_protein *pt2);

float	jz_rms_fitf(float *P,float *Q,int N,float *C,float *R);

#ifdef	__cplusplus
}
#endif

#endif


