/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
misc.h		solve some basic math problems.
=============================================================================*/

#ifndef	MISC_H
#define	MISC_H

#include "basic.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
basic operations. the final result always goes to the first argument.
-----------------------------------------------------------------------------*/

#define	JZ_VEC3_SET(V,a)						\
	((V)[0]=(V)[1]=(V)[2]=(a))

#define	JZ_VEC3_CPY(dest,src)						\
	((dest)[0]=(src)[0],(dest)[1]=(src)[1],(dest)[2]=(src)[2])

#define	JZ_VEC3_SUB(a,b,c)						\
	((a)[0]=(b)[0]-(c)[0],(a)[1]=(b)[1]-(c)[1],(a)[2]=(b)[2]-(c)[2])

#define	JZ_VEC3_NSUB(r,p)						\
	((r)[0]=(p)[0]-(p)[3],(r)[1]=(p)[1]-(p)[4],(r)[2]=(p)[2]-(p)[5])

#define	JZ_VEC3_ADD(a,b,c)						\
	((a)[0]=(b)[0]+(c)[0],(a)[1]=(b)[1]+(c)[1],(a)[2]=(b)[2]+(c)[2])
	
#define	JZ_VEC3_NADD(r,p)						\
	((r)[0]=(p)[0]+(p)[3],(r)[1]=(p)[1]+(p)[4],(r)[2]=(p)[2]+(p)[5])

#define	JZ_VEC3_MUL(p,q,a)						\
	((p)[0]=(q)[0]*(a),(p)[1]=(q)[1]*(a),(p)[2]=(q)[2]*(a))

#define	JZ_VEC3_DIV(p,q,a)						\
	((p)[0]=(q)[0]/(a),(p)[1]=(q)[1]/(a),(p)[2]=(q)[2]/(a))
	
#define	JZ_VEC3_MULA(p,q,r,a)						\
	((p)[0]=(q)[0]*a+(r)[0],(p)[1]=(q)[1]*a+(r)[1],(p)[2]=(q)[2]*a+(r)[2])
	
#define	JZ_VEC3_MULS(p,q,r,a)						\
	((p)[0]=(q)[0]*a-(r)[0],(p)[1]=(q)[1]*a-(r)[1],(p)[2]=(q)[2]*a-(r)[2])

#define	JZ_VEC3_DOT(d,a,b)						\
	((d)=(((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])))
	
#define	JZ_VEC3_NDOT(d,a)						\
	((d)=(((a)[0]*(a)[3]+(a)[1]*(a)[4]+(a)[2]*(a)[5])))
	
#define	JZ_VEC3_SQR(d,a)						\
	((d)=(((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2])))

#define	JZ_VEC3_CROS(cross,a,b)					\
	((cross)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1],\
	(cross)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2],\
	(cross)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0])
	
#define	JZ_VEC3_NCROS(cross,a)					\
	((cross)[0]=(a)[1]*(a)[5]-(a)[2]*(a)[4],\
	(cross)[1]=(a)[2]*(a)[3]-(a)[0]*(a)[5],\
	(cross)[2]=(a)[0]*(a)[4]-(a)[1]*(a)[3])
	
#define	JZ_VEC3_SUBSQ(d,p,q)						\
	((d)=(((p)[0]-(q)[0])*((p)[0]-(q)[0])\
	+((p)[1]-(q)[1])*((p)[1]-(q)[1])\
	+((p)[2]-(q)[2])*((p)[2]-(q)[2])))
	
#define	JZ_VEC3_NSUBSQ(d,p)						\
	((d)=((p)[0]-(p)[3])*((p)[0]-(p)[3])\
	+((p)[1]-(p)[4])*((p)[1]-(p)[4])\
	+((p)[2]-(p)[5])*((p)[2]-(p)[5]))
	
#define	JZ_VEC3_N2SUBSQ(d,p,q)						\
	((d)= ((p)[0]-(q)[0])*((p)[0]-(q)[0])\
	+((p)[1]-(q)[1])*((p)[1]-(q)[1])\
	+((p)[2]-(q)[2])*((p)[2]-(q)[2])\
	+((p)[3]-(q)[3])*((p)[3]-(q)[3])\
	+((p)[4]-(q)[4])*((p)[4]-(q)[4])\
	+((p)[5]-(q)[5])*((p)[5]-(q)[5]))

#define	JZ_VEC3_ROT(A,B,rot)						\
	((A)[0]=(B)[0]*(rot)[0]+(B)[1]*(rot)[1]+(B)[2]*(rot)[2],\
	(A)[1]=(B)[0]*(rot)[3]+(B)[1]*(rot)[4]+(B)[2]*(rot)[5],\
	(A)[2]=(B)[0]*(rot)[6]+(B)[1]*(rot)[7]+(B)[2]*(rot)[8])
	
#define	JZ_VEC3_RMUL(A,mat,B)						\
	JZ_VEC3_ROT(A,B,mat)
	
#define	JZ_VEC3_COS(c,p,q,d)						\
	do\
	{\
		JZ_VEC3_DOT((d),(q),(q));\
		JZ_VEC3_DOT((c),(p),(p));\
		(d)=(d)*(c);\
		(d)=sqrt(d);\
		JZ_VEC3_DOT((c),(p),(q));\
		(c)/=(d);\
	}while(0)


int	jz_select_i(int *,int,int);

int	jz_dp_local_affine(int *scores,int M,int N,
	int gap_a,int gap_b,jz_pair_int *align,
	int align_gap,int *score);

int	jz_dp_local_general_affine(int *scores,int M,int N,
	int gap_a,int gap_b,int gap_c,jz_pair_int *align,
	int align_gap,int *score);

void	jz_dp_local_affine_potential(int *P,int *scores,
	int M,int N,int gap_a,int gap_b);

void	jz_dp_local_affine_diag(int *A,int *scores,int M,int N,
	int gap_a,int gap_b);

int32_t	jz_intarray_max(int32_t *array,size_t N);

int	jz_charray_count_zero(char *array,size_t N);

void	jz_intarray_set(int32_t *array,int32_t value,size_t N);

#define	jz_charray_set	memset

#ifdef	__cplusplus
}
#endif

#endif

