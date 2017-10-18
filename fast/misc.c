/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <limits.h>

#define	SWAP(a,b,c)						\
	((c)=(a),(a)=(b),(b)=(c))

int	jz_select_i(int *array,int n,int pos)
{
	int	i,j,k,l,r;
	int	buf,s;
	
	l=0;r=n-1;
	while(l<r)
	{
		i = l+(r-l)/2;
		if (array[i]>array[r])
			SWAP(array[i],array[r],buf);
		SWAP(array[l],array[i],buf);
		i=l;
		k=r;
		s=array[l];
		for(;;)
		{
			while (array[(++i)]<s);
			while (array[--k]>s);
			if(i<k)SWAP(array[i],array[k],buf);
			else break;
		}
		SWAP(array[l],array[k],buf);
		j=k-l+1;
		if(pos<=j)r=k;
		else
		{
			l=k+1;
			pos-=j;
		}
	}
	return array[l];
}

/******************************************************************************
jz_dp_local__affine.c		DP with generalized affine gap 
******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "basic.h"

/*-----------------------------------------------------------------------------
 * Local alignment dynamic programming, affine gaps.
 * This is an implementation by definition, hence might be inefficient.
 * See research note pp5
-----------------------------------------------------------------------------*/

int	jz_dp_local_affine(int *scores,int M,int N,int gap_a,
	int gap_b,jz_pair_int *align,int align_gap,int *score)
{
	int		**A,**B,**C,**P,a,max_score;
	int		i,j,k,max_i,max_j;
	
	if(!scores||M<=0||N<=0)return -1;
	
	/* Allocate DP matrices */
	JZ_MATRIX_INITL(A,-1,M-1,-1,N-1);
	JZ_MATRIX_INITL(B,-1,M-1,-1,N-1);
	JZ_MATRIX_INITL(C,-1,M-1,-1,N-1);
	
	/* Initialize DP matrices */
	for(k=-1;k<N;k++)
		A[-1][k]=B[-1][k]=C[-1][k]=0;
	for(k=-1;k<M;k++)
		A[k][-1]=B[k][-1]=C[k][-1]=0;

	/* build the DP matrices */
	max_score=0;
	max_i=max_j=-1;
	for(i=0;i<M;i++)
	for(j=0;j<N;j++)
	{
		a=A[i-1][j-1];
		if(B[i-1][j-1]>a)a=B[i-1][j-1];
		if(C[i-1][j-1]>a)a=C[i-1][j-1];
		a+=*scores++;
		if(a<0)a=0;
		A[i][j]=a;
		if(a>max_score)
		{
			max_score=a;
			max_i=i;
			max_j=j;
		}

		a=A[i][j-1]+gap_a;
		if(B[i][j-1]>a)a=B[i][j-1];
		if(C[i][j-1]>a)a=C[i][j-1];
		a+=gap_b;
		if(a<0)a=0;
		B[i][j]=a;

		a=A[i-1][j]+gap_a;
		if(B[i-1][j]>a)a=B[i][j-1];
		if(C[i-1][j]>a)a=C[i-1][j];
		a+=gap_b;
		if(a<0)a=0;
		C[i][j]=a;
	}

	/* trace back */
	i=max_i,j=max_j;k=M+N-1;
	for(P=A;P[i][j]>0;)
	{
		if(P==A)
		{
			align[k][0]=i--;
			align[k--][1]=j--;
		}
		else if(P==B)
		{
			if(align_gap)
			{
				align[k][0]=-1;
				align[k--][1]=j;
			}
			j--;
		}
		else if(P==C)
		{
			if(align_gap)
			{
				align[k][0]=i;
				align[k--][1]=-1;
			}
			i--;
		}
		a=A[i][j];
		if(P!=A)a+=gap_a,P=A;
		if(B[i][j]>a)a=B[i][j],P=B;
		if(C[i][j]>a)a=C[i][j],P=C;
	}

	memmove(align,align+k+1,(M+N-k-1)*sizeof(jz_pair_int));
	JZ_MATRIX_FREEL(A,-1);
	JZ_MATRIX_FREEL(B,-1);
	JZ_MATRIX_FREEL(C,-1);
	*score=max_score;
	k=M+N-k-1;
	return k;
}


/*-----------------------------------------------------------------------------
Potential matrix P = A + A' - S, where
A  is the forward diagonal dynamic programming matrix
A' is the backward diagonal dynamic programming matrix
S  is the scoring matrix

P[i][j] is the maximum total score if i is paired to j.
-----------------------------------------------------------------------------*/

void	jz_dp_local_affine_potential(int *P,int *scores,
	int M,int N,int gap_a,int gap_b)
{
	int	*B,*S;
	int		i,j;
	
	if(!scores||M<=0||N<=0)return;
	JZ_ARRAY_INIT(S,M*N);
	JZ_ARRAY_INIT(B,M*N);

	for(i=0,j=M*N-1;i<=j;i++,j--)
		S[i]=scores[j],S[j]=scores[i];

	jz_dp_local_affine_diag(P,scores,M,N,gap_a,gap_b);
	jz_dp_local_affine_diag(B,S,M,N,gap_a,gap_b);
	for(i=0,j=M*N-1;j>=0;i++,j--,scores++)P[i]+=B[j]-*scores;
	
	JZ_ARRAY_FREE(S);
	JZ_ARRAY_FREE(B);
}


/*-----------------------------------------------------------------------------
 * Local alignment dynamic programming, with affine gap penalties.
 * This is to get the Diagonal matrix A.
-----------------------------------------------------------------------------*/

void	jz_dp_local_affine_diag(int *A,int *scores,int M,int N,
	int gap_a,int gap_b)
{
	int	*T,*B,*C,a;
	int		j;
	
	if(!scores||M<=0||N<=0)return;
	
	/* Allocate DP matrices */
	j=M*N;
	JZ_ARRAY_INIT(T,j+j);
	B=T;C=B+j;
	
	/* build the DP matrices */

	/* A[0][0] */
	a=*scores++;
	if(a<0)a=0;
	*A=a;
	*B=0;
	*C=0;

	/* A[0][j] */
	for(j=N-1;j;j--)
	{
		a=*scores++;
		if(a<0)a=0;
		A[1]=a;
		
		a=A[0]+gap_a;
		if(B[0]>a)a=B[0];
		if(C[0]>a)a=C[0];
		a+=gap_b;
		if(a<0)a=0;
		A++,B++,C++;
		B[0]=a;

		C[0]=0;
	}

	/* A[i][j], where i>0 */
	for(--M;M;M--)
	{
		/* A[i][0] */
		a=*scores++;
		if(a<0)a=0;
		A++,B++,C++;
		*A=a;
		*B=0;
		
		a=A[-N]+gap_a;
		if(B[-N]>a)a=B[-N];
		if(C[-N]>a)a=C[-N];
		a+=gap_b;
		if(a<0)a=0;
		C[0]=a;
		
		for(j=N-1;j;j--)
		{
			a=A[-N];
			if(B[-N]>a)a=B[-N];
			if(C[-N]>a)a=C[-N];
			a+=*scores++;
			if(a<0)a=0;
			A[1]=a;

			a=A[0]+gap_a;
			if(B[0]>a)a=B[0];
			if(C[0]>a)a=C[0];
			a+=gap_b;
			if(a<0)a=0;
			A++,B++,C++;
			B[0]=a;

			a=A[-N]+gap_a;
			if(B[-N]>a)a=B[-N];
			if(C[-N]>a)a=C[-N];
			a+=gap_b;
			if(a<0)a=0;
			C[0]=a;
		}
	}	

	JZ_ARRAY_FREE(T);
}


int32_t	jz_intarray_max(int32_t *array,size_t N)
{
	size_t		i,n;
	int32_t		m;
	
	m=INT_MIN;
	
	n=N&~3;
	for(i=0;i<n;i+=4)
	{
		if(array[i]>m)m=array[i];
		if(array[i+1]>m)m=array[i+1];
		if(array[i+2]>m)m=array[i+2];
		if(array[i+3]>m)m=array[i+3];
	}
	for(;i<N;i++)
		if(array[i]>m)m=array[i];

	return m;
}

/*-----------------------------------------------------------------------------
get the number of zero's in a char array
-----------------------------------------------------------------------------*/

int	jz_charray_count_zero(char *array,size_t N)
{
	int		k;
	int		count=0;
	
	for(k=0;k<N;k++)
		if(array[k]==0)count++;
	
	return count;
}


void	jz_intarray_set(int32_t *array,int32_t value,size_t N)
{
	size_t		k;
	
	for(k=N>>2;k;k--,array+=4)
		array[0]=array[1]=array[2]=array[3]=value;
	if(N&1)
	{
		array[0]=value;
		array++;
	}	
	if(N&2)array[0]=array[1]=value;
}


/*-----------------------------------------------------------------------------
 * Local alignment dynamic programming, with generalized affine gaps.
 * This is an implementation by definition, hence might be inefficient.
 * See research note pp5
-----------------------------------------------------------------------------*/

int	jz_dp_local_general_affine(int *scores,int M,int N,
	int gap_a,int gap_b,int gap_c,jz_pair_int *align,
	int align_gap,int *score)
{
	/* gap cost = a + b|k1-k2| + c min(k1,k2) where k1 and k2 are gaps
	 * from both sequences respectively. vector k2 and k2 stores
	 * the number of gaps for the current row. */

	int	**A,**B,**C,**D,**P,a,max_score;
	int		i,j,k,max_i,max_j;
	
	if(!scores||M<=0||N<=0)return -1;
	
	/* Allocate DP matrices */
	JZ_MATRIX_INITL(A,-1,M-1,-1,N-1);
	JZ_MATRIX_INITL(B,-1,M-1,-1,N-1);
	JZ_MATRIX_INITL(C,-1,M-1,-1,N-1);
	JZ_MATRIX_INITL(D,-1,M-1,-1,N-1);
	
	/* Initialize DP matrices */
	for(k=-1;k<N;k++)
		A[-1][k]=B[-1][k]=C[-1][k]=D[-1][k]=0;
	for(k=-1;k<M;k++)
		A[k][-1]=B[k][-1]=C[k][-1]=D[k][-1]=0;

	/* build the DP matrices */
	max_score=0;
	max_i=max_j=-1;
	for(i=0;i<M;i++)
	for(j=0;j<N;j++)
	{
		a=A[i-1][j-1];
		if(B[i-1][j-1]>a)a=B[i-1][j-1];
		if(C[i-1][j-1]>a)a=C[i-1][j-1];
		if(D[i-1][j-1]>a)a=D[i-1][j-1];
		a+=*scores++;
		if(a<0)a=0;
		A[i][j]=a;
		if(a>max_score)
		{
			max_score=a;
			max_i=i;
			max_j=j;
		}

		a=A[i][j-1]+gap_a;
		if(B[i][j-1]>a)a=B[i][j-1];
		if(C[i][j-1]>a)a=C[i][j-1];
		if(D[i][j-1]>a)a=D[i][j-1];
		a+=gap_b;
		if(a<0)a=0;
		B[i][j]=a;

		a=A[i-1][j]+gap_a;
		if(B[i][j-1]>a)a=B[i][j-1];
		if(C[i-1][j]>a)a=C[i-1][j];
		if(D[i-1][j]>a)a=D[i-1][j];
		a+=gap_b;
		if(a<0)a=0;
		C[i][j]=a;

		a=A[i-1][j-1]+gap_a;
		if(B[i-1][j-1]>a)a=B[i-1][j-1];
		if(C[i-1][j-1]>a)a=C[i-1][j-1];
		if(D[i-1][j-1]>a)a=D[i-1][j-1];
		a+=gap_c;
		if(a<0)a=0;
		D[i][j]=a;
	}

	/* trace back */
	i=max_i,j=max_j;k=M+N-1;
	for(P=A;P[i][j]>0;)
	{
		if(P==A)
		{
			align[k][0]=i--;
			align[k--][1]=j--;
		}
		else if(P==B)
		{
			if(align_gap)
			{
				align[k][0]=-1;
				align[k--][1]=j;
			}
			j--;
		}
		else if(P==C)
		{
			if(align_gap)
			{
				align[k][0]=i;
				align[k--][1]=-1;
			}
			i--;
		}
		else
		{
			if(align_gap)
			{
				align[k][0]=i;
				align[k--][1]=-2;
				align[k][0]=-2;
				align[k--][1]=j;
			}
			i--,j--;
		}
		a=A[i][j];
		if(P!=A)a+=gap_a,P=A;
		if(B[i][j]>a)a=B[i][j],P=B;
		if(C[i][j]>a)a=C[i][j],P=C;
		if(D[i][j]>a)a=D[i][j],P=D;
	}
	memmove(align,align+k+1,(M+N-k-1)*sizeof(jz_pair_int));
	JZ_MATRIX_FREEL(A,-1);
	JZ_MATRIX_FREEL(B,-1);
	JZ_MATRIX_FREEL(C,-1);
	JZ_MATRIX_FREEL(D,-1);
	*score=max_score;
	k=M+N-k-1;
	return k;
}


