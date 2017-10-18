/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "basic.h"
#include "chore.h"
#include "jz_amino.h"
#include "jz_protein.h"
#include "Alignment.h"

/*-----------------------------------------------------------------------------
Create as invalid
-----------------------------------------------------------------------------*/

void	alignment_create(Alignment *al)
{
	al->pairs=NULL;
	al->length=0;
	al->score=0.0;
	al->pvalue=1.0;
}

/*-----------------------------------------------------------------------------
Free dynamically allocated memory, and initialize to reuse
-----------------------------------------------------------------------------*/

void	alignment_clear(Alignment *al)
{
	if(al->pairs)free(al->pairs);
	al->pairs=NULL;
	al->length=0;
	al->score=0.0;
	al->pvalue=1.0;
}

/*-----------------------------------------------------------------------------
Free memory
-----------------------------------------------------------------------------*/

void	alignment_destroy(Alignment *al)
{
	if(al->pairs)free(al->pairs);
}

/*-----------------------------------------------------------------------------
Adjust the length
-----------------------------------------------------------------------------*/

void	alignment_resize(Alignment *al,int size)
{
	jz_resid_pair	*pairs;

	if(al->pairs)JZ_ARRAY_FREE(al->pairs);
	if(size<=0)pairs=NULL;
	else JZ_ARRAY_INIT(al->pairs,size);
	al->length=0;
	al->score=0;
	al->pvalue=1.0;
}

/*-----------------------------------------------------------------------------
Make a copy of the alignment data
-----------------------------------------------------------------------------*/

void	alignment_copy(Alignment *dst,Alignment *src)
{
	*dst=*src;
	if(dst->pairs)
		JZ_ARRAY_DUP(dst->pairs,src->pairs,dst->length);
}

/*-----------------------------------------------------------------------------
Display the alignment in a simple format
-----------------------------------------------------------------------------*/

void	alignment_display_simple(Alignment *al,jz_protein *p1,jz_protein *p2,
	FILE *fout)
{
	int		k;
	int		m,n,c1,c2;
	jz_resid_pair	*pairs;
	
	if(!fout)return;

	pairs=al->pairs;
	if(al->length<=0||pairs==NULL)return;
	
	for(k=0;k<al->length;k++)
	{
		m=p1->residues[pairs[k][0]].num;
		n=p2->residues[pairs[k][1]].num;
		c1=m>>12;
		if(c1)c1='A'-1+c1;else c1=' ';
		c2=n>>12;
		if(c2)c2='A'-1+c2;else c2=' ';

		fprintf(fout,"%5d%c %s %c  %c %s %5d%c\n",
			m&4095,c1,
			jz_amino_name[p1->residues[pairs[k][0]].name],
			jz_amino_CHAR[p1->residues[pairs[k][0]].name],
			jz_amino_CHAR[p2->residues[pairs[k][1]].name],
			jz_amino_name[p2->residues[pairs[k][1]].name],
			n&4095,c2);
	}
}

/*-----------------------------------------------------------------------------
display alignment in an line by line format.
-----------------------------------------------------------------------------*/

void	alignment_display_rpl(Alignment *al,jz_protein *pt1,jz_protein *pt2,
	FILE *fout,int rpl)
{
	int		i,j,k,m,n,M,N;
	jz_resid_pair	*pairs;
	char		*A,*B;
	
	if(!fout)return;
	pairs=al->pairs;
	if(al->length<=0||pairs==NULL)return;

	M=pt1->nres;
	N=pt2->nres;
	JZ_ARRAY_INIT(A,2*(M+N+8));
	B=A+M+N+4;
	
	for(m=n=-1,i=k=0;i<al->length;i++)
	{
		j=m+1;m=pairs[i][0];
		for(;j<m;j++)
		{
			A[k]=jz_amino_CHAR[pt1->residues[j].name];
			B[k++]='-';
		}

		j=n+1;n=pairs[i][1];
		for(;j<n;j++)
		{
			A[k]='-';
			B[k++]=jz_amino_CHAR[pt2->residues[j].name];
		}
		
		A[k]=jz_amino_CHAR[pt1->residues[m].name];
		B[k++]=jz_amino_CHAR[pt2->residues[n].name];
	}
	
	for(j=m+1;j<M;j++)
	{
		A[k]=jz_amino_CHAR[pt1->residues[j].name];
		B[k++]='-';
	}

	for(j=n+1;j<N;j++)
	{
		A[k]='-';
		B[k++]=jz_amino_CHAR[pt2->residues[j].name];
	}
	
	A[k]='*';
	B[k]='*';
	k++;
	A[k]=0;
	B[k]=0;

	/* now display */
	for(i=0;i<k;i+=rpl)
	{
		if(i+rpl<=k)j=rpl;else j=k-i;
		fputs("\n 1:\t",fout);
		fwrite(A+i,1,j,fout);
		fputs("\n 2:\t",fout);
		fwrite(B+i,1,j,fout);
		fputc('\n',fout);
	}
	fputc('\n',fout);

	JZ_ARRAY_FREE(A);
}

/*-----------------------------------------------------------------------------
Display basic information of the alignment
-----------------------------------------------------------------------------*/

void	alignment_display_info(Alignment *al,jz_protein *p1,jz_protein *p2,
	FILE *fout)
{
	if(!fout)return;
	fprintf(fout,"L=%-d SX=%-.3e SN=%-.3e L1=%d L2=%d RMSD=%-.3f\n",
		al->length,al->score,al->pvalue,p1->nres,p2->nres,al->rmsd);
}

/*-----------------------------------------------------------------------------
Remove gaps from an alignment
-----------------------------------------------------------------------------*/

void	alignment_trim_gaps(Alignment *al)
{
	int		i,j;
	jz_resid_pair	*pairs;

	if(!(pairs=al->pairs))return;
	for(i=j=0;j<al->length;j++)
	{
		if(pairs[j][0]<0||pairs[j][1]<0)continue;
		pairs[i][0]=pairs[j][0];
		pairs[i++][1]=pairs[j][1];
	}
	al->length=i;
}

/*-----------------------------------------------------------------------------
Create an alignment from a pairs array. The dynamic programming routine
return array of int pairs. This routine convert the arrary into array
of jz_resid_pairs
-----------------------------------------------------------------------------*/

void	alignment_build(Alignment *al,int *pairs,int len)
{
	int		k;

	if(al->pairs)
	{
		JZ_ARRAY_FREE(al->pairs);
		al->pairs=NULL;
	}
	al->length=0;
	al->score=0;
	al->pvalue=1.0;
	if(len<=0)return;
	JZ_ARRAY_INIT(al->pairs,len);
	al->length=len;
	for(k=0;k<len;k++)
	{
		al->pairs[k][0]=pairs[k+k];
		al->pairs[k][1]=pairs[k+k+1];
	}
}

