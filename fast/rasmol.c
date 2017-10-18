/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
rasmol.c	Generate rasmol scripts to view protein structure alignments
=============================================================================*/

#include "basic.h"
#include "rasmol.h"

/*-----------------------------------------------------------------------------
 * Print header for the rasmol script.
 * This routine was copied from Joe's cartoon_rasmol.cc
-----------------------------------------------------------------------------*/

void	vote_rasmol_boilerplate(FILE *fout,int skinny)
{
	fprintf(fout, "load inline\nbackground white\nset stereo off\n\
set ambient 60\n\n");
	fprintf(fout, "echo ###### FORMAT ######\n");
	fprintf(fout, "echo thick backbone denotes structurally\n");
	fprintf(fout, "echo equivalent regions.\nselect all\n");
	fprintf(fout, "color bonds none\ncolor backbone none\n\
color hbonds none\n");
	fprintf(fout, "color ssbonds none\ncolor ribbons none\n\
spacefill off\n");
	fprintf(fout, "wireframe off\nribbons off\nbackbone %d\n", skinny);
	fprintf(fout, "select :A\nbackbone %d\ncolor red\nselect :B\n\
backbone %d\ncolor cyan\n\n", skinny, skinny);
}

/*-----------------------------------------------------------------------------
 * highlight
-----------------------------------------------------------------------------*/

void	vote_rasmol_alignment_highlight(jz_pair_int *pairs,int len,
	jz_protein_residue *reses,FILE *fout,int chain)
{
	int		i,j,k,k_start;
	int		seq1,seq2,c1,c2;
	
	for(k_start=0,k=1;k<len;k++)
	{
		i=pairs[k][0];
		j=pairs[k][1];
		if((i-pairs[k-1][0]!=1)||(j-pairs[k-1][1]!=1))
		{
			if(chain=='A')
			{
				seq1=reses[pairs[k_start][0]].num;
				seq2=reses[pairs[k-1][0]].num;
			}
			else
			{
				seq1=reses[pairs[k_start][1]].num;
				seq2=reses[pairs[k-1][1]].num;
			}
			c1=seq1>>12;
			if(c1)c1='A'-1+c1;
			else c1=' ';
			c2=seq2>>12;
			if(c2)c2='A'-1+c2;
			else c2=' ';
			seq1&=4095;
			seq2&=4095;

			fprintf(fout,"select %d",seq1);
			if(c1!=' ')fputc(c1,fout);
			fprintf(fout,"-%d",seq2);
			if(c2!=' ')fputc(c2,fout);
			fprintf(fout,":%c\nbackbone 80\n",chain);

			k_start=k;
		}
	}

	k=len-1;
	if(chain=='A')
	{
		seq1=reses[pairs[k_start][0]].num;
		seq2=reses[pairs[k-1][0]].num;
	}
	else
	{
		seq1=reses[pairs[k_start][1]].num;
		seq2=reses[pairs[k-1][1]].num;
	}
	c1=seq1>>12;
	if(c1)c1='A'-1+c1;
	else c1=' ';
	c2=seq2>>12;
	if(c2)c2='A'-1+c2;
	else c2=' ';
	seq1&=4095;
	seq2&=4095;

	fprintf(fout,"select %d",seq1);
	if(c1!=' ')fputc(c1,fout);
	fprintf(fout,"-%d",seq2);
	if(c2!=' ')fputc(c2,fout);
	fprintf(fout,":%c\nbackbone 80\n",chain);
}

/*-----------------------------------------------------------------------------
 * print a rasmol script to display an Residue-Residue alignment.
-----------------------------------------------------------------------------*/

void	vote_rasmol_pairs(jz_pair_int *pairs,int len,
	jz_protein *pt1,jz_protein *pt2,int skinny,FILE *fout)
{
	jz_protein	clone2,*ppt2;
	float		ori[3],rot[9];

	extern	int	vote_rms_fit(jz_pair_int *pairs,int N,
			jz_protein *pt1,jz_protein *pt2,float *C,float *R);

	jz_protein_null(&clone2);
	if(vote_rms_fit(pairs,len,pt1,pt2,ori,rot))ppt2=pt2;
	else
	{
		jz_protein_copy(&clone2,pt2);
		ppt2=&clone2;
		jz_protein_transform(ppt2,pt2,ori,rot,255);
	}
	vote_rasmol_boilerplate(fout,skinny);
	vote_rasmol_alignment_highlight(pairs,len,pt1->residues,fout,'A');
	vote_rasmol_alignment_highlight(pairs,len,pt2->residues,fout,'B');
	fprintf(fout,"\nexit\n");
	jz_protein_pseudo_atom_lines(pt1,fout,'A',0);
	jz_protein_pseudo_atom_lines(ppt2,fout,'B',0);
	/*if(ppt2==&clone2)jz_protein_destroy(ppt2);*/
}

/*-----------------------------------------------------------------------------
 * print a rasmol script to display an Residue-Residue alignment.
-----------------------------------------------------------------------------*/

void	vote_rasmol_alignment(Alignment *al,jz_protein *pt1,jz_protein *pt2,
	FILE *fout)
{
	jz_pair_int	*pairs;
	int		k,len;
	
	if(!al||al->length<=0)return;
	len=al->length;
	JZ_ARRAY_INIT(pairs,len);
	for(k=0;k<len;k++)
		pairs[k][0]=al->pairs[k][0],
		pairs[k][1]=al->pairs[k][1];
	vote_rasmol_pairs(pairs,len,pt1,pt2,10,fout);
	JZ_ARRAY_FREE(pairs);
}

/*-----------------------------------------------------------------------------
 * print an aligned protein2 file
-----------------------------------------------------------------------------*/

void	vote_output_coords(Alignment *al,jz_protein *pt1,jz_protein *pt2,
	FILE *fout)
{
	jz_pair_int	*pairs;
	int		k,len;
	jz_protein	clone2,*ppt2;
	float		ori[3],rot[9];

	extern	int	vote_rms_fit(jz_pair_int *pairs,int N,
			jz_protein *pt1,jz_protein *pt2,float *C,float *R);

	if(!al||al->length<=0)return;
	len=al->length;
	JZ_ARRAY_INIT(pairs,len);
	for(k=0;k<len;k++)
		pairs[k][0]=al->pairs[k][0],
		pairs[k][1]=al->pairs[k][1];

	jz_protein_null(&clone2);
	if(vote_rms_fit(pairs,len,pt1,pt2,ori,rot))ppt2=pt2;
	else
	{
		jz_protein_copy(&clone2,pt2);
		ppt2=&clone2;
		jz_protein_transform(ppt2,pt2,ori,rot,255);
	}

	jz_protein_output_atom_lines(ppt2,fout,'B',0);

	JZ_ARRAY_FREE(pairs);
}


