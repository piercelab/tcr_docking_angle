/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#ifndef	__ALIGNMENT_H
#define	__ALIGNMENT_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "jz_protein.h"

typedef	short		residue_pair[2];

typedef	struct	Alignment
{
	residue_pair	*pairs;
	int		length;
	float		score;
	float		pvalue;
	float		rmsd;
}Alignment;

void	alignment_create(Alignment *al);
void	alignment_destroy(Alignment *al);
void	alignment_clear(Alignment *al);
void	alignment_resize(Alignment *al,int size);
void	alignment_copy(Alignment *dst,Alignment *src);
void	alignment_display_info(Alignment *al,jz_protein *p1,jz_protein *p2,
	FILE *fout);
void	alignment_display_simple(Alignment *al,jz_protein *p1,
	jz_protein *p2,FILE *fout);
void	alignment_display_rpl(Alignment *al,jz_protein *p1,
	jz_protein *p2,FILE *fout,int rpl);
void	alignment_trim_gaps(Alignment *al);
void	alignment_build(Alignment *al,int *pairs,int len);

#ifdef	__cplusplus
}
#endif

#endif


