/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
rasmol.h	Generate rasmol scripts to view protein structure alignments
=============================================================================*/

#ifndef	__RASMOL_H__
#define	__RASMOL_H__

#include "Alignment.h"

#ifdef	__cplusplus
extern "C" {
#endif

void	vote_rasmol_alignment(Alignment *al,jz_protein *pt1,jz_protein *pt2,
	FILE *fout);

void    vote_output_coords(Alignment *al,jz_protein *pt1,jz_protein *pt2,
	FILE *fout);


#ifdef	__cplusplus
}
#endif

#endif

