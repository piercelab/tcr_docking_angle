/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
amino.h			definitions of amino acid residues
=============================================================================*/

#ifndef	__JZ_AMINO_H
#define	__JZ_AMINO_H

#include <basic.h>

#ifdef	__cplusplus
extern	"C" {
#endif

/*-----------------------------------------------------------------------------
amino acids representation
-----------------------------------------------------------------------------*/

#define	JZ_AMINO_ALA	0
#define	JZ_AMINO_CYS	1
#define	JZ_AMINO_ASP	2
#define	JZ_AMINO_GLU	3
#define	JZ_AMINO_PHE	4

#define	JZ_AMINO_GLY	5
#define	JZ_AMINO_HIS	6
#define	JZ_AMINO_ILE	7
#define	JZ_AMINO_LYS	8
#define	JZ_AMINO_LEU	9

#define	JZ_AMINO_MET	10
#define	JZ_AMINO_ASN	11
#define	JZ_AMINO_PRO	12
#define	JZ_AMINO_GLN	13
#define	JZ_AMINO_ARG	14

#define	JZ_AMINO_SER	15
#define	JZ_AMINO_THR	16
#define	JZ_AMINO_VAL	17
#define	JZ_AMINO_TRP	18
#define	JZ_AMINO_TYR	19

#define	JZ_AMINO_NOP	20
#define	JZ_AMINO_GAP	21

extern	char	jz_amino_CHAR[];
extern	char	jz_amino_char[];
extern	char	*jz_amino_name[];

extern	int8_t		jz_amino_atom_count[];
int	jz_amino_decode_name(char *name);
int	jz_amino_decode_char(int base);

#ifdef	__cplusplus
}
#endif

#endif


