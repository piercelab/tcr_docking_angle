/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
jz_protein.h		single chain protein 3D structure
=============================================================================*/

#ifndef	__JZ_PROTEIN_H
#define	__JZ_PROTEIN_H

#include <stdio.h>

#include "basic.h"

#ifdef	__cplusplus
extern "C" {
#endif

typedef	int16_t				jz_resid_pair[2];

#define	JZ_PROTEIN_ATOM			1
#define	JZ_PROTEIN_RESIDUE		2

#ifdef	JZ_PROTEIN_USE_SSE
#define	JZ_PROTEIN_SSE			4
#define	JZ_PROTEIN_CHAIN		8
#else
#define	JZ_PROTEIN_CHAIN		4
#endif

/*-----------------------------------------------------------------------------
jz_protein_atom corresponds to a ATOM record in a PDB file.
.atom is the atom code as given in /doc/atom_code.txt
.residue is the residue in which this atom is embeded, as defined in Amino.h.
This information is the same as the .name field in jz_proteinResidue
-----------------------------------------------------------------------------*/

typedef	struct		jz_protein_atom
{
	float		coord[3];	/* coordinate of the atom	*/
	char		name[2];	/* atom name			*/
	unsigned char	pad[2];		/* not used			*/
  char            first[31]; /* first part of the ATOM line */
  char            last[27];  /* last part of the ATOM line */
}jz_protein_atom;

typedef	jz_protein_atom	jz_atom;

/*-----------------------------------------------------------------------------
jz_protein_residue corresponds to an amino acid residue in a peptide chain.
-----------------------------------------------------------------------------*/

typedef	struct		jz_protein_residue
{
	float		center[3];	/* coordinate of the center	*/
	unsigned char	name;		/* residue code, 0-19		*/
	unsigned char	natom;		/* number of atoms		*/
	uint16_t	num;		/* original sequence number	*/
	jz_protein_atom	*atoms;		/* pointer to atoms		*/
}jz_protein_residue;

typedef	jz_protein_residue	jz_residue;

/*-----------------------------------------------------------------------------
jz_protein is a single chain protein.
-----------------------------------------------------------------------------*/

typedef	struct		jz_protein
{
	jz_protein_atom	*atoms;		/* all atoms			*/
	jz_protein_residue*residues;	/* all residues			*/
	int16_t		natom;		/* number of all atoms		*/
	int16_t		nres;		/* number of all residues	*/
}jz_protein;

/*-----------------------------------------------------------------------------
Basic operations on jz_protein and the components
-----------------------------------------------------------------------------*/

int	jz_protein_load(jz_protein *pt,char *fn,int mode);
  int	bp_protein_load(jz_protein *pt,const char *fn,int mode, char chn, int low_res, int high_res);
void	jz_protein_null(jz_protein *pt);
void	jz_protein_destroy(jz_protein *pt);
void	jz_protein_clear(jz_protein *pt);
void	jz_protein_copy(jz_protein *dst,jz_protein *src);
void	jz_protein_transform(jz_protein *dst,jz_protein *src,float *ori,
	float *rot,int option);
void	jz_protein_pseudo_atom_lines(jz_protein *pt,FILE *fout,int chain,
	int start);
void	jz_protein_output_atom_lines(jz_protein *pt,FILE *fout,int chain,
	int start);

/*-----------------------------------------------------------------------------
Macros to determine the type of a PDB line
-----------------------------------------------------------------------------*/

#define	JZ_PROTEIN_PDB_HEADER(line)	((line)[0]=='H'&&(line)[2]=='A')
#define	JZ_PROTEIN_PDB_HELIX(line)	((line)[0]=='H'&&(line)[2]=='L')
#define	JZ_PROTEIN_PDB_SHEET(line)	((line)[0]=='S'&&(line)[1]=='H')
#define	JZ_PROTEIN_PDB_TURN(line)	((line)[0]=='T'&&(line)[1]=='U')
#define	JZ_PROTEIN_PDB_ATOM(line)	((line)[0]=='A'&&(line)[1]=='T')
#define	JZ_PROTEIN_PDB_HETATM(line)	((line)[0]=='H'&&(line)[2]=='T')
#define	JZ_PROTEIN_PDB_MODEL(line)	((line)[0]=='M'&&(line)[1]=='O')
#define	JZ_PROTEIN_PDB_ENDMDL(line)	((line)[0]=='E'&&(line)[3]=='M')
#define	JZ_PROTEIN_PDB_TER(line)	((line)[0]=='T'&&(line)[1]=='E')
#define	JZ_PROTEIN_PDB_CHAIN(line)	((line)[0]=='C'&&(line)[1]=='H')	
#define	JZ_PROTEIN_PDB_CA(line)		((line)[13]=='C'&&(line)[14]=='A')


#ifdef	__cplusplus
}
#endif

#endif


