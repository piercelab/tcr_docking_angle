/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "jz_amino.h"

/* names as one-letter abbreviations */
char	jz_amino_char[]="acdefghiklmnpqrstvwy?-";

/* names as one-letter abbreviations */
char	jz_amino_CHAR[]="ACDEFGHIKLMNPQRSTVWY?-";

/* convert a single-letter residue name to 0-19 */
int	jz_amino_decode_char(int base)
{
	if(base>'o')
		if(base>'u')
			if(base>'x')return base-'a'-5;
			else return base-'a'-4;
		else return base-'a'-3;
	else if(base>'j')return base-'a'-2;
	else if(base>'b')return base-'a'-1;
	else if(base=='a')return 0;
	else if(base=='-')return 21;
	else return 20;
}

/* convert a three-letter residue name to 0-19 */
int	jz_amino_decode_name(char *name)
{
	char	a,b;
	
	a=name[0];
	if(a<='L')
	{
		switch(a)
		{
		case 'A':
			b=name[2];
			switch(b)
			{
			case 'A':return 0;
			case 'P':return 2;
			case 'N':return 11;
			case 'G':return 14;
			default:return 20;
			}
		case 'C':return 1;
		case 'G':
			b=name[2];
			switch(b)
			{
			case 'U':return 3;
			case 'Y':return 5;
			case 'N':return 13;
			default:return 20;
			}
		case 'H':return 6;
		case 'I':return 7;
		case 'L':
			b=name[2];
			if(b=='S')return 8;
			else if (b=='U')return 9;
		default:return 20;
		}
	}
	else
	{
		switch(a)
		{
		case 'M':return 10;
		case 'P':
			b=name[2];
			if(b=='E')return 4;
			else if(b=='O')return 12;
			else return 20;
		case 'S':return 15;
		case 'T':
			switch(b=name[1])
			{
			case 'H':return 16;
			case 'R':return 18;
			case 'Y':return 19;
			default:return 20;
			}
		case 'V':return 17;
		default:return 20;
		}
	}
}

/* three letter abbreviations */
char	*jz_amino_name[]=
	{
		"ALA","CYS","ASP","GLU","PHE",
		"GLY","HIS","ILE","LYS","LEU",
		"MET","ASN","PRO","GLN","ARG",
		"SER","THR","VAL","TRP","TYR",
		"?","-"
	};

/*-----------------------------------------------------------------------------
Number of core atoms (H excluded) in amino acid residues
Use this to check if a PDB record is complete
-----------------------------------------------------------------------------*/

int8_t	jz_amino_atom_count[]=
{
	 5, /* ALA */
	 6, /* CYS */
	 8, /* ASP */
	 9, /* GLU */
	11, /* PHE */
	 4, /* GLY */
	10, /* HIS */
	 8, /* ILE */
	 9, /* LYS */
	 8, /* LEU */
	 8, /* MET */
	 8, /* ASN */
	 7, /* PRO */
	 9, /* GLN */
	11, /* ARG */
	 6, /* SER */
	 7, /* THR */
	 7, /* VAL */
	14, /* TRP */
	12  /* TYR */
};


