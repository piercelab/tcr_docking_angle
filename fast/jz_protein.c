/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <string.h>

#include "basic.h"
#include "chore.h"
#include "misc.h"
#include "jz_amino.h"
#include "jz_protein.h"


/*-----------------------------------------------------------------------------
clear
-----------------------------------------------------------------------------*/

void	jz_protein_clear(jz_protein *pt)
{
	jz_protein_destroy(pt);
	memset(pt,0,sizeof(jz_protein));
}

/*-----------------------------------------------------------------------------
make a copy of protein src to dst
-----------------------------------------------------------------------------*/

void	jz_protein_copy(jz_protein *dst,jz_protein *src)
{
	int		k;
	int		d;
	char		**p;
	
	/* copy sizes */
	dst->natom=src->natom;
	dst->nres=src->nres;

	/* duplicate atoms */
	if((dst->atoms=src->atoms))
		JZ_ARRAY_DUP(dst->atoms,src->atoms,dst->natom);
	
	/* duplicate residues */
	if((dst->residues=src->residues))
	{	
		JZ_ARRAY_DUP(dst->residues,src->residues,(int)dst->nres);
		if(dst->atoms)
		{
			d=(char*)dst->atoms-(char*)src->atoms;
			p=(char**)(&(dst->residues->atoms));
			for(k=dst->nres;k;k--)
			{
				if(*p)*p+=d;
				p=(char**)((char*)p
					+sizeof(jz_protein_residue));
			}
		}
	}
}


/*-----------------------------------------------------------------------------
destroy a protein instance
-----------------------------------------------------------------------------*/

void	jz_protein_destroy(jz_protein *pt)
{
	if(pt->atoms)JZ_ARRAY_FREE(pt->atoms);
	if(pt->residues)JZ_ARRAY_FREE(pt->residues);
}

/*-----------------------------------------------------------------------------
Macros to parse int numbers (as for residue sequence numbers) and
8.3f float numbers (atom coordinates). 

p point to the end of number field, num and k are int parameter, modified
by both macros. the result is passed to the first parameter.
-----------------------------------------------------------------------------*/

#define	JZ_PROTEIN_PDB_SCAN_INT4(n,p,num)				\
	do\
	{\
		for(n=0,num=1;*p>='0'&&*p<='9';p--,num*=10)n+=(*p-'0')*num;\
		if(*p=='-')n=-n;\
	}while(0)
	
#define	JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k)				\
	do\
	{\
		num=((int)p[-2]*10+(int)p[-1])*10+p[0]-111*(int)'0';\
		p-=4;\
		for(k=1000;*p>='0'&&*p<='9';p--,k*=10)num+=(*p-'0')*k;\
		if(*p=='-')num=-num;\
		d=num*0.001;\
	}while(0)
	
/*-----------------------------------------------------------------------------
load a protein given a pdb file name. 
-----------------------------------------------------------------------------*/

int	jz_protein_load(jz_protein *protein,char *fn,int mode)
{
	char			**LINES,**lines,*line,*p;
	int			nline=1024;
	jz_protein_atom		*atoms=NULL,*patom;
	jz_protein_residue	*residues=NULL,*pres;
	int			nres=0,natom=0;
	int			size_res=0;
	int			res_num=-9999999,insert_code='*';
	int			k,num,n;
	float			d;
	
	/* Load the lines from file */
	lines=LINES=jz_file_read_lines(fn,&nline,0);
	if(!lines)return 1;
	
	/* Seek the first ATOM line */
	for(;;lines++)
	{
		line=*lines;
		if(!line)
		{
			jz_str_lines_free(LINES);
			return 2;
		}		
		if(JZ_PROTEIN_PDB_ATOM(line))break;
	}

	/* Initialize the residues array. Residues array is used always */
	JZ_ARRAY_INIT(residues,1024);
	size_res=1024;
	
	/* Load atoms */
	if(mode&JZ_PROTEIN_ATOM)
	{

	  JZ_ARRAY_INIT(atoms,nline+nline+32);
	
		/* process ATOM lines */
		for(;(line=*lines);lines++)
		{
			/* stop on TER line or end of file */
			/*if(JZ_PROTEIN_PDB_TER(line))break;*/
			if(JZ_PROTEIN_PDB_ENDMDL(line))break;
			
			/* skip HETATM lines */
			if(!JZ_PROTEIN_PDB_ATOM(line))continue;
			
			/* skip alternative locations */
			if(line[16]!=' ')continue;
		
			/* get num=residue sequence number */
			p=line+25;
			JZ_PROTEIN_PDB_SCAN_INT4(num,p,k);
			
			/* if start of a new residue */
			if(num!=res_num||line[26]!=insert_code)
			{
				JZ_ARRAY_NEW(residues,pres,size_res,nres,2);
				pres->name=jz_amino_decode_name(line+17);
				JZ_CASTR(pres->atoms,int)=natom;
				pres->num=((uint16_t)num)&8191;
				if(line[26]!=' ')pres->num
					|=((uint16_t)(line[26]-'A'+1)<<13);
				res_num=num;
				insert_code=line[26];
			}
			
			/* add the atom to the residue */
			patom=atoms+natom++;
			/* to leave space for terminal atoms */
			patom->name[0]=line[13];
			patom->name[1]=line[14];
			p=(char*)line+37;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			patom->coord[0]=d;
			p=(char*)line+45;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			patom->coord[1]=d;
			p=(char*)line+53;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			patom->coord[2]=d;

			strncpy(patom->first, &line[0], 30);
			patom->first[31] = '\0';

			strncpy(patom->last, &line[54], 26);
			patom->last[27] = '\0';
		}

		/* free the lines here */
		jz_str_lines_free(LINES);

		/* resize atoms */
		if(natom>0)JZ_ARRAY_RESIZE(atoms,natom,patom);
		else
		{
			JZ_ARRAY_FREE(atoms);
			atoms=NULL;
		}

		/* link residues to atoms */
		for(k=0;k<nres;k++)
		{
			residues[k].atoms=atoms
				+JZ_CASTR(residues[k].atoms,int);
			JZ_VEC3_CPY(residues[k].center,
				residues[k].atoms[1].coord);
		}
	}

	/* otherwise load residues only, using ATOM CA lines */
	else
	{
		for(;(line=*lines);lines++)
		{
			/*if(JZ_PROTEIN_PDB_TER(line))break;*/
			if(JZ_PROTEIN_PDB_ENDMDL(line))break;
			if(!JZ_PROTEIN_PDB_ATOM(line))continue;
			if(!JZ_PROTEIN_PDB_CA(line))continue;
			if(line[16]!=' ')continue;
			JZ_ARRAY_NEW(residues,pres,size_res,nres,2);
			p=(char*)line+37;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			pres->center[0]=d;
			p=(char*)line+45;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			pres->center[1]=d;
			p=(char*)line+53;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			pres->center[2]=d;
			pres->atoms=NULL;
			pres->name=jz_amino_decode_name(line+17);
			p=line+25;
			JZ_PROTEIN_PDB_SCAN_INT4(num,p,k);
			pres->num=((uint16_t)num)&8191;
			if(line[26]!=' ')
				pres->num|=((uint16_t)(line[26]-'A'+1)<<13);
		}
		
		/* free the lines here */
		jz_str_lines_free(LINES);
	}
	
	/* now resize residues */
	if(nres&&residues)JZ_ARRAY_RESIZE(residues,nres,pres);
	else
	{
		JZ_ARRAY_FREE(residues);
		residues=NULL;
		nres=0;
	}
	
	/* fix natom in residues */
	if(residues)
	{
		if(atoms)
		{
			n=nres-1;pres=residues;
			for(k=0;k<n;k++,pres++)
				pres->natom=pres[1].atoms-pres->atoms;
			if(nres>0)pres->natom=
				natom-(pres[-1].atoms-atoms)-pres[-1].natom;
			else pres->natom=natom;
		}
		else
		{
			for(k=0;k<nres;k++)
				residues[k].natom=0;
		}
	}

	protein->atoms=atoms;
	protein->residues=residues;
	protein->natom=natom;
	protein->nres=nres;

	return 0;
}


/*-----------------------------------------------------------------------------
load a protein given a pdb file name. 
-----------------------------------------------------------------------------*/

int	bp_protein_load(jz_protein *protein,const char *fn,int mode, char chn, int low_res, int high_res)
{
	char			**LINES,**lines,*line,*p;
	int			nline=1024;
	jz_protein_atom		*atoms=NULL,*patom;
	jz_protein_residue	*residues=NULL,*pres;
	int			nres=0,natom=0;
	int			size_res=0;
	int			res_num=-9999999,insert_code='*';
	int			k,num,n;
	float			d;
	
	/* Load the lines from file */
	lines=LINES=jz_file_read_lines(fn,&nline,0);
	if(!lines)return 1;
	
	/* Seek the first ATOM line */
	for(;;lines++)
	{
		line=*lines;
		if(!line)
		{
			jz_str_lines_free(LINES);
			return 2;
		}		
		if(JZ_PROTEIN_PDB_ATOM(line))break;
	}

	/* Initialize the residues array. Residues array is used always */
	JZ_ARRAY_INIT(residues,1024);
	size_res=1024;
	
	/* Load atoms */
	if(mode&JZ_PROTEIN_ATOM)
	{

	  JZ_ARRAY_INIT(atoms,nline+nline+32);
	
		/* process ATOM lines */
		for(;(line=*lines);lines++)
		{
			/* stop on TER line or end of file */
			/*if(JZ_PROTEIN_PDB_TER(line))break;*/
			if(JZ_PROTEIN_PDB_ENDMDL(line))break;
			
			/* skip HETATM lines */
			if(!JZ_PROTEIN_PDB_ATOM(line))continue;
			
			/* skip alternative locations */
			if(line[16]!=' ')continue;
		
			/* skip if chain does not match */
			if(line[21]!=chn)continue;

			/* get num=residue sequence number */
			p=line+25;
			JZ_PROTEIN_PDB_SCAN_INT4(num,p,k);
			
			/* skip if the residue number is out of the specified range */
			if ((num < low_res) || (num > high_res)) continue;

			/* if start of a new residue */
			if(num!=res_num||line[26]!=insert_code)
			{
				JZ_ARRAY_NEW(residues,pres,size_res,nres,2);
				pres->name=jz_amino_decode_name(line+17);
				JZ_CASTR(pres->atoms,int)=natom;
				pres->num=((uint16_t)num)&8191;
				if(line[26]!=' ')pres->num
					|=((uint16_t)(line[26]-'A'+1)<<13);
				res_num=num;
				insert_code=line[26];
			}
			
			/* add the atom to the residue */
			patom=atoms+natom++;
			/* to leave space for terminal atoms */
			patom->name[0]=line[13];
			patom->name[1]=line[14];
			p=(char*)line+37;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			patom->coord[0]=d;
			p=(char*)line+45;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			patom->coord[1]=d;
			p=(char*)line+53;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			patom->coord[2]=d;

			strncpy(patom->first, &line[0], 30);
			patom->first[31] = '\0';

			strncpy(patom->last, &line[54], 26);
			patom->last[27] = '\0';
		}

		/* free the lines here */
		jz_str_lines_free(LINES);

		/* resize atoms */
		if(natom>0)JZ_ARRAY_RESIZE(atoms,natom,patom);
		else
		{
			JZ_ARRAY_FREE(atoms);
			atoms=NULL;
		}

		/* link residues to atoms */
		for(k=0;k<nres;k++)
		{
			residues[k].atoms=atoms
				+JZ_CASTR(residues[k].atoms,int);
			JZ_VEC3_CPY(residues[k].center,
				residues[k].atoms[1].coord);
		}
	}

	/* otherwise load residues only, using ATOM CA lines */
	else
	{
		for(;(line=*lines);lines++)
		{
			/*if(JZ_PROTEIN_PDB_TER(line))break;*/
			if(JZ_PROTEIN_PDB_ENDMDL(line))break;
			if(!JZ_PROTEIN_PDB_ATOM(line))continue;
			if(!JZ_PROTEIN_PDB_CA(line))continue;
			if(line[16]!=' ')continue;

			/* skip if chain does not match */
			if(line[21]!=chn)continue;
			
			/* get num=residue sequence number */
			p=line+25;
			JZ_PROTEIN_PDB_SCAN_INT4(num,p,k);
			
			/* skip if the residue number is out of the specified range */
			if ((num < low_res) || (num > high_res)) continue;

			
			JZ_ARRAY_NEW(residues,pres,size_res,nres,2);
			p=(char*)line+37;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			pres->center[0]=d;
			p=(char*)line+45;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			pres->center[1]=d;
			p=(char*)line+53;
			JZ_PROTEIN_PDB_SCAN_FLOAT83(d,p,num,k);
			pres->center[2]=d;
			pres->atoms=NULL;
			pres->name=jz_amino_decode_name(line+17);
			p=line+25;
			JZ_PROTEIN_PDB_SCAN_INT4(num,p,k);
			pres->num=((uint16_t)num)&8191;
			if(line[26]!=' ')
				pres->num|=((uint16_t)(line[26]-'A'+1)<<13);
		}
		
		/* free the lines here */
		jz_str_lines_free(LINES);
	}
	
	/* now resize residues */
	if(nres&&residues)JZ_ARRAY_RESIZE(residues,nres,pres);
	else
	{
		JZ_ARRAY_FREE(residues);
		residues=NULL;
		nres=0;
	}
	
	/* fix natom in residues */
	if(residues)
	{
		if(atoms)
		{
			n=nres-1;pres=residues;
			for(k=0;k<n;k++,pres++)
				pres->natom=pres[1].atoms-pres->atoms;
			if(nres>0)pres->natom=
				natom-(pres[-1].atoms-atoms)-pres[-1].natom;
			else pres->natom=natom;
		}
		else
		{
			for(k=0;k<nres;k++)
				residues[k].natom=0;
		}
	}

	protein->atoms=atoms;
	protein->residues=residues;
	protein->natom=natom;
	protein->nres=nres;

	return 0;
}

/*-----------------------------------------------------------------------------
fill a protein instance with all zeros
-----------------------------------------------------------------------------*/

void	jz_protein_null(jz_protein *pt)
{
	memset(pt,0,sizeof(jz_protein));
}


/*-----------------------------------------------------------------------------
jz_protein_pseudo_atom_lines	print pseudo ATOM lines. CA atoms only.
-----------------------------------------------------------------------------*/

void	jz_protein_pseudo_atom_lines(jz_protein *pt,FILE *fout,int chain,
	int start)
{
	int		i,k,m,n,count,seq,ch;
	jz_protein_residue	*reses;
	jz_protein_atom		*atoms;
	
	if(!pt)return;
	atoms=pt->atoms;
	if(pt->natom<=0)atoms=NULL;
	
	reses=pt->residues;
	n=pt->nres;
	
	if(atoms)
	{
	  for(count=k=0;k<n;k++)
		{	
			atoms=reses[k].atoms;
			m=reses[k].natom;
			if(!atoms||m<=0)continue;
			seq=reses[k].num&4095;
			ch=seq>>12;
			if(ch)ch='A'-1+ch;
			else ch=' ';
			for(i=0;i<m;i++)
			{
			  fprintf(fout,"ATOM  %5d  %c%c  %s %c%4d%c   %8.3f%8.3f%8.3f\n",
				  count,atoms[i].name[0],atoms[i].name[1],
				  jz_amino_name[reses[k].name],chain,seq,ch,atoms[i].coord[0],
				  atoms[i].coord[1],atoms[i].coord[2]);
			  count++;
			}
		}	
		
		return;	
	}

	for(k=0;k<n;k++)
	{
		seq=reses[k].num&4095;
		ch=seq>>12;
		if(ch)ch='A'-1+ch;
		else ch=' ';
		fprintf(fout,"ATOM  %5d  CA  %s %c%4d%c   %8.3f%8.3f%8.3f\n",
			k+start,jz_amino_name[reses[k].name],chain,seq,ch,
			reses[k].center[0],reses[k].center[1],
			reses[k].center[2]);

	}
}

void	jz_protein_output_atom_lines(jz_protein *pt,FILE *fout,int chain,
	int start)
{
	int		i,k,m,n,count,seq,ch;
	jz_protein_residue	*reses;
	jz_protein_atom		*atoms;
	
	if(!pt)return;
	atoms=pt->atoms;
	if(pt->natom<=0)atoms=NULL;
	
	reses=pt->residues;
	n=pt->nres;
	
	if(atoms)
	{
	  for(count=k=0;k<n;k++)
		{	
			atoms=reses[k].atoms;
			m=reses[k].natom;
			if(!atoms||m<=0)continue;
	
			for(i=0;i<m;i++)
			{
			  fprintf(fout,"%s%8.3f%8.3f%8.3f%s\n",
				  atoms[i].first,atoms[i].coord[0],atoms[i].coord[1],
				  atoms[i].coord[2],atoms[i].last);
			  count++;
			}
		}	
		
		return;	
	}

	for(k=0;k<n;k++)
	{
		seq=reses[k].num&4095;
		ch=seq>>12;
		if(ch)ch='A'-1+ch;
		else ch=' ';
		fprintf(fout,"ATOM  %5d  CA  %s %c%4d%c   %8.3f%8.3f%8.3f\n",
			k+start,jz_amino_name[reses[k].name],chain,seq,ch,
			reses[k].center[0],reses[k].center[1],
			reses[k].center[2]);

	}
}


/*-----------------------------------------------------------------------------
translate the structure, dst has been cloned from src. the behavior is
determined by the option parameter. 
-----------------------------------------------------------------------------*/

void	jz_protein_transform(jz_protein *dst,jz_protein *src,float *ori,
	float *rot,int option)
{
	float		*p,*q,foo[3];
	int		k;

	/* transform all residues */
	if(option&4&&dst->residues)
	{		
		p=dst->residues->center;
		q=src->residues->center;
		for(k=dst->nres;k;k--)
		{
			JZ_VEC3_SUB(foo,q,ori);
			JZ_VEC3_ROT(p,foo,rot);
			p=(float*)((char*)p+sizeof(jz_protein_residue));
			q=(float*)((char*)q+sizeof(jz_protein_residue));
		}
	}

	/* transform all atoms */
	if(option&8&&dst->atoms)
	{
		p=dst->atoms->coord;
		q=src->atoms->coord;
		for(k=dst->natom;k;k--)
		{
			JZ_VEC3_SUB(foo,q,ori);
			JZ_VEC3_ROT(p,foo,rot);
			p=(float*)((char*)p+sizeof(jz_protein_atom));
			q=(float*)((char*)q+sizeof(jz_protein_atom));
		}
	}
}



