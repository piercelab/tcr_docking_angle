/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
fast.c		pairwise alignment main program
=============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "basic.h"
#include "chore.h"
#include "misc.h"
#include "jz_amino.h"
#include "jz_protein.h"
#include "Alignment.h"
#include "rasmol.h"
#include "vote.h"

extern	int	refine_vote_count;

int	main(int argc,char **argv)
{
	char		*usage="\n\
FAST - FAST Alignment and Search Tool, is a protein structure\n\
alignment algorithm developed by Jianhua Zhu at the Bioinformatics\n\
Program of Boston University. Please visit the FAST homepage\n\
at http://biowulf.bu.edu/FAST/.\n\
\n\
USAGE: $ fast protein1 protein2 [-r rasmol-script]\n\
Try 'fast --help' for more information.\n\
\n";

	char		*help="\n\
FAST - FAST Alignment and Search Tool, is a protein structure\n\
alignment algorithm developed by Jianhua Zhu at the Bioinformatics\n\
Program of Boston University. Please visit the FAST homepage\n\
at http://biowulf.bu.edu/FAST/.\n\
\n\
USAGE: $ fast protein1 protein2 [-r rasmol-script]\n\
\n\
Description:\n\
%1		fn	protein #1\n\
%2		fn	protein #2\n\
-o|--out	fn	output file name [stdout]\n\
-r|--rasmol	fn	generate a rasmol script for the alignment\n\
-f|--fit        fn      output fitted protein2 coords to file\n\
-p|--path	fn	path to protein structures\n\
-h|--help		display this message\n\
\n\
To visualize alignment under unix/linux:\n\
	$ rasmol -script <rasmol-script>\n\
\n\
Email to jianhua@bu.edu for bug report.\n\
\n";
	
	/* variables to hold line arguments */

	char		*protein1_name=NULL;
	char		*protein2_name=NULL;
	char		*out_fn=NULL;
	
	char            *fit_fn=NULL;
	char		*rasmol_fn=NULL;
	char		*path=NULL;
	
	/* command line description */
	jz_cmd		cmds[]=
	{
			{2,JZ_TSTR,&protein1_name,NULL},
			{2,JZ_TSTR,&protein2_name,NULL},
			{3,JZ_TSTR,&out_fn,"o|out"},
			{3,JZ_TSTR,&fit_fn,"f|fit"},
			{3,JZ_TSTR,&rasmol_fn,"r|R|ras|rasmol"},
			{3,JZ_TSTR,&path,"p|path"},
			{1,JZ_TMESSAGE,&help,"h|help"},
	};
	
	/* working variables */
	int		ret=0;
	FILE		*fout=NULL;
	FILE            *fit_file=NULL;
	FILE		*rasmol_file=NULL;
	jz_protein	protein1,*p1;
	jz_protein	protein2,*p2;
	Alignment	alignment;

	char		*error_messages[]=
	{
			"OK",
			"Bad protein name(s)",
			"Loading structure #1 failed",
			"Loading structure #2 failed",
	};

	char		full_fn1[256],full_fn2[256];
	int		path_len=0;
	int		err;
	char		path_ch='/';
	
	if(argc==1)
	{
		printf(usage);
		return 1;
	}

	if(JZ_CMD_PARSE(cmds))
	{
		jz_cmd_message(stdout);
		return 1;
	}

	/* open output file */
	if(!JZ_CMD_USED(cmds[2]))fout=stdout;
	else if(!(fout=fopen(out_fn,"wt")))
	{
		printf("fast: cannot open output file: %s\n",out_fn);
		return 2;
	}
	
	/* initialize the alignment space */
	alignment_create(&alignment);

	/* prepare full file names */
	if(path) /* --path */
	{
		if(strchr(path,'\\'))path_ch='\\';
		else if(strchr(path,'/'))path_ch='/';
		
		path_len=strlen(path);
		strcpy(full_fn1,path);
		strcpy(full_fn2,path);
		if(path[path_len-1]!=path_ch)
		{
			full_fn1[path_len]=path_ch;
			full_fn1[path_len+1]=0;
			full_fn2[path_len++]=path_ch;
			full_fn2[path_len]=0;
		}
	}
	
	/* open file to save rasmol script */
	if(rasmol_fn&&!(rasmol_file=fopen(rasmol_fn,"wt")))
	{
		printf("fast: cannot create rasmol file: %s\n",rasmol_fn);
		ret=6;goto __EXIT;
	}
	
	/* open file to save fitted coords */
	if(fit_fn&&!(fit_file=fopen(fit_fn,"wt")))
	{
		printf("fast: cannot create PDB file: %s\n",fit_fn);
		ret=6;goto __EXIT;
	}
	

	/* initialize protein strucutres */
	jz_protein_null(&protein1);
	jz_protein_null(&protein2);

	/* prepare protein1_name and protein2_name */
	if(!protein2_name)protein2_name=protein1_name;
	if(!protein1_name)JZ_ERROR_GOTO(err=1,__DONE);
	if(strlen(protein1_name)>128||strlen(protein2_name)>128)
		JZ_ERROR_GOTO(err=1,__DONE);

	err=0;
	
	/* display the pair of proteins */

	fprintf(fout,"FAST ALIGNMENT: %s %s",protein1_name,protein2_name);

	/* determine the protein structures p1 and p2 */
	if(strchr(protein1_name,'/')||strchr(protein1_name,'\\'))
	{
		if(jz_protein_load(&protein1,protein1_name,255))p1=NULL;
		else p1=&protein1;
	}
	else if(path)
	{
		strcpy(full_fn1+path_len,protein1_name);
		if(jz_protein_load(&protein1,full_fn1,255))p1=NULL;
		else p1=&protein1;
	}
	else
	{
		if(jz_protein_load(&protein1,protein1_name,255))p1=NULL;
		else p1=&protein1;
	}
	if(!p1)JZ_ERROR_GOTO(err=2,__DONE);
	
	if(strchr(protein2_name,'/')||strchr(protein2_name,'\\'))
	{
		if(jz_protein_load(&protein2,protein2_name,255))p2=NULL;
		else p2=&protein2;
	}
	else if(path)
	{
		strcpy(full_fn2+path_len,protein2_name);
		if(jz_protein_load(&protein2,full_fn2,255))p2=NULL;
		else p2=&protein2;
	}
	else
	{
		if(jz_protein_load(&protein2,protein2_name,255))p2=NULL;
		else p2=&protein2;
	}
	if(!p2)JZ_ERROR_GOTO(err=3,__DONE);

	/* the alignment is done here */
	vote_pairwise(p1,p2,&alignment);

	/* display residue alignment */
	/*fprintf(fout,"\nRES_ALIGN: ");*/
	fprintf(fout,"\n");
	alignment_display_info(&alignment,p1,p2,fout);
	alignment_display_rpl(&alignment,p1,p2,fout,60);
	alignment_display_simple(&alignment,p1,p2,fout);
	fflush(fout);
	
	/* save rasmol script */
	if(rasmol_file)
	{	
		vote_rasmol_alignment(&alignment,p1,p2,rasmol_file);
		fclose(rasmol_file);
		rasmol_file=NULL;
	}

	/* save fitted protein structure */
	if(fit_file)
	{	
		vote_output_coords(&alignment,p1,p2,fit_file);
		fclose(fit_file);
		fit_file=NULL;
	}

__DONE:

	if(err)
	{	
		fprintf(fout,"\tError: %s\n\n",error_messages[err]);
	}
	else	fprintf(fout,"\n\n");

__EXIT:

	jz_protein_clear(&protein1);
	jz_protein_clear(&protein2);
	fclose(fout);
	return ret;
}


