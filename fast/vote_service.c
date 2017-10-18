/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>
#include <limits.h>

#include "vote.h"
#include "basic.h"
#include "chore.h"
#include "jz_zim.h"

/*-----------------------------------------------------------------------------
 * print integer matrix to be viewed in matlab
-----------------------------------------------------------------------------*/

int	vote_save_matrix(int *MAT,int M,int N,int LD,char *fn) 
{
	jz_color_gray	*img;
	int		ret;

	img=jz_zim_int2gray(MAT,M,N,LD);
	ret=jz_zim_gray_save(img,M,N,N,fn);
	JZ_ARRAY_FREE(img);
	return ret;
}

/*-----------------------------------------------------------------------------
 * print integer matrix to be viewed in matlab with load_image.m
-----------------------------------------------------------------------------*/

void	vote_print_matrix(int *MAT,int M,int N,int LD,char *fn)
{
	FILE		*f;
	int		*P;
	int		i,j;

	f=fopen(fn,"wt");
	fprintf(f,"%d\n%d\n",M,N);

	for(P=MAT,i=0;i<M;i++,P+=LD)
	for(j=0;j<N;j++)
		fprintf(f,"%d\n%d\n%d\n",P[j],P[j],P[j]);
	fclose(f);
}

/*-----------------------------------------------------------------------------
* View matrix
-----------------------------------------------------------------------------*/

void	vote_view_matrix(int *Mat,char *msg,int M,int N,int option)
{
	char	str[256];
	int	k,cutoff;
	int	*p,*T;

	fprintf(stderr,"%s:",msg);
	fgets(str,256,stdin);
	jz_str_trim(str);
	JZ_ARRAY_INIT(T,M*N);
	if(option>=0&&option<4)JZ_ARRAY_COPYF(T,Mat,M*N);
	
	/* 0: directly print the matrix */
	if(option==0);

	/* 1: make it binary */
	else if(option==1)
	{
		for(k=M*N,p=T;k;k--,p++)
			if(*p)*p=1;else *p=0;
	}

	/* 2: remove negatives for LGC */
	else if(option==2)
	{
		for(k=M*N,p=T;k;k--,p++)
			if(*p<0)*p=0;
	}
	
	/* 3: filter by median */
	else if(option==3)
	{
		int 	*buf;
		int	*q;

		JZ_ARRAY_INIT(buf,M*N);
		JZ_ARRAY_COPYF(buf,Mat,M*N);

		q=buf;	
		for(k=M*N,p=T;k;k--,p++)
			if(*p>0)*q++=*p;
		k=q-buf;
		cutoff=jz_select_i(buf,k,k/2);
		for(k=M*N,p=T;k;k--,p++)
			if(*p<cutoff)*p=cutoff;
		JZ_ARRAY_FREE(buf);
	}

	/* pause to allow visualization of the file */
	if(vote_save_matrix(T,M,N,N,str))
		fprintf(stderr,"FAILED, press any key to continue...\n");
	else fprintf(stderr,"Saved, presee any key to continue...\n");	
	JZ_MATRIX_FREE(T);
	getchar();
}

/*-----------------------------------------------------------------------------
Print an alignment with given color
-----------------------------------------------------------------------------*/

void	vote_view_align(jz_pair_int *pairs,char *prompt,int M,int N,
	int len,jz_color_rgb color)
{
	char		str[256];
	int		i,j,k;
	jz_color_gray	*p;
	jz_color_rgb	**T;

	fprintf(stderr,"%s:",prompt);
	fgets(str,256,stdin);
	jz_str_trim(str);
	JZ_MATRIX_INIT(T,M,N);
	for(k=M*N*3,p=&T[0][0][0];k;k--)*p++=0;
	
	for(k=0;k<len;k++)
	{
		i=pairs[k][0];
		j=pairs[k][1];
		p=T[i][j];
		p[0]=color[0];
		p[1]=color[1];
		p[2]=color[2];
	}
	jz_zim_rgb_save(T[0],M,N,N,str);
	JZ_MATRIX_FREE(T);

	/* pause to allow visualization of the file */
	fprintf(stderr,"Saved, presee any key to continue...\n");	
	getchar();
}


