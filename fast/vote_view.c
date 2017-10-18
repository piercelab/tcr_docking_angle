/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <stdio.h>

#include "basic.h"
#include "chore.h"
#include "jz_zim.h"

/*-----------------------------------------------------------------------------
 * print integer matrix to be viewed in matlab with load_image.m
-----------------------------------------------------------------------------*/

int	vote_print_matrix(int *MAT,int M,int N,int LD,char *fn)
{
	FILE		*f;
	int		*P;
	int		i,j;

	f=fopen(fn,"wt");
	if(!f)return 1;
	fprintf(f,"ZIM ASCII PF8BIT\n%d\n%d\n",M,N);

	for(P=MAT,i=0;i<M;i++,P+=LD)
	for(j=0;j<N;j++)
		fprintf(f,"%d\n",P[j]);
	fclose(f);
	return 0;
}

/*-----------------------------------------------------------------------------
* View matrix
-----------------------------------------------------------------------------*/

void	vote_view_matrix(int *Mat,char *msg,int M,int N,int option)
{
	char	str[256];
	int	k,n=M*N;
	int	*T;
	int	m1,m2;
	float	a,b;

	fprintf(stderr,"%s:",msg);
	fgets(str,256,stdin);
	jz_str_trim(str);
	JZ_ARRAY_INIT(T,M*N);
	JZ_ARRAY_COPYF(T,Mat,M*N);
	
	/* 0: directly print the matrix */
	
	/* 1: remove negatives */
	if(option==1)
	{
		for(k=0;k<n;k++)
			if(T[k]<0)T[k]=0;
	}
	
	/* reverse black-white and normalize */
	m1=m2=T[0];
	for(k=0;k<n;k++)
	{
		if(T[k]<m1)m1=T[k];
		if(T[k]>m2)m2=T[k];
	}

	a=255.0/(m1-m2);
	b=-a*m2;
	
	for(k=0;k<n;k++)
		T[k]=(int)(a*T[k]+b);
	
	/* pause to allow visualization of the file */
	if(vote_print_matrix(T,M,N,N,str))
		fprintf(stderr,"FAILED, press any key to continue...\n");
	else fprintf(stderr,"Saved, presee any key to continue...\n");	
	JZ_MATRIX_FREE(T);
	getchar();
}

/*-----------------------------------------------------------------------------
* View matrix
-----------------------------------------------------------------------------*/

void	vote_view_allow(char *Mat,char *msg,int M,int N)
{
	char	str[256];
	int	k,n=M*N;
	int	*T;

	fprintf(stderr,"%s:",msg);
	fgets(str,256,stdin);
	jz_str_trim(str);
	JZ_ARRAY_INIT(T,M*N);
	
	for(k=0;k<n;k++)
		if(Mat[k]==0)T[k]=1;
		else T[k]=0;
	
	/* pause to allow visualization of the file */
	if(vote_print_matrix(T,M,N,N,str))
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


