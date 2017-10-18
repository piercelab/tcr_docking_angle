#include "jz_zim.h"
#include <stdio.h>

jz_color_rgb	jz_color_red={255,0,0};
jz_color_rgb	jz_color_green={0,255,0};
jz_color_rgb	jz_color_blue={0,0,255};
jz_color_rgb	jz_color_cgray={128,128,128};
jz_color_rgb	jz_color_white={255,255,255};

/*-----------------------------------------------------------------------------
Save a 2 D byte array, given an FILE pointer.
-----------------------------------------------------------------------------*/

int	jz_zim_gray_fwrite(jz_color_gray *im,int M,int N,int L,FILE *fout)
{
	int32_t		i;

	if(!im)return -1;
	if(!fout)return -2;
	i=0;
	if(fwrite(&i,4,1,fout)!=1)return -2;
	i=M;
	if(fwrite(&i,4,1,fout)!=1)return -2;
	i=N;
	if(fwrite(&i,4,1,fout)!=1)return -2;
	if(N==L)
	{
		if(fwrite(im,M*N*sizeof(jz_color_gray),1,fout)!=1)
			return -2;
	}
	else
	{
		for(;M;M--,im+=L)
			if(fwrite(im,N*sizeof(jz_color_gray),1,fout)!=1)
				return -2;
	}

	return 0;
}

/*-----------------------------------------------------------------------------
Save 2D image, given a file name
-----------------------------------------------------------------------------*/

int	jz_zim_gray_save(jz_color_gray *im,int M,int N,int L,char *fn)
{
	FILE		*fout;
	int		ret;

	if(!im)return -1;
	fout=fopen(fn,"w");
	if(!fout)return -2;
	ret=jz_zim_gray_fwrite(im,M,N,L,fout);
	if(ret)ret-=2;
	fclose(fout);
	return ret;
}

/*-----------------------------------------------------------------------------
Save a color picture, as a jz_color_rgb array.
-----------------------------------------------------------------------------*/

int	jz_zim_rgb_fwrite(jz_color_rgb *im,int M,int N,int L,FILE *fout)
{
	int32_t		i;

	if(!im)return -1;
	if(!fout)return -2;
	i=1;
	if(fwrite(&i,4,1,fout)!=1)return -2;
	i=M;
	if(fwrite(&i,4,1,fout)!=1)return -2;
	i=N;
	if(fwrite(&i,4,1,fout)!=1)return -2;
	if(N==L)
	{
		if(fwrite(im,M*N*sizeof(jz_color_rgb),1,fout)!=1)
			return -2;
	}
	else
	{
		for(;M;M--,im+=L)
			if(fwrite(im,N*sizeof(jz_color_rgb),1,fout)!=1)
				return -2;
	}

	return 0;
}

/*-----------------------------------------------------------------------------
Save a color picture, as a jz_color_rgb array.
-----------------------------------------------------------------------------*/

int	jz_zim_rgb_save(jz_color_rgb *im,int M,int N,int L,char *fn)
{
	FILE		*fout;
	int		ret;

	if(!im)return -1;
	fout=fopen(fn,"w");
	if(!fout)return -2;
	ret=jz_zim_rgb_fwrite(im,M,N,L,fout);
	if(ret)ret-=2;
	fclose(fout);
	return ret;
}

/*-----------------------------------------------------------------------------
convert an integer matrix into a gray array
-----------------------------------------------------------------------------*/

jz_color_gray*
jz_zim_int2gray(int *mat,int M,int N,int L)
{
	jz_color_gray	*ret;
	int		i,k;
	int		mx,mn;
	float		factor;

	JZ_ARRAY_INIT(ret,M*N);
	if(!ret)return NULL;

	mx=mn=mat[0];
	if(L<N)L=N;
	if(N==L)
	{
		for(k=M*N-1;k;k--)
			if(mat[k]<mn)mn=mat[k];
			else if(mat[k]>mx)mx=mat[k];
	}
	else
	{
		for(--mat,i=M;i;i--,mat+=L)
		{
			for(k=N;k;k--)
				if(mat[k]<mn)mn=mat[k];
				else if(mat[k]>mx)mx=mat[k];
		}
		mat=mat-M*L+1;
	}
	
	factor=255.0/(mx-mn);
	
	if(N==L)
	{
		for(--ret,--mat,k=M*N;k;k--)
			ret[k]=(jz_color_gray)((mat[k]-mn)*factor);
		return ret+1;
	}

	else
	{
		for(--mat,--ret,i=M;i;i--,ret+=N,mat+=L)
		{
			for(k=N;k;k--)
			{
				ret[k]=(jz_color_gray)((mat[k]-mn)*factor);
			}
		}
		return ret-M*N;
	}
}


