#ifndef	__JZ_ZIM_H
#define	__JZ_ZIM_H

#ifdef	__cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
The gray color system
-----------------------------------------------------------------------------*/

typedef	unsigned char		jz_color_gray;

/*-----------------------------------------------------------------------------
The rgb color system
-----------------------------------------------------------------------------*/

typedef	jz_color_gray		jz_color_rgb[3];

/*-----------------------------------------------------------------------------
The hsv color system
-----------------------------------------------------------------------------*/

typedef	float			jz_color_hsv[3];

extern	jz_color_rgb		jz_color_red,jz_color_green,jz_color_blue,
				jz_color_white,jz_color_cgray;

#define	JZ_ZIM_GRAY		0
#define	JZ_ZIM_RGB		1

/*-----------------------------------------------------------------------------
Save a 2 D byte array, given an FILE pointer.
-----------------------------------------------------------------------------*/

int	jz_zim_gray_fwrite(jz_color_gray *im,int M,int N,int L,FILE *fout);

/*-----------------------------------------------------------------------------
Save 2D image, given a file name
-----------------------------------------------------------------------------*/

int	jz_zim_gray_save(jz_color_gray *im,int M,int N,int L,char *fn);

/*-----------------------------------------------------------------------------
Save a color picture, as a jz_color_rgb array.
-----------------------------------------------------------------------------*/

int	jz_zim_rgb_fwrite(jz_color_rgb *im,int M,int N,int L,FILE *fout);

/*-----------------------------------------------------------------------------
Save a color picture, as a jz_color_rgb array.
-----------------------------------------------------------------------------*/

int	jz_zim_rgb_save(jz_color_rgb *im,int M,int N,int L,char *fn);


/*-----------------------------------------------------------------------------
Convert an integer matrix into an gray picture
-----------------------------------------------------------------------------*/

jz_color_gray*
jz_zim_int2gray(int *mat,int M,int N,int L);

#ifdef	__cplusplus
}
#endif

#endif


