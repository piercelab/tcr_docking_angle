/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "vote.h"

/*-----------------------------------------------------------------------------
Given an alignment, get the rotation matrix, using rms_fit
-----------------------------------------------------------------------------*/

int	vote_rms_fit(jz_pair_int *pairs,int N,jz_protein *pt1,jz_protein *pt2,
	float *C,float *R)
{
	float		*P,*Q,*pp,*qq,*p,*q,ss;
	int		k;
	jz_protein_residue	*R1,*R2;
	
	JZ_ARRAY_INIT(P,N*6);
	Q=P+N*3;
	R1=pt1->residues;
	R2=pt2->residues;
	pp=P,qq=Q;
	for(k=0;k<N;k++)
	{
		p=R1[pairs[k][0]].center;
		JZ_VEC3_CPY(pp,p);
		q=R2[pairs[k][1]].center;
		JZ_VEC3_CPY(qq,q);
		pp+=3;
		qq+=3;
	}

	ss=jz_rms_fitf(P,Q,N,C,R);
	if(ss<0.0)return ss=-ss;;
	JZ_ARRAY_FREE(P);
	return 0;
}


#include <stdio.h>
#include <math.h>
#include "vote.h"

/*-----------------------------------------------------------------------------
P and Q each is a float[N][3] array. get the center C and and the
rotation matrix of the optimal transformation to fit Q to P. On return,
P is unchanged, Q is modified accoring to the best fitting. C is the
origin of the new coordinate system as in the old one. R is the
rotation matrix. if C is NULL, then assume that origin is correct,
only rotation is performed.

This function was manually translated from Zhiping's MATCH.F. Joe's f2c 
version was read for reference. 
-----------------------------------------------------------------------------*/

#define	RMS_SHIFT	(sizeof(float)*3)

float	jz_rms_fitf(float *P,float *Q,int N,float *C,float *R)
{
	float	*pp,*qq;
	float	QB[8],B[8];
	float	A[16],F[16];
	float	alpha,beta,gamma;
	float	alphsave,betasave,gammsave;
	float	alphanum,betanum,gammanum;
	float	alphaden,betaden,gammaden;
	float	sa,ca,sb,cb,sc,cc;
	float	d,x,y,z,step;
	float	pi,ss,ss0,sumsqs;
	int		i,j,k;
	/*int		count=0;*/

	alpha=beta=gamma=alphsave=betasave=gammsave=0.0;
	if(C)
	{	
		/* COMPUTE THE CENTROIDS AND TRANSLATE BOTH SETS OF POINTS 
		 * TO HAVE CENTROIDS AT (0,0,0). */
	
		pp=P,qq=Q;
		for(i=N;i;i--)
		{
			alpha+=pp[0];
			beta+=pp[1];
			gamma+=pp[2];
			alphsave+=qq[0];
			betasave+=qq[1];
			gammsave+=qq[2];
			pp=(float*)((char*)pp+RMS_SHIFT);
			qq=(float*)((char*)qq+RMS_SHIFT);
		}
		x=1.0/N;
		alpha*=x;
		beta*=x;
		gamma*=x;
		alphsave*=x;
		betasave*=x;
		gammsave*=x;
		C[0]=alpha;
		C[1]=beta;
		C[2]=gamma;
		QB[0]=alphsave;
		QB[1]=betasave;
		QB[2]=gammsave;
	}
	pp=P,qq=Q;
	ss=0.0;
	for(i=N;i;i--)
	{
		pp[0]-=alpha;
		pp[1]-=beta;
		pp[2]-=gamma;
		qq[0]-=alphsave;
		qq[1]-=betasave;
		qq[2]-=gammsave;
		d=pp[0]-qq[0];
		ss+=d*d;
		d=pp[1]-qq[1];
		ss+=d*d;
		d=pp[2]-qq[2];
		ss+=d*d;
		pp=(float*)((char*)pp+RMS_SHIFT);
		qq=(float*)((char*)qq+RMS_SHIFT);
	}
	ss/=N;
	ss0=sqrt(ss);

/*
__RETRY:	
*/

	/* COMPUTE THE MATRIX A = SUMM pi*QI-TRANSPOSE AND THE MAX SUM 
	 * OF SQUARES  PP + QQ */

	pp=P,qq=Q;
	d=0.0;
	A[0]=A[1]=A[2]=A[3]=A[4]=A[5]=A[6]=A[7]=A[8]=d;
	for(i=N;i;i--)
	{
		d=pp[0];
		A[0]+=d*qq[0];
		A[1]+=d*qq[1];
		A[2]+=d*qq[2];
		d=pp[1];
		A[3]+=d*qq[0];
		A[4]+=d*qq[1];
		A[5]+=d*qq[2];
		d=pp[2];
		A[6]+=d*qq[0];
		A[7]+=d*qq[1];
		A[8]+=d*qq[2];
		pp=(float*)((char*)pp+RMS_SHIFT);
		qq=(float*)((char*)qq+RMS_SHIFT);
	}

	pi=3.14159265358979323846;
	sumsqs=0.0;
	alphsave=betasave=gammsave=sumsqs;

	/* SEARCH FOR BETA IN [-pi/2,pi/2], GAMMA IN [0,2*pi), BEST ALPHA
	 * FOR THE LARGEST VALUE OF (P-TRANSPOSE)*R*Q.  THE PROCEDURE COULD
	 * BE ACCELERATED BY COMPUTING AND STORING SINE AND COSINE VALUES. */
	
	step=pi/24.0;
	for(i=0;i<48;i++)
	{
		gamma=i*step;
		sc=sin(gamma);
		cc=cos(gamma);
		for(j=-12;j<=12;j++)
		{
			beta=j*step;
			sb=sin(beta);
			cb=cos(beta);
			alphanum=A[3]*cb-A[4]*sb*sc-A[1]*cc
				-A[5]*sb*cc+A[2]*sc;
			alphaden =A[0]*cb-A[1]*sb*sc+A[4]*cc
				-A[2]*sb*cc-A[5]*sc;

			/* ALPHA IS CHOSEN SO THAT THE FIRST DERIVATIVE 
			 * OF P-TRANSPOSE*R*Q IS ZERO AND THE SECOND 
			 * DERIVATIVE IS NEGATIVE.  THUS P-TRANSPOSE*R*Q
			 * IS MAXIMUM FOR THE VALUES OF GAMMA AND BETA THAT 
			 * ARE BEING TESTED. */

			if(alphaden==0.0)
			{
				if(alphanum>=0.0)alpha=0.5*pi;
				else alpha=-0.5*pi;
			}
			else alpha=atan(alphanum/alphaden);
			if(alphaden<0.0)alpha+=pi;

			/* NOW TEST THE SUM OF SQUARES (PP + QQ 
			 * - 2*P-trans RQ = PP + QQ - 2*SS) 
			 * FOR THIS ROTATION OF Q */

			sa=sin(alpha);
			ca=cos(alpha);

			R[0]=ca*cb;
			R[1]=-ca*sb*sc-sa*cc;
			R[2]=-ca*sb*cc+sa*sc;
			R[3]=sa*cb;
			R[4]=-sa*sb*sc+ca*cc;
			R[5]=-sa*sb*cc-ca*sc;
			R[6]=sb;
			R[7]=cb*sc;
			R[8]=cb*cc;
			
			ss=0.0;
			for(k=0;k<9;k++)ss+=A[k]*R[k];
			if(ss>sumsqs)
			{
				sumsqs=ss;
				alphsave=alpha;
				betasave=beta;
				gammsave=gamma;
			}
		}
	}
	alpha=alphsave;
	beta=betasave;
	gamma=gammsave;
	sa=sin(alpha);
	ca=cos(alpha);
	sb=sin(beta);
	cb=cos(beta);
	sc=sin(gamma);
	cc=cos(gamma);

	/* NOW DO 10 NEWTON'S METHOD ITERATIONS TO IMPROVE THE SOLUTION */

	for(i=10;i;i--)
	{
		alphanum=A[3]*cb-A[4]*sb*sc-A[1]*cc-A[5]*sb*cc+A[2]*sc;
		alphaden=A[0]*cb-A[1]*sb*sc+A[4]*cc-A[2]*sb*cc-A[5]*sc;
		betanum=A[6]-A[1]*ca*sc-A[4]*sa*sc-A[2]*ca*cc-A[5]*sa*cc;
		betaden=A[0]*ca+A[3]*sa+A[7]*sc+A[8]*cc;
		gammanum=-A[1]*ca*sb-A[4]*sa*sb+A[7]*cb+A[2]*sa-A[5]*ca;
		gammaden=-A[1]*sa+A[4]*ca-A[2]*ca*sb-A[5]*sa*sb+A[8]*cb;
		B[0]=-sa*alphaden+ca*alphanum;
		B[1]=-sb*betaden+cb*betanum;
		B[2]=-sc*gammaden+cc*gammanum;
		F[0]=-ca*alphaden-sa*alphanum;
		F[1]=-sa*(-A[0]*sb-A[1]*cb*sc-A[2]*cb*cc)
			+ca*(-A[3]*sb-A[4]*cb*sc-A[5]*cb*cc);
		F[2]=-sa*(-A[1]*sb*cc-A[4]*sc+A[2]*sb*sc-A[5]*cc)
			+ca*(-A[4]*sb*cc+A[1]*sc+A[2]*cc+A[5]*sb*sc);
		F[3]=F[1];
		F[4]=-cb*betaden-sb*betanum;
		F[5]=-sb*(A[7]*cc-A[8]*sc)+cb*(-A[1]*ca*cc
			-A[4]*sa*cc+A[2]*ca*sc+A[5]*sa*sc);
		F[6]=F[2];
		F[7]=F[5];
		F[8]=-cc*gammaden-sc*gammanum;
		
		/* solve equations inlined ==> B */
		x=F[4]*F[8]-F[5]*F[7];
		y=F[2]*F[7]-F[1]*F[8];
		z=F[1]*F[5]-F[2]*F[4];
		d=F[0]*x+F[3]*y+F[6]*z;
		d=1.0/d;
		x=B[0]*x+B[1]*y+B[2]*z;
		y=B[0]*(F[5]*F[6]-F[3]*F[8])+B[1]*(F[0]*F[8]-F[2]*F[6])
			+B[2]*(F[2]*F[3]-F[0]*F[5]);
		z=B[0]*(F[3]*F[7]-F[4]*F[6])+B[1]*(F[1]*F[6]-F[0]*F[7])
			+B[2]*(F[0]*F[4]-F[1]*F[3]);
		B[0]=x*d;
		B[1]=y*d;
		B[2]=z*d;
		
		alpha=alphsave-B[0];
		beta=betasave-B[1];
		gamma=gammsave-B[2];

		if(alpha>pi*2||alpha<-pi*2)break;
		if(beta>pi*2||beta<-pi*2)break;
		if(gamma>pi*2||gamma<-pi*2)break;

		sa=sin(alpha);
		ca=cos(alpha);
		sb=sin(beta);
		cb=cos(beta);
		sc=sin(gamma);
		cc=cos(gamma);

		R[0]=ca*cb;
		R[1]=-ca*sb*sc-sa*cc;
		R[2]=-ca*sb*cc+sa*sc;
		R[3]=sa*cb;
		R[4]=-sa*sb*sc+ca*cc;
		R[5]=-sa*sb*cc-ca*sc;
		R[6]=sb;
		R[7]=cb*sc;
		R[8]=cb*cc;
		ss=0.0;
		for(j=0;j<9;j++)ss+=A[j]*R[j];
		if(ss>sumsqs)
		{
			sumsqs=ss;
			alphsave=alpha;
			betasave=beta;
			gammsave=gamma;
		}
		if(ss<0)i+=10;
	}
	
	/* ROTATE Q THROUGH THE ANGLES ALPH, BETA, GAMMA */
	qq=Q;
	for(i=N;i;i--,qq=(float*)((char*)qq+RMS_SHIFT))
	{
		x=R[0]*qq[0]+R[1]*qq[1]+R[2]*qq[2];
		y=R[3]*qq[0]+R[4]*qq[1]+R[5]*qq[2]; 
		z=R[6]*qq[0]+R[7]*qq[1]+R[8]*qq[2]; 
		qq[0]=x;
		qq[1]=y;
		qq[2]=z;
	}
	ss=0.0;
	pp=P;
	qq=Q;
	for(i=N;i;i--)
	{
		d=pp[0]-qq[0];
		ss+=d*d;
		d=pp[1]-qq[1];
		ss+=d*d;
		d=pp[2]-qq[2];
		ss+=d*d;
		pp=(float*)((char*)pp+RMS_SHIFT);
		qq=(float*)((char*)qq+RMS_SHIFT);
	}
	ss/=N;
	ss=sqrt(ss);
/*	
	if(ss>ss0)
	{
		count++;
		if(count<3)goto __RETRY;
		else return -ss0;
	}
*/	

	/* TRANSLATE BACK TO THE ORIGINAL P CENTROID */
	if(C)
	{	
		x=C[0];
		y=C[1];
		z=C[2];
		pp=P,qq=Q;
		for(i=N;i;i--)
		{
			pp[0]+=x;
			qq[0]+=x;
			pp[1]+=y;
			qq[1]+=y;
			pp[2]+=z;
			qq[2]+=z;
			pp=(float*)((char*)pp+RMS_SHIFT);
			qq=(float*)((char*)qq+RMS_SHIFT);
		}

		/* solve the correct C */
		B[0]=R[0]*QB[0]+R[1]*QB[1]+R[2]*QB[2]-x;
		B[1]=R[3]*QB[0]+R[4]*QB[1]+R[5]*QB[2]-y;
		B[2]=R[6]*QB[0]+R[7]*QB[1]+R[8]*QB[2]-z;
		x=R[4]*R[8]-R[5]*R[7];
		y=R[2]*R[7]-R[1]*R[8];
		z=R[1]*R[5]-R[2]*R[4];
		d=R[0]*x+R[3]*y+R[6]*z;
		d=1.0/d;
		x=B[0]*x+B[1]*y+B[2]*z;
		y=B[0]*(R[5]*R[6]-R[3]*R[8])+B[1]*(R[0]*R[8]-R[2]*R[6])
			+B[2]*(R[2]*R[3]-R[0]*R[5]);
		z=B[0]*(R[3]*R[7]-R[4]*R[6])+B[1]*(R[1]*R[6]-R[0]*R[7])
			+B[2]*(R[0]*R[4]-R[1]*R[3]);
		C[0]=x*d;
		C[1]=y*d;
		C[2]=z*d;
	}
	
	return ss;
}


