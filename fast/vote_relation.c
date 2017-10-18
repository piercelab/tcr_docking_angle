/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <math.h>

#include "vote.h"

/*-----------------------------------------------------------------------------
Get unit vectors. VECS[i] is the unit vector from residue i-1 to residue i.
Size of the array is N+1, in which:
VECS[0] is defined as the same as VECS[1],
VECS[N] is defined as the same as VECS[N-1].
-----------------------------------------------------------------------------*/

void	vote_vector(int *VECS,int *PTD,int N)
{
	int		d;

	JZ_VEC3_SUB(VECS,PTD+3,PTD);
	VOTE_UNIT_VECTOR(VECS,d);
	for(VECS+=3,--N;N;N--,VECS+=3,PTD+=3)
	{
		JZ_VEC3_SUB(VECS,PTD+3,PTD);
		VOTE_UNIT_VECTOR(VECS,d);
	}
	JZ_VEC3_CPY(VECS,VECS-3);
}

/*-----------------------------------------------------------------------------
For now, let me skip masks
-----------------------------------------------------------------------------*/

void	vote_relation(Relation *RR,int *PTD,int N,int option)
{
	Relation	*prr,*qrr;
	int		*p,*q,*r,*s;
	int		i,j,d,X[3];
	int		*VECS;
	int16_t		*acos_lut;
	uint16_t	*distance_lut;
	uint16_t	*angle_lut,*angle_lut2;

	JZ_ARRAY_INIT(VECS,3*(N+2));
	vote_vector(VECS,PTD,N);
	acos_lut=vote_acos_lut+VOTE_ACOS_LUT_SHIFT;
	if(option)
	{
		distance_lut=vote_distance_mask_lut;
		angle_lut=vote_angle_mask_lut;
		angle_lut2=vote_angle_mask_lut2;
	}	
	else
	{
		distance_lut=vote_distance_bit_lut;
		angle_lut=vote_angle_bit_lut;
		angle_lut2=vote_angle_bit_lut;
	}	

	/* get relations between residue i and j */
	for(i=0,p=VECS,r=PTD;i<N;i++,p+=3,r+=3)
	{
		prr=RR+i*(N+1);
		qrr=prr;
		q=p;
		s=r;
		
		prr->mask[0]=0;
		prr++;
		qrr+=N;
		s+=3;
		q+=3;

		for(j=i+1;j<N;j++,prr++,qrr+=N,s+=3)
		{
			prr->mask[7]=VC_RR_ORDER_LT;
			qrr->mask[7]=VC_RR_ORDER_GT;
			
			/* distance */
			JZ_VEC3_SUB(X,s,r);
			JZ_VEC3_DOT(d,X,X);
			if(d<0)d=INT_MAX;
			d=(int)(sqrt(d)+1);
			prr->distance=qrr->distance=d;
			prr->mask[0]=qrr->mask[0]=distance_lut[d>>2];

			if(d>VC_RR_DISTANCE_MAX)
			{
				q+=3;
				prr->mask[0]=qrr->mask[0]=0;
				continue;
			}

			/* make X unit vector */
			JZ_VEC3_MUL(X,X,VOTE_DISCRETIZE_UNIT_VECTOR);
			JZ_VEC3_DIV(X,X,d);
	
			/* angles: see 010603-16:02 for definition */
			/* redefined in 012303-19:05 */
			JZ_VEC3_DOT(d,p,q);
			d>>=VOTE_ACOS_SHIFT;
			/*
			if(d<-2048)d=-2048;
			if(d>2048)d=2048;
			*/
			d=acos_lut[d];
			prr->angles[4]=d;
			prr->mask[5]=angle_lut2[d];

			JZ_VEC3_DOT(d,p,X);
			d>>=VOTE_ACOS_SHIFT;
			if(d<-2048)d=-2048;
			if(d>2048)d=2048;
			d=acos_lut[d];
			prr->angles[0]=d;
			prr->mask[1]=angle_lut[d];

			JZ_VEC3_DOT(d,q,X);
			d>>=VOTE_ACOS_SHIFT;
			if(d<-2048)d=-2048;
			if(d>2048)d=2048;
			d=acos_lut[d];
			prr->angles[1]=d;
			prr->mask[2]=angle_lut[d];

			q+=3;
			JZ_VEC3_DOT(d,p+3,q);
			d>>=VOTE_ACOS_SHIFT;
			if(d<-2048)d=-2048;
			if(d>2048)d=2048;
			d=acos_lut[d];
			prr->angles[5]=d;
			prr->mask[6]=angle_lut2[d];

			JZ_VEC3_DOT(d,p+3,X);
			d>>=VOTE_ACOS_SHIFT;
			if(d<-2048)d=-2048;
			if(d>2048)d=2048;
			d=acos_lut[d];
			prr->angles[2]=d;
			prr->mask[3]=angle_lut[d];

			JZ_VEC3_DOT(d,q,X);
			d>>=VOTE_ACOS_SHIFT;
			if(d<-2048)d=-2048;
			if(d>2048)d=2048;
			d=acos_lut[d];
			prr->angles[3]=d;
			prr->mask[4]=angle_lut[d];
			
			/* copy angles to qrr */
			qrr->angles[0]=prr->angles[1];
			qrr->angles[1]=prr->angles[0];
			qrr->angles[2]=prr->angles[3];
			qrr->angles[3]=prr->angles[2];
			qrr->angles[4]=prr->angles[4];
			qrr->angles[5]=prr->angles[5];
			qrr->mask[1]=prr->mask[2];
			qrr->mask[2]=prr->mask[1];
			qrr->mask[3]=prr->mask[4];
			qrr->mask[4]=prr->mask[3];
			qrr->mask[5]=prr->mask[5];
			qrr->mask[6]=prr->mask[6];
		}
	}

	JZ_ARRAY_FREE(VECS);
}


