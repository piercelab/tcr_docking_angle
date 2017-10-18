/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/* calculate the rmsd of two structures given an alignment */
float	vote_rmsd(int *pairs,int len,jz_protein *pt1,jz_protein *pt2)
{
	float		*P,*Q,*pp,*qq,*p,*q,ss;
	float		C[32],R[32];
	jz_residue	*R1,*R2;
	int		k;

	JZ_ARRAY_INIT(P,len*6);
	Q=P+len*3;
	R1=pt1->residues;
	R2=pt2->residues;
	pp=P;qq=Q;
	for(k=0;k<len;k++)
	{
		p=R1[pairs[k+k]].center;
		JZ_VEC3_CPY(pp,p);
		q=R2[pairs[k+k+1]].center;
		JZ_VEC3_CPY(qq,q);
		pp+=3;
		qq+=3;
	}
	
	ss=jz_rms_fitf(P,Q,len,C,R);
	if(ss<0.0)ss=-ss;
	JZ_ARRAY_FREE(P);
	return ss;
}


