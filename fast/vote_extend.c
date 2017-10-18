/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "basic.h"

/*-----------------------------------------------------------------------------
Extend an alignment, then call shink()
-----------------------------------------------------------------------------*/

int	vote_extend(jz_pair_int *pairs,int len,int M,int N)
{
	jz_pair_int	*new_pairs;
	int		k,L,flag=1;
	int		last_x=-1;
	int		last_y=-1;

	if(len<=0)return len;
	JZ_ARRAY_INIT(new_pairs,len*8+4);

	/* extend to right */
	for(L=k=0;k<len;k++)
	{
		if(flag&&(pairs[k][0]-1>last_x&&pairs[k][1]-1>last_y))
		{
			new_pairs[L][0]=pairs[k][0]-1;
			new_pairs[L][1]=pairs[k][1]-1;
			last_x=new_pairs[L][0];
			last_y=new_pairs[L][1];
			L++;
		}	
		flag=0;
		new_pairs[L][0]=pairs[k][0];
		new_pairs[L][1]=pairs[k][1];
		last_x=new_pairs[L][0];
		last_y=new_pairs[L][1];
		L++;
		if(pairs[k][0]>=M-1||pairs[k][1]>=N-1)continue;
		if(k<len-1)
		{
			if(pairs[k+1][0]==pairs[k][0]+1)continue;
			if(pairs[k+1][1]==pairs[k][1]+1)continue;
		}
		new_pairs[L][0]=pairs[k][0]+1;
		new_pairs[L][1]=pairs[k][1]+1;
		last_x=new_pairs[L][0];
		last_y=new_pairs[L][1];
		L++;
		flag=1;
	}

	JZ_ARRAY_COPYF(pairs,new_pairs,L);
	JZ_ARRAY_FREE(new_pairs);
	return L;
}

