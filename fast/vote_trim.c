/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
Remove short aligned segments. The current trim routine has a bug probably.
-----------------------------------------------------------------------------*/

int	vote_trim(jz_pair_int *pairs,int len,int limit)
{
	jz_pair_int	*new_pairs;
	int		i,k,N;
	int		left,right,flag,inner;
	int		d1,d2;

	if(len<=0)return len;
	JZ_ARRAY_INIT(new_pairs,len+4);
	for(N=k=0;k<len;k++)
	{
		inner=flag=left=right=0;
		d1=d2=999;
		for(i=k-1;i>=0;i--)
		{
			d1=pairs[i+1][0]-pairs[i][0];
			d2=pairs[i+1][1]-pairs[i][1];
			if(d1==1&&d2==1)inner++;
			else break;
			if(inner>=limit)break;
		}

		if(d1+d2==3||d1+d2==4)
		{
			flag|=1;
			for(--i;i>=0;i--)
			{
				d1=pairs[i+1][0]-pairs[i][0];
				d2=pairs[i+1][1]-pairs[i][1];
				if(d1==1&&d2==1)left++;
				else break;
				if(left>=limit)break;
			}
		}
		
		d1=d2=999;
		for(i=k+1;i<len;i++)
		{
			d1=pairs[i][0]-pairs[i-1][0];
			d2=pairs[i][1]-pairs[i-1][1];
			if(d1==1&&d2==1)inner++;
			else break;
			if(inner>=limit)break;
		}
		if(d1+d2==3||d1+d2==4)
		{
			flag|=2;
			for(++i;i<len;i++)
			{
				d1=pairs[i][0]-pairs[i-1][0];
				d2=pairs[i][1]-pairs[i-1][1];
				if(d1==1&&d2==1)right++;
				else break;
				if(right>=limit)break;
			}
		}	
		
		inner++;
		if(inner<limit)
		{
			if(flag!=3)continue;
			if(inner<limit-1)continue;
			if(left<limit||right<limit)continue;
		}	
		new_pairs[N][0]=pairs[k][0];
		new_pairs[N++][1]=pairs[k][1];
	}

	JZ_ARRAY_COPYF(pairs,new_pairs,N);
	JZ_ARRAY_FREE(new_pairs);
	return N;
}

/*-----------------------------------------------------------------------------
Remove short aligned segments. The alignment must be trimed (no gap)
-----------------------------------------------------------------------------*/

int	vote_drop_short_segments(jz_pair_int *pairs,int len,int len_min)
{
	jz_pair_int		*new_pairs,*ppair;
	int			i,k,m,n;

	JZ_ARRAY_INIT(new_pairs,len);
	ppair=new_pairs;
	
	for(k=0;k<len;k++)
	{
		for(m=0,i=k-1;i>=0&&m<len_min;i--,m++)
			if(pairs[i][0]!=pairs[i+1][0]-1)break;
		for(i=k+1;i<len&&m<len_min;i++,m++)
			if(pairs[i][0]!=pairs[i-1][0]+1)break;
		
		for(n=0,i=k-1;i>=0&&n<len_min;i--,n++)
			if(pairs[i][1]!=pairs[i+1][1]-1)break;
		for(i=k+1;i<len&&n<len_min;i++,n++)
			if(pairs[i][1]!=pairs[i-1][1]+1)break;
		if(m>=len_min&&n>=len_min)
		{	
			(*ppair)[0]=pairs[k][0];
			(*ppair)[1]=pairs[k][1];
			ppair++;
		}		
	}

	len=ppair-new_pairs;
	JZ_ARRAY_COPYF(pairs,new_pairs,len);
	JZ_ARRAY_FREE(new_pairs);
	return len;
}

/*-----------------------------------------------------------------------------
 * Drop short aligned segments: there are at least four neighbors 
-----------------------------------------------------------------------------*/

int	vote_trim_old(jz_pair_int *pairs,int len)
{
	jz_pair_int	*new_pairs;
	int		i,k,m,n,N;
	int		left[2],right[2];

	if(len<=0)return len;
	JZ_ARRAY_INIT(new_pairs,len);
	for(N=k=0;k<len;k++)
	{
		/* go left to find out left[0] and left[1] */
		for(left[0]=left[1]=0,n=0,i=k-1;i>=0;i--)
		{
			m=pairs[i+1][0]+pairs[i+1][1]-pairs[i][0]-pairs[i][1];
			if(m==2)left[n]++;
			else if(m<=5)
			{
				n++;
				if(n<2)left[n]++;
				else break;
			}
			else break;
			if(left[n]>4)break;		
		}
		left[1]+=left[0];
		
		/* go right and get right[0] and right[1] */
		for(right[0]=right[1]=0,n=0,i=k+1;i<len;i++)
		{
			m=pairs[i][0]+pairs[i][1]-pairs[i-1][0]-pairs[i-1][1];
			if(m==2)right[n]++;
			else if(m<=5)
			{
				n++;
				if(n<2)right[n]++;
				else break;
			}
			else break;
			if(right[n]>4)break;
		}
		right[1]+=right[0];
		/*if(left[0]+right[0]<3)continue;*/
		if(left[1]+right[0]<4&&left[0]+right[1]<4)continue;
		new_pairs[N][0]=pairs[k][0];
		new_pairs[N++][1]=pairs[k][1];
	}

	JZ_ARRAY_COPYF(pairs,new_pairs,N);
	JZ_ARRAY_FREE(new_pairs);
	return N;
}


