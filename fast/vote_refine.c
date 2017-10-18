/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
Refine alignment with VOTE
-----------------------------------------------------------------------------*/

void	vote_refine_allow(char *A,int *VOTE,int *pairs,int len,
	int *pt1,int *pt2,int M,int N)
{
	int		*pp,*p,*q,*pg;
	int		CA[12],CB[12],U[3],V[3];
	int		D1_max[4],D1_min[4];
	int		*D2,*pd;
	int		i,j,k,m,n;
	int		d,a;
	char		*pa;
	int		*LA,*LB;
	int		LB_max,LB_min;
	int		bad_score,worse_score;

	/* initialize the Allow matrix */
	jz_charray_set(A,1,M*N);
	
	/* define LA and LB as in 122103-17:56 */
	JZ_ARRAY_INIT(LA,M+N);
	LB=LA+M;
	for(a=d=k=0;k<len;k++)
	{
		m=pairs[k+k];
		n=pairs[k+k+1];
		for(i=a;i<=m;i++)LA[i]=k;
		for(j=d;j<=n;j++)LB[j]=k;
		a=m+1;d=n+1;
	}
	for(i=a;i<M;i++)LA[i]=len;
	for(j=d;j<N;j++)LB[j]=len;
	
	/* the alignment itself is always allowed */
	for(k=0;k<len;k++)
	{
		i=pairs[k+k];
		j=pairs[k+k+1];
		m=i*N+j;
		A[m]=0;
		/*if(i>0)A[m-N]=0;*/
		/*if(j>0)A[m-1]=0;*/
		if(i>0&&j>0)A[m-N-1]=0;
		/*if(i<M-1)A[m+N]=0;*/
		/*if(j<N-1)A[m+1]=0;*/
		if(i<M-1&&j<N-1)A[m+N+1]=0;
	}	
	if(len<8)return;
	
	/* calculate the centers */
	for(i=0,pp=pairs;i<4;i++)
	{
		if(i==3)n=len-(len/4)*3;
		else n=len/4;
		JZ_VEC3_SET(U,0);
		JZ_VEC3_SET(V,0);
		for(j=n;j;j--,pp+=2)
		{
			p=pt1+pp[0]*3;
			JZ_VEC3_ADD(U,U,p);
			q=pt2+pp[1]*3;
			JZ_VEC3_ADD(V,V,q);
		}
		j=i*3;
		JZ_VEC3_CPY(CA+j,U);
		JZ_VEC3_CPY(CB+j,V);
	}

	/* define three pairs of centers */
	i=(len/4)*2;
	j=len-i;
	p=CA,q=CA+6;
	JZ_VEC3_ADD(U,p,q);
	JZ_VEC3_DIV(U,U,i);
	q=CA+3;
	JZ_VEC3_ADD(p,p,q);
	JZ_VEC3_DIV(p,p,i);
	p+=3,q+=3;
	JZ_VEC3_ADD(p,p,q);
	JZ_VEC3_DIV(p,p,i);
	p+=3,q+=3;
	JZ_VEC3_ADD(p,p,q);
	JZ_VEC3_DIV(p,p,j);
	JZ_VEC3_CPY(q,U);
	p=CB,q=CB+6;
	JZ_VEC3_ADD(U,p,q);
	JZ_VEC3_DIV(U,U,i);
	q=CB+3;
	JZ_VEC3_ADD(p,p,q);
	JZ_VEC3_DIV(p,p,i);
	p+=3,q+=3;
	JZ_VEC3_ADD(p,p,q);
	JZ_VEC3_DIV(p,p,i);
	p+=3,q+=3;
	JZ_VEC3_ADD(p,p,q);
	JZ_VEC3_DIV(p,p,j);
	JZ_VEC3_CPY(q,U);

	/* calculate distance to centers, in protein pt2 */
	JZ_ARRAY_INIT(D2,N*4);
	for(p=pt2,pd=D2,i=N;i;i--,p+=3,pd+=4)
	{
		JZ_VEC3_SUBSQ(d,p,CB);
		pd[0]=(int)sqrt(d);
		JZ_VEC3_SUBSQ(d,p,CB+3);
		pd[1]=(int)sqrt(d);
		JZ_VEC3_SUBSQ(d,p,CB+6);
		pd[2]=(int)sqrt(d);
		JZ_VEC3_SUBSQ(d,p,CB+9);
		pd[3]=(int)sqrt(d);
	}

	/* Find the longest center-center vector */
	

	
	
	
	/* check if allowed */
	bad_score=-16*len;
	worse_score=VOTE[N-1];
	for(i=0,p=pt1,pa=A,pg=VOTE;i<M;i++,pa+=N,pg+=N,p+=3)
	{
		JZ_VEC3_SUBSQ(d,p,CA);
		a=(int)sqrt(d);
		D1_max[0]=VC_REFINE_ALLOW_MAX(a);
		D1_min[0]=VC_REFINE_ALLOW_MIN(a);
		JZ_VEC3_SUBSQ(d,p,CA+3);
		a=(int)sqrt(d);
		D1_max[1]=VC_REFINE_ALLOW_MAX(a);
		D1_min[1]=VC_REFINE_ALLOW_MIN(a);
		JZ_VEC3_SUBSQ(d,p,CA+6);
		a=(int)sqrt(d);
		D1_max[2]=VC_REFINE_ALLOW_MAX(a);
		D1_min[2]=VC_REFINE_ALLOW_MIN(a);
		JZ_VEC3_SUBSQ(d,p,CA+9);
		a=(int)sqrt(d);
		D1_max[3]=VC_REFINE_ALLOW_MAX(a);
		D1_min[3]=VC_REFINE_ALLOW_MIN(a);
		
		LB_max=LA[i]+VC_REFINE_SEQ_TH(len);
		LB_min=LA[i]-VC_REFINE_SEQ_TH(len);
		
		for(j=0,pd=D2;j<N;j++,pd+=4)
		{
			/*if(pg[j]<=bad_score&&pg[j]!=worse_score)continue;*/
			if(LB[j]>LB_max||LB[j]<LB_min)continue;
			if(pd[0]>D1_max[0])continue;
			if(pd[0]<D1_min[0])continue;
			if(pd[1]>D1_max[1])continue;
			if(pd[1]<D1_min[1])continue;
			if(pd[2]>D1_max[2])continue;
			if(pd[2]<D1_min[2])continue;
			if(pd[3]>D1_max[3])continue;
			if(pd[3]<D1_min[3])continue;
			pa[j]=0;
		}
	}

	/* trim isolated */
	n=N+1;
	m=n+n;
	pa=A;
	for(i=3;i<M-3;i++,pa+=N)
	for(j=3;j<N-3;j++)
	{
		k=0;
		if(!pa[j])k++;
		if(!pa[j-n])k++;
		if(!pa[j-m])k++;
		if(!pa[j-m-n])k++;
		if(!pa[j+n])k++;
		if(!pa[j+m])k++;
		if(!pa[j+m+n])k++;
		if(k<4)pa[j]=1;
		else if(pa[j])pa[j]=2;
	}
	for(i=0,n=M*N;i<n;i++)
		if(A[i]==2)A[i]=0;

	JZ_ARRAY_FREE(D2);
}

/*-----------------------------------------------------------------------------
Vote with the current alignment
-----------------------------------------------------------------------------*/

void	vote_refine_vote(int *VOTE,char *ALLOW,int *PAIRS,int len,
	Relation *RR1,Relation *RR2,int M,int N,int ***cache,
	int ***new_cache,int *map)
{
	Relation	*R2,*RRP1,*RRP2;
	int		i,j,k;
	int		a,b,s,d,sum;
	char		*pa;
	int		*pv;
	int		**pci,**pnci;
	int		*pcij,*pncij;

	/* Initialize the vote matrix */
	k=VC_REFINE_BAD_SCORE*len;
	jz_intarray_set(VOTE,k,M*N);
	
	for(i=0,pa=ALLOW,pv=VOTE;i<M;i++,pa+=N,pv+=N,RR1+=M)
	for(pci=cache[i],pnci=new_cache[i],j=0,R2=RR2;j<N;j++,R2+=N)
	{
		if(pa[j])
		{
			pnci[j]=NULL;
			continue;
		}
		JZ_ARRAY_INIT(pnci[j],len);
		pncij=pnci[j];
		pcij=pci[j];
	
		/* the entry is not active during last round */
		if(pcij==NULL)
		{
			for(sum=k=0;k<len;k++)
			{
				if((PAIRS[k+k]==i)||(PAIRS[k+k+1]==j))
				{
					pncij[k]=0;
					continue;
				}
					
				RRP1=RR1+PAIRS[k+k];
				RRP2=R2+PAIRS[k+k+1];
				VOTE_SCORE_PAIR(s,RRP1,RRP2,a,b,d);
				pncij[k]=s;
				sum+=s;
			}	
		}
		
		else
		{
			for(sum=k=0;k<len;k++)
			{
				/* the aligned pair is in cache */
				if(map[k]>=0)
					s=pcij[map[k]];
				else
				{
					if((PAIRS[k+k]==i)
						||(PAIRS[k+k+1]==j))
					{
						pncij[k]=0;
						continue;
					}
					
					RRP1=RR1+PAIRS[k+k];
					RRP2=R2+PAIRS[k+k+1];
					VOTE_SCORE_PAIR
						(s,RRP1,RRP2,a,b,d);	
				}
				pncij[k]=s;
				sum+=s;	
			}
		}
		pv[j]=sum+len*8;
	}

	VOTE[N-1]=VC_REFINE_BAD_SCORE*len;
}

/*-----------------------------------------------------------------------------
The refinement algorithm
-----------------------------------------------------------------------------*/

int	vote_refine(int *PAIRS,int len,int *VOTE,
	int *pt1,int *pt2,Relation *RR1,Relation *RR2,
	int M,int N,int *Score)
{
	int		*PAIRS_SAVE;
	int		k,n;
	int		gap_a,gap_b;
	int		len_save;
	int		score_save;
	int		score;
	char		*ALLOW;
	int		*PAIRS_OLD;
	int		len_old;
	int		*map;
	int		***cache1,***cache2;
	int		***cache,***cache_new,***cache_tmp;
	
	/* initialize refine_allow matrix */
	JZ_ARRAY_INIT(ALLOW,M*N);
	
	/* initialize caches */
	JZ_MATRIX_INIT(cache1,M,N);
	JZ_MATRIX_INIT(cache2,M,N);
	cache=cache1;
	cache_new=cache2;
	for(k=0,n=M*N;k<n;k++)cache1[0][k]=NULL;
	JZ_ARRAY_COPYF(cache2[0],cache1[0],M*N);

	/* initialize map */
	JZ_ARRAY_INIT(map,M+N);
	
	/* make a copy to the old alignment */
	JZ_ARRAY_INIT(PAIRS_OLD,M+N);
	JZ_ARRAY_COPYF(PAIRS_OLD,PAIRS,len+len);
	len_old=len;

	/* get initial cache */

	vote_cache_initial(cache,RR1,RR2,PAIRS,len,M,N);
	
	/* unit map */
	for(k=0;k<len;k++)map[k]=k;
	
	/* shrink */
	len=vote_shrink(PAIRS,len,M,N,RR1,RR2,cache1,map,&score);

	/* save the shrinked alignment, as the starting point */
	JZ_ARRAY_INIT(PAIRS_SAVE,M+N);
	JZ_ARRAY_COPYF(PAIRS_SAVE,PAIRS,len+len);
	len_save=len;
	score_save=len;

	jz_intarray_set(VOTE,INT_MAX,M*N);
	
	/* refinement loop */
	for(k=0;k<VC_REFINE_ROUND;k++)
	{
		if(len<=8)break;

		/* find the current map */
		vote_map(map,PAIRS,len,PAIRS_OLD,len_old);

		/* find alignable pairs given the current alignment */
		vote_refine_allow(ALLOW,VOTE,PAIRS,len,pt1,pt2,M,N);
	
/*
fprintf(stderr,"refine round %d, allowd entry = %d, alignment length = %d\n",
	k,jz_charray_count_zero(ALLOW,M*N),len);
*/	

		/* vote with the current alignment */
		vote_refine_vote(VOTE,ALLOW,PAIRS,len,RR1,RR2,M,N,
			cache,cache_new,map);
		
		/* cleanup cache */
		vote_cache_clean(cache[0],M,N);
		
		/* swap the cache variable, renew alignment */
		cache_tmp=cache;
		cache=cache_new;
		cache_new=cache_tmp;
		JZ_ARRAY_COPYF(PAIRS_OLD,PAIRS,len+len);
		len_old=len;
		
		/* DP to find out the best sequential alignment */
		gap_a=VC_REFINE_GAPA*len;
		gap_b=VC_REFINE_GAPB*len;
		len=jz_dp_local_affine(VOTE,M,N,gap_a,gap_b,
			(jz_pair_int *)PAIRS,0,&score);
	
		len=vote_extend((jz_pair_int*)PAIRS,len,M,N);
		/* include non-sequential pairs */
		
		/* shrink the alignment to remove bad pairs */
		vote_map(map,PAIRS,len,PAIRS_OLD,len_old);

		len=vote_shrink(PAIRS,len,M,N,RR1,RR2,cache,map,&score);
		len=vote_trim((jz_pair_int*)PAIRS,len,2);
		vote_score(PAIRS,len,M,N,RR1,RR2);

		if(len==len_save&&score==score_save)break;
		
		/* test if score is improved */
		if(score+VC_SCORE_LEN*len>score_save+VC_SCORE_LEN*len_save)
		{
			JZ_ARRAY_COPYF(PAIRS_SAVE,PAIRS,len+len);
			len_save=len;
			score_save=score;

			/*
			fprintf(stderr,"len=%d, score=%d\n",len,score);
			vote_view_align(PAIRS_SAVE,"ALIGN",M,N,len,
				jz_color_red);
			*/	
		}
	}

	JZ_ARRAY_COPYF(PAIRS,PAIRS_SAVE,len_save+len_save);
	*Score=score_save;
	JZ_ARRAY_FREE(PAIRS_SAVE);
	JZ_ARRAY_FREE(PAIRS_OLD);
	JZ_ARRAY_FREE(map);
	JZ_ARRAY_FREE(ALLOW);
	vote_cache_clean(cache1[0],M,N);
	vote_cache_clean(cache2[0],M,N);
	JZ_MATRIX_FREE(cache1);
	JZ_MATRIX_FREE(cache2);
	return len_save;
}


