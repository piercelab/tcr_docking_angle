/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
basic.c		Fundamental data structure support
=============================================================================*/


#include "basic.h"
#include <limits.h>

void    jz_array_list_free(void** array_list,int N)
{
	int k;
	if(!array_list)return;
	for(k=0;k<N;k++)
		if(array_list[k])JZ_LIST_FREE(array_list[k]);
	JZ_FREE(array_list);
}

void	*jz_array_resize(void **array,int size,int element)
{
	void	*p;
	
	if(*array==NULL)return NULL;
	if(size<=0)
	{
		JZ_FREE(*array);
		return (*array=NULL);
	}
	p=(void*)JZ_REALLOC(*array,element*size);
	return (p==NULL)?NULL:(*array=p);
}

int jz_chunk_grow(jz_chunk *chunk);

void    *jz_chunk_alloc(jz_chunk *chunk,int size)
{
	void	*ptr;
	
	if(chunk->available<size)
		if(jz_chunk_grow(chunk))return NULL;
	chunk->available-=size;
	ptr=(void*)(chunk->top);
	chunk->top+=size;
	return ptr;
}

void    jz_chunk_free(jz_chunk *chunk)
{
	jz_chunk *p;
	
	for(;chunk;)
	{
		p=(jz_chunk*)(chunk->next);
		free(chunk);
		chunk=p;
	}
}

int	jz_chunk_grow(jz_chunk *chunk)
{
	char *mem;
	mem=(char*)JZ_MALLOC(chunk->capacity+sizeof(void*));
	if(mem==NULL)return 1;
	*((void**)mem)=chunk->next;
	chunk->next=(void*)mem;
	chunk->available=chunk->capacity;
	chunk->top=mem+sizeof(void*);
	return 0;
}

jz_chunk*jz_chunk_init(int capacity)
{
	jz_chunk	*chunk;
	
	if(capacity<1)return NULL;
	chunk=(jz_chunk*)JZ_MALLOC(sizeof(jz_chunk)+capacity);
	if(!chunk)return NULL;
	chunk->next=NULL;
	chunk->top=(char*)(chunk)+sizeof(jz_chunk);
	chunk->capacity=capacity;
	chunk->available=capacity;
	return chunk;
}

int	jz_comp_int(const void *a,const void *b)
{
	return *((int*)(a))-*((int*)(b));
}

int	jz_comp_ptr(const void *a,const void *b)
{
	return (int)(a)-(int)(b);
}

int	jz_hash_primes[35]=
{
	11,	19,	37,	73,	109,
	163,	251,	367,	557,	823,
	1237,	1861,	2777,	4177,	6247,
	9371,	14057,	21089,	31627,	47431,
	71143,	106721,	160073,	240101,	360163,
	540217,	810343,	1215497,1823231,2734867,
	4102283,6153409,9230113,13845163,
};

static	const int jz_hash_prime_count
	=sizeof(jz_hash_primes)/sizeof(jz_hash_primes[0]);

int	jz_hash_closest_prime(int n)
{
	int k;
	for(k=0;k<jz_hash_prime_count;k++)
		if(jz_hash_primes[k]>n)return jz_hash_primes[k];
	return jz_hash_primes[jz_hash_prime_count-1];
}

int	jz_hashf_int(const void *ptr)
{
	return *((int*)(ptr));
}

int	jz_hashf_ptr(const void *ptr)
{
	return (int)(ptr);
}

int	jz_hashf_str(const char *s)
{
	int	ret=0;
	int	ch;
	
	for(;;)
	{
		ch=*s++;
		if(!ch)break;
		ret=(ret<<5)^(ch|ch<<7);
	}	
	return ret;
}

int	jz_hashf_strn_len=20;

int	jz_hashf_strn(const char *s)
{
	int	ret=0;
	int	ch;
	int	count;
	
	for(count=0;(ch=*s)&&(count<jz_hashf_strn_len);s++,count++)
		ret=(ret<<5)^(ch|ch<<7);
	return ret;
}

void	jz_hashl_each(void *set,void (*op)(void *),int element)
{
	int		k,slot,size=JZ_HASHL_SIZE(set);
	int		*lut=JZ_HASHL_LUT(set);
	int		*chain=JZ_HASHL_CHAIN(set);
	
	for(k=0;k<size;k++)
	{
		if((slot=lut[k])<0)continue;
		for(;;)
		{
			op(((char*)set)+slot*element);
			if((slot=chain[slot])<0)break;
		}
	}
}

/*----------------------------------------------------------------------------
initialize a hash list. 
----------------------------------------------------------------------------*/

void	*jz_hashl_init(jz_hashf hashf,jz_comp comp,int size,
	float expand,float shrink,float factor,int element)
{
	jz_hashl	*h;
	int		*p,k;
	
	size=jz_hash_closest_prime(size);
	if(expand<0)expand=JZ_HASHL_EXPAND_DEFAULT;
	if(shrink<0)shrink=JZ_HASHL_SHRINK_DEFAULT;
	if(factor<0)factor=JZ_HASHL_FACTOR_DEFAULT;
	k=(int)(size*expand);
	h=(jz_hashl*)JZ_MALLOC(sizeof(jz_hashl)+k*(element+sizeof(int))
		+size*sizeof(int)+3);
	if(h==NULL)return NULL;
	h->lut=(int*)((((long int)((char*)(h+1)+k*element))+3)^3);
	h->chain=h->lut+size;
	h->mincount=(int)(size*shrink)+1;
	if(h->mincount<=jz_hash_primes[0])h->mincount=0;
	h->maxcount=k;
	h->empty=-1;
	h->expand=expand;
	h->shrink=shrink;
	h->factor=factor;
	h->count=0;
	h->size=size;
	h->hashf=hashf?hashf:(jz_hashf)jz_hashf_int;
	h->comp=comp?comp:(jz_comp)jz_comp_int;
	for(p=h->lut,k=0;k<size;k++)*p++=-1;
	return h+1;
}

/*-----------------------------------------------------------------------------
read a hashl from a stream. return NULL on failure. the loaded 
maxcount is always equal to the count.
-----------------------------------------------------------------------------*/

void	*jz_hashl_read(FILE *fin,jz_hashf hashf,jz_comp comp,int element)
{
	jz_hashl	*hashl,header;
	int		size,err=-1;
	
	if(fin==NULL)return NULL;
	if(fread(&header,sizeof(jz_hashl),1,fin)!=1)return NULL;
	size=sizeof(jz_hashl)+element*header.maxcount;
	if(fseek(fin,-sizeof(jz_hashl),SEEK_CUR))return NULL;
	hashl=(jz_hashl*)JZ_MALLOC(size);
	if(hashl==NULL)return NULL;
	hashl->lut=NULL;
	hashl->chain=NULL;
	if(fread(hashl,size,1,fin)!=1)goto __EXIT;
	size=sizeof(int)*header.size;
	hashl->lut=(int*)JZ_MALLOC(size);
	if(hashl->lut==NULL)goto __EXIT;
	if(fread(hashl->lut,size,1,fin)!=1)goto __EXIT;
	size=sizeof(int)*header.maxcount;
	hashl->chain=(int*)JZ_MALLOC(size);
	if(hashl->chain==NULL)goto __EXIT;
	if(fread(hashl->chain,size,1,fin)!=1)goto __EXIT;
	err=0;
__EXIT:
	if(err)
	{
		if(!hashl)return NULL;
		if(hashl->lut)JZ_FREE(hashl->lut);
		if(hashl->chain)JZ_FREE(hashl->lut);
		return NULL;
	}
	hashl->hashf=hashf;
	hashl->comp=comp;
	return hashl+1;
}

/*-----------------------------------------------------------------------------
relink pointers to create a hashl, no checking for validity. end will points
to the end on success.
-----------------------------------------------------------------------------*/

void	*jz_hashl_relink(void *mem,void **end,int element)
{
	jz_hashl	*h;
	
	if(mem==NULL)return NULL;
	h=(jz_hashl*)mem;
	h->lut=(int*)((char*)h+sizeof(jz_hashl)+element*h->maxcount);
	h->chain=h->lut+h->size;
	*end=(void*)(h->chain+h->maxcount);
	return h+1;
}

/*----------------------------------------------------------------------------
resize a hash table. option == 0 for set and 1 for map. 
----------------------------------------------------------------------------*/

int	jz_hashl_resize(void **set,int option,int element)
{
	jz_hashl	*h;
	void		*news,*old;
	int		k,old_size,size,slot,count,*lut,*chain;
	char		*p;
	jz_hashf	hashf;
	jz_comp		comp;
	
	old=*set;
	count=JZ_HASHL_COUNT(old);
	old_size=JZ_HASHL_SIZE(old);
	size=jz_hash_closest_prime((int)(count*JZ_HASHL_FACTOR(old)));
	if(size==old_size)return 0;
	k=(int)(size*JZ_HASHL_EXPAND(old));
	h=(jz_hashl*)JZ_MALLOC(sizeof(jz_hashl)+k*(element+sizeof(int))
		+size*sizeof(int)+3);
	if(h==NULL)return 1;
	h->maxcount=k;
	h->lut=(int*)((((long int)((char*)(h+1)+k*element))+3)^3);
	h->chain=h->lut+size;
	h->mincount=(int)((float)size*(h->shrink=JZ_HASHL_SHRINK(old)))+1;
	if(h->mincount<=jz_hash_primes[0])h->mincount=0;
	h->empty=-1;
	h->expand=JZ_HASHL_EXPAND(old);
	h->factor=JZ_HASHL_FACTOR(old);
	h->count=count;
	h->size=size;
	h->hashf=hashf=JZ_HASHL_HASHF(old);
	h->comp=comp=JZ_HASHL_COMP(old);
	for(lut=h->lut,k=0;k<size;k++)*lut++=-1;
	news=h+1;
	
	if(JZ_HASHL_EMPTY(old)<0)
	{
		if(news!=old)memcpy(news,old,count*element);
	}
	else
	{
		p=(char*)news;
		lut=JZ_HASHL_LUT(old);
		chain=JZ_HASHL_CHAIN(old);
		if(element==sizeof(void*))
		{
			for(k=0;k<old_size;k++)
			{
				if((slot=lut[k])<0)continue;
				for(;;)
				{
					*((void**)p)=*((void**)((char*)old
						+slot*sizeof(void*)));
					p+=sizeof(void*);
					if((slot=chain[slot])<0)break;
				}
			}
		}
		else if(element==(sizeof(void*)*2))
		{
			void	**q;
			for(k=0;k<old_size;k++)
			{
				if((slot=lut[k])<0)continue;
				for(;;)
				{
					q=(void**)((char*)old+slot*
						(sizeof(void*)+sizeof(void*)));
					((void**)p)[0]=q[0];
					((void**)p)[1]=q[1];
					p+=sizeof(void*)+sizeof(void*);
					if((slot=chain[slot])<0)break;
				}
			}
		}
		else
		{
			for(k=0;k<old_size;k++)
			{
				if((slot=lut[k])<0)continue;
				for(;;)
				{
					memcpy(p,(char*)old+slot*element,
						element);
					p+=element;
					if((slot=chain[slot])<0)break;
				}
			}
		}
	}
	p=(char*)news;
	lut=h->lut;
	chain=h->chain;
	for(k=0;k<size;k++)lut[k]=-1;
	if(option==0)
		for(k=0;k<count;k++,p+=element)
		{
			slot=abs(hashf(p))%size;
			chain[k]=lut[slot];
			lut[slot]=k;
		}
	else
		for(k=0;k<count;k++,p+=element)
		{
			slot=abs(hashf(*((void**)p)))%size;
			chain[k]=lut[slot];
			lut[slot]=k;
		}
	*set=news;
	JZ_FREE(JZ_HASHL(old));
	return 0;
}

/*----------------------------------------------------------------------------
delete a key from the set. return the slot number of the deleted key, or -1
if not found.
----------------------------------------------------------------------------*/

int	jz_hashl_set_delete(void **set,void *node,int element)
{
	int		k;
	int		*chain,*p;
	void		*s=*set;
	jz_comp		comp;

	if((JZ_HASHL_COUNT(s))<JZ_HASHL_MINCOUNT(s))
	{
		jz_hashl_resize(set,0,element);
		s=*set;
	}
	p=JZ_HASHL_LUT(s)+abs(JZ_HASHL_HASHF(s)(node))%JZ_HASHL_SIZE(s);
	chain=JZ_HASHL_CHAIN(s);
	comp=JZ_HASHL_COMP(s);
	for(k=*p;k>=0;p=chain+k,k=*p)
		if(!comp(node,(char*)s+k*element))goto __FOUND;
	return -1;
__FOUND:
	*p=chain[k];
	chain[k]=JZ_HASHL_EMPTY(s);
	JZ_HASHL_EMPTY(s)=k;
	JZ_HASHL_COUNT(s)--;
	return k;
}

/*----------------------------------------------------------------------------
search a key in the set. return the index number or -1 if not found.
----------------------------------------------------------------------------*/

int	jz_hashl_set_find(void *set,void *node,int element)
{
	int		k;
	int		*chain;
	jz_comp		comp;
	
	k=JZ_HASHL_LUT(set)[abs(JZ_HASHL_HASHF(set)(node))%JZ_HASHL_SIZE(set)];
	if(k>=0)
	{
		chain=JZ_HASHL_CHAIN(set);
		comp=JZ_HASHL_COMP(set);
		do{}while(comp((char*)set+k*element,node)&&((k=chain[k])>=0));
	}
	return k;
}

/*----------------------------------------------------------------------------
insert a key into the set. return slot number if successiful, if 
a matched node is found, it is returned with match. return negative 
if the hash list is full but no new slot available.
----------------------------------------------------------------------------*/

int	jz_hashl_set_insert(void **set,void *node,int *match,int element)
{
	int		k,j,m,slot;
	int		*chain;
	jz_comp		comp;
	void		*s=*set;

	m=abs(JZ_HASHL_HASHF(s)(node));
	k=m%JZ_HASHL_SIZE(s);
	chain=JZ_HASHL_CHAIN(s);
	j=JZ_HASHL_LUT(s)[k];
	if(j<0)goto __NEXT;
	comp=JZ_HASHL_COMP(s);
	for(;;)
	{
		if(comp(node,(char*)s+j*element)==0)
		{
			*match=1;
			return j;
		}
		if((j=chain[j])<0)break;
	}
__NEXT:
	*match=0;
	if(JZ_HASHL_COUNT(s)>=JZ_HASHL_MAXCOUNT(s))
	{
		if(jz_hashl_resize(set,0,element))return -1;
		s=*set;
		chain=JZ_HASHL_CHAIN(s);
		k=m%JZ_HASHL_SIZE(s);
	}
	slot=JZ_HASHL_EMPTY(s);
	if(slot>=0)JZ_HASHL_EMPTY(s)=chain[slot];
	else slot=JZ_HASHL_COUNT(s);
	chain[slot]=JZ_HASHL_LUT(s)[k];
	JZ_HASHL_LUT(s)[k]=slot;
	JZ_HASHL_COUNT(s)++;
	return slot;
}

/*----------------------------------------------------------------------------
remove a key from the hash table. do not consider shrinking.
----------------------------------------------------------------------------*/

int	jz_hashl_set_remove(void *set,void *node,int element)
{
	int		k;
	int		*chain,*p;
	
	p=JZ_HASHL_LUT(set)+abs(JZ_HASHL_HASHF(set)(node))%JZ_HASHL_SIZE(set);
	chain=JZ_HASHL_CHAIN(set);
	for(k=*p;k>=0;p=chain+k,k=*p)
		if(JZ_HASHL_COMP(set)(node,(char*)set+k*element)==0)
			goto __FOUND;
	return -1;
__FOUND:
	*p=chain[k];
	JZ_HASHL_EMPTY(set)=k;
	JZ_HASHL_COUNT(set)--;
	return k;
}

/*----------------------------------------------------------------------------
search a key in the set. the matched node become the head of the bucket.
----------------------------------------------------------------------------*/

int	jz_hashl_set_search(void *set,void *node,int element)
{
	int		k,*entry,*p;
	int		*chain;
	jz_comp		comp;
	
	p=entry=JZ_HASHL_LUT(set)+abs(JZ_HASHL_HASHF(set)(node))
		%JZ_HASHL_SIZE(set);
	if((k=*entry)<0)return -1;
	chain=JZ_HASHL_CHAIN(set);
	comp=JZ_HASHL_COMP(set);
	for(;;)
	{
		if(!comp(node,(char*)set+k*element))
		{
			if(p!=entry)
			{
				*p=chain[k];
				chain[k]=*entry;
				*entry=k;
			}
			return k;
		}
		p=chain+k;
		k=*p;
		if(k<0)break;
	}
	return -1;
}

/*-----------------------------------------------------------------------------
write a hashl into a stream. return total number of bytes written, 
return -1 on failure. if concise!=0, then let maxcount=count;
-----------------------------------------------------------------------------*/

int	jz_hashl_write(void *hashl,FILE *fout,int concise,int element)
{
	size_t		size,actual;
	jz_hashl	header;
	
	header=*(JZ_HASHL(hashl));
	if(concise)header.maxcount=header.count;
	if(fwrite(&header,sizeof(jz_hashl),1,fout)!=1)return -1;
	size=header.maxcount*element;
	if(size>0&&fwrite(hashl,size,1,fout)!=1)return -1;
	actual=size;
	size=header.size*sizeof(int);
	if(fwrite(header.lut,size,1,fout)!=1)return -1;
	actual+=size;
	size=sizeof(int)*header.maxcount;
	if(size>0&&fwrite(header.chain,size,1,fout)!=1)return -1;
	actual+=size;
	return actual;
}

float	jz_list_factor=1.51;

void    *jz_list_addn(jz_list **list,int element,void *src,int count)
{
	jz_list *m=*list-1;
	int	size;
	char	*p;

	if(!src||(count<=0))return NULL;
	if((size=JZ_LIST_COUNT(*list)+count)>JZ_LIST_SIZE(*list))
	{
		size=(int)(size*jz_list_factor)+1;
		m=(jz_list*)JZ_REALLOC(m,element*size+sizeof(jz_list));
		if(m==NULL)return NULL;
		m->size=size;
		*list=m+1;
	}
	p=(char*)(m+1);
	p+=m->count*element;
	memcpy(p,src,count*element);
	m->count+=count;
	return m+1;
}

void	*jz_list_grow_(jz_list **list,int element,int *size)
{
	jz_list	*m;

	*size=(int)(*size*jz_list_factor)+1;
	m=(jz_list*)JZ_REALLOC(JZ_LIST(*list),element**size+sizeof(jz_list));
	if(m==NULL)return NULL;
	return *list=m+1;
}

void    *jz_list_grow(jz_list **list,int element)
{
	jz_list		*m;
	int		size;

	size=JZ_LIST_SIZE(*list);
	size=(int)(size*jz_list_factor)+1;
	m=(jz_list*)realloc(JZ_LIST(*list),element*size+sizeof(jz_list));
	if(m==NULL)return NULL;
	m->size=size;
	return (*list=(m+1));
}

void    *jz_list_init(int size,int element)
{
	jz_list	*h;
	
	if((size<1))return NULL;
	h=(jz_list*)JZ_MALLOC(element*size+sizeof(jz_list));
	if(!h)return NULL;
	h->size=size;
	h->count=0;
	return h+1;
}

void    jz_list_list_free(void **list_list)
{
	int k,count;
	if(!list_list)return;
	count=JZ_LIST_COUNT(list_list);
	for(k=0;k<count;k++)if(list_list[k])JZ_LIST_FREE(list_list[k]);
	JZ_LIST_FREE(list_list);
}

void    *jz_list_resize(jz_list **list,int size,int element)
{
	jz_list	*m;
	
	if(size<1)return NULL;
	m=(jz_list*)JZ_REALLOC(JZ_LIST(*list),element*size+sizeof(jz_list));
	if(m==NULL)
	{
		fprintf(stderr,"can not alloc\n");
		return NULL;
	}
	m->size=size;
	if(m->count>size)m->count=size;
	return *list=m+1;
}

void	*jz_matrix_init(int nrow,int ncol,int element)
{
	int	i,size;
	char	**m,*p;

	if((nrow<1)||(ncol<1))return NULL;
	m=(char**)JZ_MALLOC(element*nrow*ncol+nrow*sizeof(void*));
	if(m==NULL)return NULL;
	size=ncol*element;
	for(i=0,p=(char*)(m+nrow);i<nrow;i++,p+=size)m[i]=p;
	return m;
}

void	*jz_matrix_init_lh(int rl,int rh,int cl,int ch,int element)
{
	int	i,size;
	char	**m,*p;

	if((rh<rl)||(ch<cl)) return NULL;
	m=(char**)JZ_MALLOC(element*(rh-rl+1)*(ch-cl+1)
		+(rh-rl+1)*sizeof(*m));
	if(m==NULL)return NULL;
	size=(ch-cl+1)*element;
	m-=rl;
	for(i=rl,p=(char*)(m+rh+1)-cl*element;i<=rh;i++,p+=size)m[i]=p;
	return m;
}

void	*jz_mem_dup(void *mem,size_t size)
{
	void	*ret;

	if(!mem||!size)return NULL;
	ret=(void*)JZ_MALLOC(size);
	if(ret)memcpy(ret,mem,size);
	return ret;
}

void	jz_mem_map(void *mem,size_t size,unsigned char *map)
{
	unsigned char	*s=(char*)mem;
	size_t		n;
	
	if(size==0)return;
	n=size&3;
	size-=n;
	for(;size;size-=4)
	{
		s[0]=map[s[0]];
		s[1]=map[s[1]];
		s[2]=map[s[2]];
		s[3]=map[s[3]];
		s+=4;
	}
	if(n&2)
	{
		s[0]=map[s[0]];
		s[1]=map[s[1]];
		if(n&1)s[2]=map[(int)s[2]];
	}
	else if(n&1)s[0]=map[s[0]];
}



