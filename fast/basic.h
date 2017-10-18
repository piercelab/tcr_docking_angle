/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
basic.h		Fundamental data structure support
=============================================================================*/

#ifndef	__JZ_BASIC_H
#define	__JZ_BASIC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <inttypes.h>

#ifdef	__cplusplus

extern	"C" {

#endif

/*-----------------------------------------------------------------------------
basic operators used as macro parameters
-----------------------------------------------------------------------------*/

#define	JZ_OP2_MIN(a,b)			(((a)<(b))?(a):(b))
#define	JZ_OP2_MAX(a,b)			(((a)>(b))?(a):(b))

typedef	void(*jz_op)(const void *node);
typedef	int(*jz_comp)(const void *node1,const void *node2);

#define	JZ_CAST(var,type)		((type)(var))
#define	JZ_CASTP(var,type)		(((type)*)(&(var)))
#define	JZ_CASTR(var,type)		(*((type *)(&(var))))

#define	JZ_MALLOC(size)			((void*)(malloc(size)))

#define	JZ_MALLOCS(ptr,size,shift)					\
	((JZ_CASTR((ptr),void*)=(void*)JZ_MALLOC((size)+shift)))\
	?(JZ_CASTR((ptr),char*)=((char*)(ptr)+(shift))):NULL)

#define	JZ_CALLOC(nmemb,size)						\
	((void*)(calloc((nmemb),(size))))

#define	JZ_REALLOC(pointer,size)					\
	((void*)realloc((void*)(pointer),(size)))

#define	JZ_FREE(ptr)			((free((void*)(ptr))),0)

#define	JZ_FREEC(ptr)			((ptr)?JZ_FREE(ptr):0)

#define	JZ_RESIZE(ptr,size,buf)						\
	((JZ_CASTR((buf),void*)=(void*)JZ_REALLOC((ptr),(size)))\
	?(JZ_CASTR((ptr),void*)=(void*)(buf)):NULL)

int	jz_comp_int(const void*,const void*);

void	jz_mem_map(void *mem,size_t size,unsigned char *map);

void	jz_mem_revsb(void *mem,size_t size);

void	*jz_mem_sub(void *mem,size_t oft,size_t size);

void	*jz_mem_dup(void *mem,size_t size);


/*=============================================================================
JZ_ARRAY	Linear array container
=============================================================================*/


#define	JZ_ARRAY_INIT(array,size)					\
        (JZ_CASTR((array),void*)=JZ_MALLOC((size)*sizeof(*(array))))

#define	JZ_ARRAY_FREE(array)		JZ_FREE(array)		

/*----------------------------------------------------------------------------
change the size and lower bound of an array.
----------------------------------------------------------------------------*/	

#define	JZ_ARRAY_RESIZE(array,size,foo)					\
	JZ_RESIZE((array),((int)((size)*sizeof(*(array)))),(foo))
	
#define	JZ_ARRAY_NEW(array_,ptr_,size_,count_,factor_)			\
	do\
	{\
		if((count_)>=(size_))\
		{\
			void	*foo;\
			int	_my_size=(int)((size_)*(factor_))+1;\
			foo=(void*)realloc((array_),\
				_my_size*sizeof(*(array_)));\
			if(foo==NULL){(ptr_)=NULL;break;}\
			JZ_CASTR((array_),void*)=foo;\
			(size_)=_my_size;\
		}\
		JZ_CASTR((ptr_),void*)=(void*)((array_)+(count_));\
		(count_)++;\
	}while(0)

/*----------------------------------------------------------------------------
array duplication and copying
----------------------------------------------------------------------------*/	
	
#define	JZ_ARRAY_DUP(dest,src,size)					\
	(JZ_CASTR(dest,void*)=jz_mem_dup((src),(size)*sizeof(*(src))))

#define	JZ_ARRAY_DUPL(dest,src,size,low)				\
	((JZ_CASTR(dest,void*)=jz_mem_dup(((src)+(low)),\
	(size)*sizeof(*(src))))?((dest)-=low):NULL)
	
#define	JZ_ARRAY_COPYF(dest,src,size)					\
	memcpy((dest),(src),(size)*sizeof(*(src)))

#define	JZ_ARRAY_SETF(array,size,val)					\
	JZ_MEM_SETF((array),(size),sizeof(*(array)),&(val))
	
/*----------------------------------------------------------------------------
JZ_ARRAY_QSORTF		quick sort function, given comp function
----------------------------------------------------------------------------*/

#define	JZ_ARRAY_QSORTF(array,size,comp)				\
	qsort((array),(size),sizeof(*(array)),\
	(int(*)(const void*,const void*))(comp))

void	*jz_array_resize(void **array,int size,int element);


/*=============================================================================
jz_chunk      non-free allocator for small but large number of pieces.

a large chunk of memory is pre-alloced. each time a small piece is allocated 
from the chunk. it is assumed that the user never free a pointer until the 
whole chunk is destroyed. jz_chunk is defined as any pointer type. to avoid
misuse, it is recommended to be defined as void*.
=============================================================================*/

typedef	struct
{
	void		*next;
	int		capacity;
	int		available;
	char		*top;
}jz_chunk;

jz_chunk*jz_chunk_init(int capacity);

void	jz_chunk_free(jz_chunk *chunk);

void	*jz_chunk_alloc(jz_chunk *chunk,int size);

#define	JZ_CHUNK(chunk)			((jz_chunk*)(chunk))

#define	JZ_CHUNK_CAPACITY(chunk)	(JZ_CHUNK(chunk)->capacity)

#define	JZ_CHUNK_AVAILABLE(chunk)	(JZ_CHUNK(chunk)->available)

#define	JZ_CHUNK_NEXT(chunk)		(JZ_CHUNK(chunk)->next)

#define	JZ_CHUNK_TOP(chunk)		(JZ_CHUNK(chunk)->top)

#define	JZ_CHUNK_INIT(chunk,capacity)					\
	(*((jz_chunk**)(&(chunk)))=jz_chunk_init(capacity))

#define	JZ_CHUNK_FREE(chunk)		jz_chunk_free(JZ_CHUNK(chunk))

#define	JZ_CHUNK_ALLOC(chunk,size)      jz_chunk_alloc((chunk),(size))

#define	JZ_CHUNK_NEW(chunk,node)					\
	(*((void**)(node))=jz_chunk_alloc((chunk),sizeof(*(node))))


/*===========================================================================
jz_hashl	hash table with an array container 
===========================================================================*/

#define	JZ_HASHL_SHRINK_DEFAULT		0.382

#define	JZ_HASHL_EXPAND_DEFAULT		1.618

#define	JZ_HASHL_FACTOR_DEFAULT		1

typedef	int(*jz_hashf)(const void *key);

typedef	struct		jz_hashl
{
	jz_hashf	hashf;
	jz_comp		comp;
	float		shrink;
	float		expand;
	float		factor;
	int		empty;
	int		*chain;
	int		*lut;
	int		size;
	int		mincount;
	int		maxcount;
	int		count;
}jz_hashl;

#define	JZ_HASHL(hash)			((jz_hashl*)(hash)-1)

#define	JZ_HASHL_HASHF(hash)		(JZ_HASHL(hash)->hashf)

#define	JZ_HASHL_COMP(hash)		(JZ_HASHL(hash)->comp)

#define	JZ_HASHL_SHRINK(hash)		(JZ_HASHL(hash)->shrink)

#define	JZ_HASHL_EXPAND(hash)		(JZ_HASHL(hash)->expand)

#define	JZ_HASHL_FACTOR(hash)		(JZ_HASHL(hash)->factor)

#define	JZ_HASHL_EMPTY(hash)		(JZ_HASHL(hash)->empty)

#define	JZ_HASHL_CHAIN(hash)		(JZ_HASHL(hash)->chain)

#define	JZ_HASHL_LUT(hash)		(JZ_HASHL(hash)->lut)

#define	JZ_HASHL_SIZE(hash)		(JZ_HASHL(hash)->size)

#define	JZ_HASHL_MINCOUNT(hash)		(JZ_HASHL(hash)->mincount)

#define	JZ_HASHL_MAXCOUNT(hash)		(JZ_HASHL(hash)->maxcount)

#define	JZ_HASHL_COUNT(hash)		(JZ_HASHL(hash)->count)

#define	JZ_HASHL_INIT(hash,hashf,comp,size,expand,shrink,factor)	\
	(JZ_CASTR((hash),void*)=jz_hashl_init((jz_hashf)(hashf),\
	(jz_comp)(comp),(size),(expand),(shrink),(factor),sizeof(*(hash))))

#define	JZ_HASHL_FREE(hash)						\
	JZ_FREE(JZ_HASHL(hash))
	
#define	JZ_HASHLS_FIND(set,key)						\
	jz_hashl_set_find((void*)(set),(void*)(key),sizeof(*(set)))

#define	JZ_HASHLS_FIND_(set,key,hashf,comp,slot,dummy)			\
	do\
	{\
		slot=JZ_HASHL_LUT(set)[abs(hashf(key))%JZ_HASHL_SIZE(set)];\
		if(slot>=0)\
		{\
			int	*chain=JZ_HASHL_CHAIN(set);\
			do{dummy=set+slot;}while((comp(dummy,key))\
				&&((slot)=chain[slot])>=0);\
		}\
	}while(0)

#define	JZ_HASHLS_SEARCH(set,key)					\
	jz_hashl_set_search((void*)(set),(void*)(key),sizeof(*(set)))

#define	JZ_HASHLS_INSERT(set,key,match)					\
	jz_hashl_set_insert((void**)(&(set)),(void*)(key),\
	(int*)(&(match)),sizeof(*(set)))
		
#define	JZ_HASHLS_INSERT1(set,key,match,slot)				\
	((((slot)=jz_hashl_set_insert((void**)(&(set)),(void*)(&(key)),\
	(int*)(&(match)),sizeof(*(set))))>=0)?((set)[slot]=(key)):NULL)
		
#define	JZ_HASHLS_DELETE(set,key)					\
	jz_hashl_set_delete((void**)(&(set)),(void*)(key),sizeof(*(set)))

#define	JZ_HASHLS_REMOVE(set,key)					\
	jz_hashl_set_remove((void*)(set),(void*)(key),sizeof(*(set)))

#define	JZ_HASHL_WRITE(hash,fout,concise)				\
	jz_hashl_write((void*)(hash),(fout),(int)(concise),sizeof(*(hash)))

#define	JZ_HASHL_READ(hash,fin,hashf,comp)				\
	(JZ_CASTR((hash),void*)=jz_hashl_read((fin),\
	(jz_hashf)(hashf),(jz_comp)(comp),sizeof(*(hash))))

#define	JZ_HASHL_RELINK(hash,mem,end)				\
	(JZ_CASTR((hash),void*)=jz_hashl_relink((void*)(mem),\
	(void**)(&(end)),sizeof(*(hash))))

void	*jz_hashl_init(jz_hashf hashf,jz_comp comp,int size,
	float expand,float shrink,float factor,int element);

int	jz_hashl_resize(void **set,int option,int element);

int	jz_hashl_set_find(void *set,void *node,int element);

int	jz_hashl_set_search(void *set,void *node,int element);

int	jz_hashl_set_insert(void **set,void *node,int *match,int element);

int	jz_hashl_set_delete(void **set,void *node,int element);

int	jz_hashl_set_remove(void *set,void *node,int element);

void	jz_hashl_each(void *set,void (*op)(void *),int element);

int	jz_hashl_map_find(void *map,void *key,int element);

int	jz_hashl_map_search(void *map,void *key,int element);

int	jz_hashl_map_insert(void **map,void *key,int *match,int element);

int	jz_hashl_map_delete(void **map,void *key,int element);

int	jz_hashl_map_remove(void *map,void *key,int element);

int	jz_hashl_write(void *hashl,FILE *fout,int concise,int element);

void	*jz_hashl_read(FILE *fin,jz_hashf hashf,jz_comp comp,int element);

void	*jz_hashl_relink(void *mem,void **end,int element);

int	jz_hash_closest_prime(int n);
int	jz_hashf_int(const void *);
int	jz_comp_int(const void *,const void *);
int	jz_hashf_str(const char *);
int	jz_hashf_strn(const char *);
extern	int jz_hashf_strn_len;


/*===========================================================================
jz_list		autoincrement array. Beginning index is 0.
===========================================================================*/

typedef	struct
{
	int	size;
	int	count;
}jz_list;

#define	JZ_LIST(list)			((jz_list*)(list)-1)

#define	JZ_LIST_COUNT(list)		(((int*)(list))[-1])

#define	JZ_LIST_SIZE(list)		(((int*)(list))[-2])

#define	JZ_LIST_INIT(list,size)						\
	((JZ_MALLOCS((list),(size)*sizeof(*(list)),sizeof(jz_list))\
	?((((int*)(list))[-1]=0,((int*)(list))[-2]=(size))):0)

#define	JZ_LIST_FREE(list)		JZ_FREE(JZ_LIST(list))

#define	JZ_LIST_GROW(list)						\
	jz_list_grow((jz_list**)(&(list)),sizeof(*(list)))

#define	JZ_LIST_GROW_(list,size)					\
	jz_list_grow_((jz_list**)(&(list)),sizeof(*(list)),(int*)(&(size)))

#define	JZ_LIST_RESIZE(list,size)					\
	jz_list_resize((jz_list**)(&(list)),sizeof(*(list)),(size))

#define	JZ_LIST_ADD(list,value)						\
	do\
	{\
		int	*pc;\
		pc=(int*)(list)-1;\
		if(*pc>=pc[-1])\
		{\
			JZ_LIST_GROW(list);\
			pc=(int*)(list)-1;\
		}\
		(list)[(*pc)++]=(value);\
	}while(0)

#define	JZ_LIST_NEW(list,ptr)						\
	do\
	{\
		register int	*pc;\
		pc=(int*)(list)-1;\
		if(*pc>=pc[-1])\
		{\
			if(JZ_LIST_GROW(list)==NULL)\
			{\
				ptr=NULL;\
				break;\
			}\
			pc=(int*)(list)-1;\
		}\
		(void*)ptr=(void*)(list+(*pc)++);\
	}while(0)

#define	JZ_LIST_ADDN(list,src,count)					\
	jz_list_addn((jz_list**)(&list),sizeof(*(list)),(src),(count))

#define	JZ_ARRAY_LIST_FREE(array_list,N)				\
	jz_array_list_free((void**)(array_list),N)

#define	JZ_LIST_LIST_FREE(list_list)					\
	jz_list_list_free((void**)list_list)

void	*jz_list_init(int size,int element);

void	*jz_list_grow(jz_list**list,int element);

void	*jz_list_resize(jz_list**list,int element,int size);

void	*jz_list_addn(jz_list**list,int element,void *src,int count);

void	jz_array_list_free(void** array_list,int N);

void	jz_list_list_free(void **list_list);

extern	float jz_list_factor;


#define	JZ_MATRIX_ARRAY(matrix)		((matrix)[0])

#define	JZ_MATRIX_INIT(matrix,nrow,ncol)				\
	(*((void**)(&(matrix)))=\
	jz_matrix_init((nrow),(ncol),sizeof(**(matrix))))

#define	JZ_MATRIX_INITL(matrix,rl,rh,cl,ch)				\
	(*((void**)(&(matrix)))=\
	jz_matrix_init_lh((rl),(rh),(cl),(ch),sizeof(**(matrix))))

#define	JZ_MATRIX_FREE(matrix)		JZ_FREE(matrix)

#define	JZ_MATRIX_FREEL(matrix,rl)	JZ_FREE((matrix)+(rl))

void	*jz_matrix_init(int nrow,int ncol,int element);
void	*jz_matrix_init_lh(int rl,int rh,int cl,int ch,int element);


typedef	struct	jz_pair
{
	void	*first;
	void	*second;
}jz_pair;

typedef	int	jz_pair_int[2];

#ifdef	__cplusplus
}
#endif

#endif


