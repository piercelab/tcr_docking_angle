/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <strings.h>

#include "basic.h"
#include "chore.h"

/*----------------------------------------------------------------------------
print command line parsing error message
----------------------------------------------------------------------------*/

void	jz_cmd_message(FILE *stream)
{
	switch(jz_cmd_err)
	{
	case	JZ_CMD_SUCCESS:break;
	/*	fprintf(stream,"success\n");break;*/
	case	JZ_CMD_ENUMBER:
		fprintf(stream,"ERROR:unexpected argument\n");
		break;
	case	JZ_CMD_EARGUMENT:
		fprintf(stream,"ERROR:Invalid argument %s\n",jz_cmd_arg);
		break;
	case	JZ_CMD_EPARAM:
		fprintf(stream,"ERROR:Parameter absent for %s\n",jz_cmd_arg);
		break;
	case	JZ_CMD_ETYPE:
		fprintf(stream,"ERROR:Invalid parameter type in '%s'\n",
			jz_cmd_param);
		break;
	case	JZ_CMD_EMULTIPLE:
		fprintf(stream,
			"ERROR:Repeated use of option %s or equivalents\n",
			jz_cmd_arg);
		break;
	}
}


int	jz_cmd_pass_param(void *addr,int type,int multiple,char *str);

int	jz_cmd_pass_param(void *addr,int type,int multiple,char *str)
{
	char *end;
	int i;
	double d;
	
	switch(type)
	{
	case	JZ_TBOOL:
		if(*((int*)addr)==0)*((int*)addr)=1;
		else *((int*)addr)=0;
		break;
	case	JZ_TSTR:
		if(multiple)
		{
			if(*((char***)addr)==NULL)
				JZ_LIST_INIT((*((char***)addr)),1);
			JZ_LIST_ADD((*((char***)addr)),str);
		}
		else *((char**)addr)=str;
		break;
	case	JZ_TINT:
		errno=0;
		i=strtol(str,&end,10);
		if(*end||errno)return JZ_CMD_ETYPE;
		if(multiple)
		{
			if(*((int**)addr)==NULL)
				JZ_LIST_INIT((*((int**)addr)),1);
			JZ_LIST_ADD((*((int**)addr)),i);
		}
		else *((int*)(addr))=i;
		break;
	case	JZ_TFLOAT:
		errno=0;
		d=strtod(str,&end);
		if(*end||errno)return JZ_CMD_ETYPE;
		if(multiple)
		{
			if(*((float**)addr)==NULL)
				JZ_LIST_INIT((*((float**)addr)),1);
			JZ_LIST_ADD((*((float**)addr)),(float)d);
		}
		else *((float*)addr)=(float)d;
		break;
	case	JZ_TDOUBLE:
		errno=0;
		d=strtod(str,&end);
		if(*end||errno)return JZ_CMD_ETYPE;
		if(multiple)
		{
			if(*((double**)addr)==NULL)
				JZ_LIST_INIT((*((double**)addr)),1);
			JZ_LIST_ADD((*((double**)addr)),d);
		}
		else *((double*)addr)=d;
		break;
	case	JZ_TCHAR:
		if(multiple)
		{
			if(*((char**)addr)==NULL)
				JZ_LIST_INIT((*((char**)addr)),1);
			JZ_LIST_ADD((*((char**)addr)),*str);
		}
		else *((char*)addr)=*str;
		break;
	case	JZ_TMESSAGE:
		if(*((char**)addr))printf(*((char**)addr));
		exit(0);
	default:return JZ_CMD_ETYPE;
	}
	return 0;
}
	
/*----------------------------------------------------------------------------
parse all command line arguments
----------------------------------------------------------------------------*/

int	jz_cmd_err=0;
char	*jz_cmd_arg=NULL;
char	*jz_cmd_param=NULL;

int	jz_cmd_parse(int argc,char **argv,jz_cmd *opts,int opt_count,
	int option,jz_pair **pairs,char ***strs)
{
	int	ret=0;
	int	k,pos,len,arg_len,opt,match;
	int	prefix_type;
	char	*item=NULL,*p,*q,*s,*param=NULL;
	int	addvar;
	
	if(pairs&&(*pairs==NULL))
		JZ_HASHL_INIT(*pairs,jz_hashf_str,strcmp,9,-1,-1,-1);
	if(strs&&(*strs==NULL))JZ_LIST_INIT(*strs,4);
	
	for(pos=1;pos<argc;)
	{
		item=argv[pos++];
		len=strlen(item);
		if(!strncmp(item,JZ_LONG_PREFIX,JZ_LONG_PREFIX_LENGTH))
			prefix_type=2;
		else if((item[0]==JZ_SHORT_PREFIX)&&!isdigit(item[1]))
			prefix_type=1;
		else prefix_type=0;
		param=NULL;
		p=item;
		if(prefix_type)
		{
			p+=((prefix_type==1)?1:JZ_LONG_PREFIX_LENGTH);
			q=strchr(p,'=');
			if(q)
			{
				param=q+1;
				*q=0;
			}
		}
		len=strlen(p);
		if(len==0)goto __EXIT;
		addvar=0;
		for(k=0;k<opt_count;k++)
		{
			opt=opts[k].option;
			if(!(opt&JZ_CMD_ARGUMENT)&&!prefix_type)
			{
				if(opts[k].addr==NULL)continue;
				if(jz_cmd_pass_param(opts[k].addr,
					opts[k].type,opt&JZ_CMD_MULTIPLE,p))
				{
					ret=JZ_CMD_ETYPE;
					goto __EXIT;
				}
				if(!(opt&JZ_CMD_MULTIPLE))opts[k].addr=NULL;
				goto __DONE;
			}
			q=opts[k].arg;
			if(q==NULL)continue;
			for(match=0;;)
			{
				s=strchr(q,'|');
				if(s==NULL)arg_len=strlen(q);
				else arg_len=s-q;
				if(len==arg_len)
				{
					if(((len==1)&&(*q==*p))
						||!strncmp(p,q,len))
					{
						match=1;
						break;
					}
				}
				if(s==NULL)break;
				q=s+1;
			}
			if(!match)continue;
			if(opts[k].addr==NULL)
			{
				ret=JZ_CMD_EMULTIPLE;
				goto __EXIT;
			}
			if(opt&JZ_CMD_PARAM)
			{
				if(param==NULL)
				{
					if(pos>=argc)
					{
						ret=JZ_CMD_EPARAM;
						goto __EXIT;
					}
					param=argv[pos++];
				}
			}
			if(jz_cmd_pass_param(opts[k].addr,opts[k].type,
				opt&JZ_CMD_MULTIPLE,param))
			{
				ret=JZ_CMD_ETYPE;
				goto __EXIT;
			}
			if(!(opt&JZ_CMD_MULTIPLE))opts[k].addr=NULL;
			if(opt&JZ_CMD_VAR)addvar=1;
			goto __DONE;
		}
		if(prefix_type==0)ret=JZ_CMD_ENUMBER;
		else
		{
			if(!(option&JZ_CMD_VAR))
			{
				ret=JZ_CMD_EARGUMENT;
				goto __EXIT;
			}
			else addvar=1;
		}
__DONE:;
		/* argument passed successfullly */
		if(addvar&&pairs)
		{
			p=strdup(p);
			if(param)q=strdup(param);else q=NULL;
			JZ_LIST_ADD(*strs,p);
			JZ_LIST_ADD(*strs,q);
		}
	}
__EXIT:	
	jz_cmd_arg=item;
	jz_cmd_param=param;
	return jz_cmd_err=ret;
}

/*-----------------------------------------------------------------------------
adjust data alignment when writing a file, at most 64 bytes. return 0 on 
success, or negative on failure. align is not required to be a power of two.
-----------------------------------------------------------------------------*/

int	jz_file_align(FILE *fout,int align)
{
	int		offset;
	size_t		size;
	char		dummy[64];

	offset=ftell(fout);
	if(offset<0)return -1;
	size=((offset-1)/align+1)*align-offset;
	if(size>0&&(fwrite(dummy,size,1,fout)!=1))return -2;
	return 0;
}

/*----------------------------------------------------------------------------
Read a line of unknown length. This is based on fgets, with additional
dirty work to allocate buffer of approriate size. If the length of the
string is less than or equal to size-1, the original buffer is not shrinked.
The '\n' at the end of the line is removed, and is not counted in length.
----------------------------------------------------------------------------*/

char	*jz_file_read_line(int *len,int size,FILE *stream)
{
	char		*ret,*buf;
	char		*str,**strs;
	int		n,m,count;

	buf=(char*)malloc(size);
	buf[0]=0;
	fgets(buf,size,stream);
	n=strlen(buf);
	if(n<=0)
	{
		*len=0;
		free(buf);
		return NULL;
	}
	if(buf[n-1]=='\n')
	{
		buf[n-1]=0;
		*len=n-1;
		return buf;
	}
	JZ_LIST_INIT(strs,8);
	m=size-1;
	for(;;)
	{
		str=(char*)malloc(size);
		str[0]=0;
		fgets(str,size,stream);
		n=strlen(str);
		if(n<=0)
		{
			free(str);
			break;
		}
		if(str[n-1]=='\n')str[--n]=0;
		m+=n;
		JZ_LIST_ADD(strs,str);
		if((n+1)<size)break;
	}
	size--;
	*len=m;
	ret=(char*)realloc(buf,m+1);
	count=JZ_LIST_COUNT(strs)-1;
	buf=ret+size;
	for(m=0;m<count;m++,buf+=size)
	{
		memcpy(buf,strs[m],size);
		free(strs[m]);
	}
	memcpy(buf,strs[m],n+1);
	free(strs[m]);
	JZ_LIST_FREE(strs);
	free(buf);
	return ret;
}

/*-----------------------------------------------------------------------------
Load a text file into memory, return a char** JZ_ARRAY of lines. must be 
freed with jz_str_lines_free(). This is faster than fgets() line by line
since the lines reside in the original buffer.
-----------------------------------------------------------------------------*/

char	**jz_file_read_lines(const char *fn,int *nline,int trim)
{
	int	fd;
	char	**ret;

	if(!fn)return NULL;
	fd=open(fn,O_RDONLY);
	if(fd<0)return NULL;
	ret=jz_file_read_lines_fd(fd,nline,trim);
	close(fd);
	return ret;
}

/*-----------------------------------------------------------------------------
Load a text file into memory, return a char** JZ_ARRAY of lines. must be 
freed with jz_str_lines_free(). This is faster than fgets() line by line
since the lines reside in the original buffer.
-----------------------------------------------------------------------------*/

char	**jz_file_read_lines_fd(int fd,int *nline,int trim)
{
	char		**ret=NULL;
	struct stat	st;
	size_t		size;
	char		*buf;
	
	/* load the text file */
	if(fstat(fd,&st))return NULL;
	size=st.st_size;

	if(size==0)
	{
		JZ_ARRAY_INIT(ret,1);
		if(ret==NULL)return NULL;
		*nline=0;
		ret[0]=NULL;
		return ret;
	}
	buf=(char*)malloc((size+6)&~3);
	if(!buf)return NULL;
	if(read(fd,buf,size)!=size)
	{
		free(buf);
		return NULL;
	}
	ret=jz_str_lines_link(buf,size,nline,trim);
	if(!ret)JZ_ARRAY_FREE(buf);
	return ret;
}

/*----------------------------------------------------------------------------
alternative to jz_file_read_lines. must be freed with jz_str_list_free().
----------------------------------------------------------------------------*/

char    **jz_file_read_list(char *fn)
{
	char	**ret;
	FILE	*file;
	int	len;
	char	*line,*str;

	if((file=fopen(fn,"rt"))==NULL)return NULL;
	JZ_LIST_INIT(ret,8);
	for(;;)
	{
		line=jz_file_read_line(&len,256,file);
		if(!line)break;
		if(len==0)
		{
			free(line);
			continue;
		}
		JZ_ARRAY_RESIZE(line,len+1,str);
		JZ_LIST_ADD(ret,line);
	}
	fclose(file);
	return ret;
}

/*----------------------------------------------------------------------------
alternative to jz_file_read_lines. must be freed with jz_str_list_free().
----------------------------------------------------------------------------*/

char    **jz_file_read_list_stream(FILE *file)
{
	char	**ret;
	int	len;
	char	*line,*str;

	JZ_LIST_INIT(ret,8);
	for(;;)
	{
		line=jz_file_read_line(&len,256,file);
		if(!line)break;
		if(len==0)
		{
			free(line);
			continue;
		}
		JZ_ARRAY_RESIZE(line,len+1,str);
		JZ_LIST_ADD(ret,line);
	}
	fclose(file);
	return ret;
}

/*----------------------------------------------------------------------------
search the input in the options list, return the index if a match is found.
return -1 if no found.
----------------------------------------------------------------------------*/

int	jz_str_choice(char **options,int n,char *input,int ncase)
{
	int	(*comp_str)(const char *,const char *);
	int	k;

	comp_str=ncase?&strcasecmp:&strcmp;
	for(k=0;k<n;k++)
		if(!comp_str(options[k],input))return k;
	return -1;	
}

void    jz_str_chomp(char* s)
{
	if(*s==0)return;
	for(;*s;s++);
	for(--s;(*s=='\n')||(*s=='\r');s--);
	s[1]=0;
}

/*----------------------------------------------------------------------------
input a chomped string from a stream
----------------------------------------------------------------------------*/

int	jz_str_input(char *s,int size,FILE *stream)
{
	int len;
	if(fgets(s,size,stream)==NULL)return -1;
	len=strlen(s);
	if((len>0)&&(s[len-1]=='\n'))s[--len]=0;
	return len;
}

/*-----------------------------------------------------------------------------
dispose a char** array returned by jz_str_lines_link.
-----------------------------------------------------------------------------*/

void	jz_str_lines_free(char **array)
{
	if(!array)return;
	if(array[0])free(array[0]);
	free(array);
}


/*-----------------------------------------------------------------------------
This function is to load a text file into memory with a single I/O call.

Given a piece of memory, create a char** JZ_ARRAY of lines, which is
ended with a NULL pointer. The array must be freed with jz_str_lines_free(). 
The number of lines is returned in *nline. An initial guess of nline is
passed through *nline. The original memory is modified to replace '\n'
with '\0''s at the end of each line. NULL return value indicates failure.
If trim is not zero, then remove spaces at both ends of each line and 
eliminate all empty lines.

On windows machines, EOL contains two characters \r\n. Please turn on
trim flag before it is redesigned for windows machine.
-----------------------------------------------------------------------------*/

#define	__IS_SPACE(ch)							\
	((ch)==' '||(ch)=='\t'||(ch)=='\v'||(ch)=='\r'||(ch)=='\f')
	
char	**jz_str_lines_link(char *mem,size_t size,int *nline,int trim)
{
	char		**ret=NULL,**line;
	int		capacity,count;
	char		*p,*q,*r,*end,ch;
	
	if(!mem)return NULL;

	/* Some files are not ended with a new line. Use an extra byte */
	end=mem+size;
	end[0]='\n';
	end[1]='Z';

	/* Initialize the returning line array */
	capacity=*nline;
	if(capacity<=0)capacity=256;
	count=0;
	JZ_ARRAY_INIT(ret,capacity);
	if(ret==NULL)return NULL;
	if(trim)goto __TRIM;

	/* build lines without trimming */
	p=mem;
	if(p>=end)goto __EXIT;
	for(;;)
	{
		/* use memchr in stead */
		/* for(q=p;*q!='\n';q++); */
		
		q=memchr(p,'\n',size);
		if(q>=end)break;
		*q=0;
		JZ_ARRAY_NEW(ret,line,capacity,count,2);
		*line=p;
		p=q+1;
	}
	goto __EXIT;
	
__TRIM:
	/* build trimed lines */
	/* let p points to the first non-space character */
	for(p=mem;;p++)
	{
		ch=*p;
		if(!(__IS_SPACE(ch)||ch=='\n'))break;
	}
	if(p>=end)goto __EXIT;
	q=memchr(p+1,'\n',size);
	/* for(q=p+1;*q!='\n';q++); */
	for(r=q-1;;r--)
	{
		ch=*r;
		if(!__IS_SPACE(ch))break;
	}
	r[1]=0;

	/* Move the first line to free the mem buffer later */
	if(p!=mem)
	{
		memmove(mem,p,r-p+1);
		p=mem;
	}

	for(;;)
	{	
		/* p points to the first non-space character;
		 * q points to the terminal '\n', add this line */
		JZ_ARRAY_NEW(ret,line,capacity,count,2);
		*line=p;
		if(q>=end)break;

		/* Find the first non-empty character */
		for(p=q+1;;p++)
		{
			ch=*p;
			if(!(__IS_SPACE(ch)||ch=='\n'))break;
		}
		if(p>=end)break;
		q=memchr(p+1,'\n',size);
		/* for(q=p+1;*q!='\n';q++); */
		for(r=q-1;;r--)
		{
			ch=*r;
			if(!__IS_SPACE(ch))break;
		}
		r[1]=0;
	}
	
__EXIT:
	*nline=count;
	JZ_ARRAY_NEW(ret,line,capacity,count,1.1);
	*line=NULL;
	JZ_ARRAY_RESIZE(ret,count+1,line);
	return ret;
}

#undef	__IS_SPACE


void	jz_str_list_free(char **str_list)
{	
	int     k,count;
	
	count=JZ_LIST_COUNT(str_list);
	for(k=0;k<count;k++)if(str_list[k])free(str_list[k]);
}

/*----------------------------------------------------------------------------
remove spaces on both sides
----------------------------------------------------------------------------*/

void	jz_str_trim(char *s)
{
	char *q;
	for(s=q=s;isspace(*q);q++);
	for(;*q;)*s++=*q++;
	for(--s;isspace(*s);)--s;
	*(++s)=0;
}


