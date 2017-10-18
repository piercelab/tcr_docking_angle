/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

/*=============================================================================
chore.h		the command line parser and assorted small utilities
=============================================================================*/

#ifndef	__JZ_CHORE_H
#define	__JZ_CHORE_H

#ifdef	__cplusplus

extern	"C" {

#endif

/*=============================================================================
jz_cmd			unix-style command line parser.

each command line item is specified by a jz_cmd record. all items are 
described by an jz_cmd array, in which order determined items (those 
without -x or --xxx prefixes) are placed in the order desired, before 
any items with a -x or --xxx prefixes. 

jz_cmd.option is a integer containing the following information, used in 
combination (or):

name		flag	description
-----------------------------------------------------------------------
JZ_CMD_PREFIX	1	option has a argument (-x or --xxx)
JZ_CMD_PARAM	2	option has a required parameter 
JZ_CMD_DEFAULT	4	argument has a default param
JZ_CMD_MULTIPLE	8	argument can be used multiple times ==> jz_list
JZ_CMD_UNIQUE	16	parameter must be unique if multiple

jz_cmd.prefix is a string containing all possible prefixes, seperated 
by |, for example: "i|g|in|input" means that short prefixes -i, -g and 
long prefixes --in, --input are regarded as the same option. 

jz_cmd.type is the data type of the parameter. the parser automatically 
check if the user supplied legal parameter of the given type. 
jz_cmd.address is the address of the variable which will be assigned the 
value following the prefix. if the parameter is of basic data type like 
JZ_TFLOAT, the variable must be defined as the same type. if the argument 
can be used multiple times, the variable must be defined as type * and 
initialized to NULL (or the user may initialize the list himself);
 
the behavior of the parser itself can be modified by a list of parameters 
in the calling variable list. NOTE: not implemented.

name		param	description
---------------------------------------------------------------------------
JZ_CMD_SEPERATE	1	allow the standard form -o <param>
JZ_CMD_ASSIGN	2	allow -o=<param> and --option=<param>
JZ_CMD_LONGSHORT	4	allow long prefixes to use short letter
JZ_CMD_SHORTLONG	8	allow short prefixes to use long letter
JZ_CMD_COMBINE	16	allow combination of short flags
JZ_CMD_FOLLOW	32	allow -o<param> for short prefixes
---------------------------------------------------------------------------

the parser return 0 on success. once there is an error, the parser stops
immediately and return a non-zero to signify an error. the returned param
can be passed to the dignostic routine jz_cmd_message() to print an error 
message.

return name		param	description
---------------------------------------------------------------------------
JZ_CMD_SUCCESS		0	no error
JZ_CMD_ENUMBER		1	number of arguments is incorrect
JZ_CMD_EPREFIX		2	prefix is not defined
JZ_CMD_EDEFAULT		3	an argument without default is missing
JZ_CMD_ETYPE		4	type of a supplied argument is incorrect
JZ_CMD_EMULTIPLE	5	option used multiple times but not allowed
---------------------------------------------------------------------------

=============================================================================*/

#define	JZ_SHORT_PREFIX			'-'
#define	JZ_LONG_PREFIX			"--"
#define	JZ_LONG_PREFIX_LENGTH		2

#define	JZ_CMD_ARGUMENT			1
#define	JZ_CMD_PARAM			2
#define	JZ_CMD_MULTIPLE			4
#define	JZ_CMD_REQUIRED			8
#define	JZ_CMD_VAR			16
#define	JZ_CMD_INDEX			32

#define	JZ_CMD_SEPERATE			1
#define	JZ_CMD_ASSIGN			2
#define	JZ_CMD_LONGSHORT		4
#define	JZ_CMD_SHORTLONG		8
#define	JZ_CMD_COMBINE			16
#define	JZ_CMD_FOLLOW			32

#define	JZ_CMD_SUCCESS			0
#define	JZ_CMD_ENUMBER			1
#define	JZ_CMD_EARGUMENT		2
#define	JZ_CMD_EPARAM			3
#define	JZ_CMD_ETYPE			4
#define	JZ_CMD_EMULTIPLE		5

#define	JZ_TNONE			0
#define	JZ_TBOOL			1
#define	JZ_TCHAR			2
#define	JZ_TINT				3
#define	JZ_TFLOAT			4
#define	JZ_TDOUBLE			5
#define	JZ_TSTR				6
#define	JZ_TMESSAGE			7
#define	JZ_TFUNC			9

typedef	struct	jz_cmd
{
	char	option;
	char	type;
	void	*addr;
	char	*arg;
}jz_cmd;

int	jz_cmd_parse(int argc,char **argv,jz_cmd *opts,int opt_count,
	int option,jz_pair **kvpairs,char ***strs);
	
void	jz_cmd_message(FILE *stream);

#define	JZ_CMD_USED(opt)		(!((opt).addr))

#define	JZ_CMD_PARSE(opts)						\
	jz_cmd_parse(argc,argv,(opts),sizeof(opts)/sizeof(jz_cmd),\
	0,NULL,NULL)

#define	JZ_CMD_PARSE_UNIX		JZ_CMD_PARSE	

#define	JZ_CMD_PARSE_VAR(opts,kvpairs,strs)				\
	jz_cmd_parse(argc,argv,(opts),sizeof(opts)/sizeof(jz_cmd),\
	JZ_CMD_VAR,&(kvpairs),&(strs))

extern	int jz_cmd_err;
extern	char *jz_cmd_arg;
extern	char *jz_cmd_param;


#define	JZ_ERR_HERE(msg)		fprintf(stderr,msg);

#define	JZ_ERROR_HERE			JZ_ERR_HERE

#define	JZ_ERROR_GOTO(expression,__ERR_EXIT)				\
	do\
	{\
		(expression);\
		goto __ERR_EXIT;\
	}while(0)

#define	JZ_ERR_GOTO			JZ_ERROR_GOTO	

#define	JZ_ERROR_EXIT(stream,code,message)				\
	do\
	{\
		fprintf((stream),(message));\
		exit(code);\
	}while(0)

void	jz_err_exit(FILE *stream,int err,char *s);


int	jz_file_align(FILE *fout,int align);

char*	jz_file_read_line(int *len,int size,FILE *stream);

char**	jz_file_read_lines_fd(int fd,int *nline,int trim);

char**	jz_file_read_lines(const char *fn,int *nline,int trim);

char**	jz_file_read_list(char *fn);

char**	jz_file_read_list_stream(FILE *stream);


#define	JZ_STR_CHOICE(haystack,needle,ncase)				\
	jz_str_choice((haystack),sizeof(haystack)/sizeof((haystack)[0]),\
	(needle),(ncase))

void    jz_str_chomp(char* s);

void	jz_str_trim(char *s);

void	jz_str_toupper(char *s);

void	jz_str_tolower(char *s);

int	jz_str_input(char *s,int size,FILE *stream);

void	jz_str_list_free(char **str_list);

void	jz_str_lines_free(char **lines);

char**	jz_str_lines_link(char *mem,size_t size,int *nline,int trim);

int	jz_str_choice(char **options,int n,char *input,int ncase);

#ifdef	__cplusplus

} /* extern "C" */

#endif

#endif

