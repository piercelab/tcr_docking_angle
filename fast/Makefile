SHELL=/bin/sh
BINDIR=${HOME}/bin
CC=gcc -O3 -Wall -I.

CCO=${CC} -c -o $@ $<

.SUFFIXES:.cc .c .cpp
.cc.o:
	${CCO}
.c.o:
	${CCO}
.cpp.o:
	${CCO}

OBJS=	fast.o Alignment.o rasmol.o vote_pairwise.o vote_discretize.o\
	vote_lgc_dist.o vote_relation.o vote_align.o vote_lgc.o\
	vote_trim_lgc.o vote_candidate.o vote_tally.o vote_eliminate.o\
	vote_alignment.o vote_shrink.o vote_refine.o vote_decay_lut.o\
	vote_matrix.o vote_score.o vote_trim.o vote_acos_lut.o\
	vote_follow.o vote_remove.o vote_rms_fit.o vote_extend.o\
	vote_relation_luts.o vote_consolidate.o basic.o chore.o misc.o\
	vote_tally_fast.o vote_consolidate2.o vote_cache.o\
	jz_amino.o jz_protein.o vote_rmsd.o

LIBOBJS=	Alignment.o rasmol.o vote_pairwise.o vote_discretize.o\
	vote_lgc_dist.o vote_relation.o vote_align.o vote_lgc.o\
	vote_trim_lgc.o vote_candidate.o vote_tally.o vote_eliminate.o\
	vote_alignment.o vote_shrink.o vote_refine.o vote_decay_lut.o\
	vote_matrix.o vote_score.o vote_trim.o vote_acos_lut.o\
	vote_follow.o vote_remove.o vote_rms_fit.o vote_extend.o\
	vote_relation_luts.o vote_consolidate.o basic.o chore.o misc.o\
	vote_tally_fast.o vote_consolidate2.o vote_cache.o\
	jz_amino.o jz_protein.o vote_rmsd.o

libs:	${LIBOBJS}
	ar r libfast.dylib ${LIBOBJS} # -lm

fast:	${OBJS}
	${CC} -o fast ${OBJS} -lm # -static


clean:
	/bin/rm -f core a.out *.o fast vote_relation_luts_prepare

install:
	/bin/cp -f fast ${BINDIR}


## LUTS ARE PREPARED HERE

vote_relation_luts_prepare:vote_relation_luts_prepare.o
	${CC} -o vote_relation_luts_prepare vote_relation_luts_prepare.o 

vote_relation_luts.c:vote_relation_luts_prepare vote.h
	./vote_relation_luts_prepare vote_relation_luts.c

## END OF LUTS



