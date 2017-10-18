/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include "vote.h"

/*-----------------------------------------------------------------------------
Consolidate the changes
-----------------------------------------------------------------------------*/

void	vote_consolidate2(Voter **VC,int nvoter)
{
	VoteEdge	**addr,*p;
	int		score,m,k;
	
	/* consolidate the changes */
	for(k=0;k<nvoter;k++)
	{
		if(VC[k]->nvote<0)continue;

		m=VC[k]->nvote;
		score=VC[k]->score;
		
		/* follow the list of supporters, remove those eliminated */
		for(addr=&(VC[k]->link);(p=*addr);)
		{
			if(p->source->nvote<0)
			{
				*addr=p->next;
				m--;
				score-=p->score;
			}	
			else addr=&(p->next);
		}	
		
		VC[k]->nvote=m;
		VC[k]->score=score;
	}
}


