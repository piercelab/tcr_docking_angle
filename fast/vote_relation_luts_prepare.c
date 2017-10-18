/******************************************************************************
This file is part of the FAST protein structure alignment package,
developed by Jianhua Zhu at the Bioinformatics Program of Boston University.
******************************************************************************/

#include <stdio.h>
#include "vote.h"

/*-----------------------------------------------------------------------------
LUT tables to calculate 
-----------------------------------------------------------------------------*/

void	vote_relation_luts_prepare(FILE *fout)
{
	int		k,b;
	uint16_t	bit;
	uint16_t	bits[VC_RR_DISTANCE_MAX+128];
	uint16_t	masks[VC_RR_DISTANCE_MAX+128];
	uint16_t	angle_bits[1025];
	uint16_t	angle_masks[1025];
	uint16_t	angle_masks2[1025];

	/* Distance LUT and MASK */

	bit=1;
	b=128;
	
	for(k=0;k<VC_RR_DISTANCE_MAX+128;k++)
	{
		if(k>=b)
		{
			bit<<=1;
			if(bit==0)break;
			b=(int)(b*1.20)+256;
		}
		bits[k]=bit;
		masks[k]=bit|(bit<<1)|(bit>>1);
	}

	fprintf(fout,"#include \"vote.h\"\n\n");

	fprintf(fout,"uint16_t\tvote_distance_bit_lut[]=\n{");
	for(k=0;k<VC_RR_DISTANCE_MAX+128;k++)
	{
		if(k&3)continue;
		if((k&31)==0)fprintf(fout,"\n");
		fprintf(fout,"\t%u,",(unsigned int)bits[k]);
	}
	fprintf(fout,"\n};\n\n");

	fprintf(fout,"uint16_t\tvote_distance_mask_lut[]=\n{");
	for(k=0;k<VC_RR_DISTANCE_MAX+128;k++)
	{
		if(k&3)continue;
		if((k&31)==0)fprintf(fout,"\n");
		fprintf(fout,"\t%u,",(unsigned int)masks[k]);
	}
	fprintf(fout,"\n};\n\n");

	/* Angle LUTS */
	b=64;
	bit=1;
	for(k=0;k<1025;k++)
	{
		if(k>=b)
		{
			bit<<=1;
			if(bit==0)break;
			b+=64;
		}
		angle_bits[k]=bit;
		angle_masks[k]=bit|(bit<<1)|(bit<<2)|(bit<<3)|(bit<<4)
			|(bit>>1)|(bit>>2)|(bit>>3)|(bit>>4);
		angle_masks2[k]=bit|(bit<<1)|(bit<<2)|(bit<<3)|(bit<<4)
			|(bit<<5)|(bit<<6)|(bit>>1)|(bit>>2)|(bit>>3)|(bit>>4)
			|(bit>>5)|(bit>>6);
	}

	angle_bits[1024]=angle_bits[1023];
	angle_masks[1024]=angle_masks[1023];
	angle_masks2[1024]=angle_masks2[1023];

	fprintf(fout,"uint16_t\tvote_angle_bit_lut[]=\n{");
	for(k=0;k<1025;k++)
	{
		if((k&7)==0)fprintf(fout,"\n");
		fprintf(fout,"\t%u,",(unsigned int)angle_bits[k]);
	}
	fprintf(fout,"\n};\n\n");

	fprintf(fout,"uint16_t\tvote_angle_mask_lut[]=\n{");
	for(k=0;k<1025;k++)
	{
		if((k&7)==0)fprintf(fout,"\n");
		fprintf(fout,"\t%u,",(unsigned int)angle_masks[k]);
	}
	fprintf(fout,"\n};\n\n");

	fprintf(fout,"uint16_t\tvote_angle_mask_lut2[]=\n{");
	for(k=0;k<1025;k++)
	{
		if((k&7)==0)fprintf(fout,"\n");
		fprintf(fout,"\t%u,",(unsigned int)angle_masks2[k]);
	}
	fprintf(fout,"\n};\n\n");

}


int	main(int argc,char **argv)
{
	FILE		*fout;

	fout=fopen(argv[1],"w");
	if(!fout)exit(1);
	vote_relation_luts_prepare(fout);
	fclose(fout);
	return 0;
}


