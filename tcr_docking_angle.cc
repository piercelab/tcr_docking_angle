#include <iostream>
#include "tcr_complex.h"

int main(int argc, char** argv)
{ 
  string tcrfile = "";
  string pepmhcfile = "";
  int mhc_type = 0;
  double xshift = 0.0;
  double crossshift = 0.0;
  double angleshift = 0.0;

  if (argc < 2)
    {
      cout << "usage:\n" <<
	"tcr_docking_angle tcr_complex_file" << endl << "OR" << endl <<
	"tcr_docking_angle tcr_file mhcpep_file" << endl << "OR" << endl <<
	"tcr_docking_angle tcr_complex_file mhc_type" << endl << "OR" << endl <<
	"tcr_docking_angle tcr_file mhcpep_file mhc_type" << endl <<
	"MHC Type: 0 = Class I; 1 = Class II; 2 = CD1d, CD1b; 3 = MR1; 4 = CD1d (4EI5); 5 = CD1c" << endl;
      exit(0);
    }
  else if (argc == 2)
    {
      tcrfile = argv[1];
    }
  else if (argc == 3)
    {
      tcrfile = argv[1];
      pepmhcfile = argv[2];
      if (pepmhcfile.length() == 1) 
	{
	  mhc_type = atoi(argv[2]);
	  pepmhcfile = "";
	}
    }
  else if (argc == 4)
    {
      tcrfile = argv[1];
      pepmhcfile = argv[2];
      mhc_type = atoi(argv[3]);
    }
  else if (argc == 7)
    {
      tcrfile = argv[1];
      pepmhcfile = argv[2];
      mhc_type = atoi(argv[3]);
      xshift = atof(argv[4]);
      crossshift = atof(argv[5]);
      angleshift = atof(argv[6]);
    }
  
  TCRComplex comp(tcrfile, pepmhcfile, mhc_type);
  if (argc == 7)
    {
      comp.SetXShift(xshift);
      comp.SetCrossShift(crossshift);
      comp.SetTiltShift(angleshift);
    }
  //comp.SetZShift(8.0); // trying this for start3 position
  comp.CalcDockingAngle();
  
  return 0;
}
