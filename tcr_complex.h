#include <string>
#include <stdio.h>

using namespace std;

static const double PI = 3.14159;

struct Atom
{
  double x;
  double y;
  double z;
  int res;
  string res_name;
  string ins_code;
  string atom;
  string line;
};

// struct Protein: store an array of atoms
struct ProteinChain
{
  Atom* Atoms;
  int num_atoms;
  string chain_id;
};

class TCRComplex
{
 public:
  TCRComplex(string tcrfile, string pepmhcfile = "", int mhc_type = 0);
  ~TCRComplex(); // destructor
  double CalcDockingAngle();
  void SetXShift(double shift) { my_tcrx_shift = shift; }
  void SetZShift(double shift) { my_tcrz_shift = shift; }
  void SetTiltShift(double shift) { my_tcr_tilt_shift = shift; }
  void SetCrossShift(double shift) { my_tcr_cross_shift = shift; }

 private:
  void LoadChainFromFile(string, string, ProteinChain*);
  void LoadTCR(string); // load proteins from file
  void LoadPepMHC(string);
  bool IsAlreadyLoaded(ProteinChain*, int, string, string);
  bool GetCoords(ProteinChain*, int, string, string, double&, double&, double&);
  bool GetCoordsRes(ProteinChain*, int, string, string, string, double&, double&, double&);
  void CalcTCRSGVector(double&, double&, double&, double&, double&, double&);
  void CalcTCRRotVector(double&, double&, double&);
  void AlignTCRRotVector(double&, double&, double&, double&, double&, double&);
  void AlignTCRToAxis(FILE*, FILE*, double, double, double, double, double, double, double, double, double, double&, double&, double&, bool shift = false);
  void AlignMHCToAxis(FILE*, FILE*, FILE*, double, double, double, double, double, double, double, double, double, double&, double&, double&);
  void AlignMolToVector(double*, double*, double*, double, double, double, double&, double&, double&,
			double, double, double, int);
  void CalcMHCVectors(double&, double&, double&, double&, double&, double&, double&, double&, double&);
  void CalcCrossProd(double, double, double, double, double, double, double&, double&, double&);
  void OutputAxis(FILE*, double, double, double, double, double, double, string chain = "G");
  double GetRotationAxisAngle(float*, double&, double&, double&);

  string my_tcra_chain;
  string my_tcrb_chain;
  string my_mhc_chain;
  string my_mhcb_chain;
  string my_pep_chain;

  bool IsClassII;
  int my_mhc_type;

  bool verbose;

  string my_tcr_file;
  string my_pepmhc_file;

  double my_z_offset;
  double my_tcrx_shift;
  double my_tcrz_shift;
  double my_tcr_tilt_shift;
  double my_tcr_cross_shift;

  int tcra_cys_shift;
  int tcrb_cys_shift;

  ProteinChain my_mhc;
  ProteinChain my_mhcb;
  ProteinChain my_pep;
  ProteinChain my_tcra;
  ProteinChain my_tcrb;
};
