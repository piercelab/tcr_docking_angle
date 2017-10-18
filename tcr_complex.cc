#include "tcr_complex.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include "fast/jz_protein.h"
#include "fast/Alignment.h"
#include "fast/vote.h"


TCRComplex::TCRComplex(string tcrfile, string pepmhcfile, int mhc_type)
{
  // initialization
  my_tcra_chain = "D";
  my_tcrb_chain = "E";
  my_mhc_chain = "A";
  my_pep_chain = "C";
  
  
  my_mhc_type = mhc_type;
  if (my_mhc_type == 1) IsClassII = true;

  my_z_offset = 25.0;
  my_tcrz_shift = 0.0;  

  // alternative TCR docking start site
  // start2
  my_tcrx_shift = 20.0;
  my_tcr_tilt_shift = 25.0;
  my_tcr_cross_shift = -45.0;

  /*
  // start3
  my_tcrx_shift = 10.0;
  my_tcr_tilt_shift = 12.5;
  my_tcr_cross_shift = -22.5;
  */

  /*
  // start4
  my_tcrx_shift = 25.0;
  my_tcr_tilt_shift = 40.0;
  my_tcr_cross_shift = -60.0;
  */

  if (pepmhcfile == "") pepmhcfile = tcrfile;
  if (IsClassII) my_mhcb_chain = "B";

  my_tcr_file = tcrfile;
  my_pepmhc_file = pepmhcfile;

  verbose = true;

  LoadTCR(my_tcr_file);
  LoadPepMHC(my_pepmhc_file);
}

TCRComplex::~TCRComplex()
{
  delete[] my_tcra.Atoms;
  delete[] my_tcrb.Atoms;
  delete[] my_mhc.Atoms;
  delete[] my_pep.Atoms;
  if (IsClassII) delete[] my_mhcb.Atoms;
}

void TCRComplex::LoadPepMHC(string filename)
{
  LoadChainFromFile(filename, my_mhc_chain, &my_mhc);
  // cout << my_mhc.num_atoms << " atoms loaded for MHC" << endl;
  LoadChainFromFile(filename, my_pep_chain, &my_pep);
  // cout << my_pep.num_atoms << " atoms loaded for peptide" << endl;
  if (my_mhc.num_atoms == 0)
    {
      cerr << "No atoms loaded for MHC!" << endl;
      exit(1);
    }
  if (IsClassII)
    {
      LoadChainFromFile(filename, my_mhcb_chain, &my_mhcb);
      if (my_mhcb.num_atoms == 0)
	{
	  cerr << "No atoms loaded for MHC B chain!" << endl;
	  exit(1);
	}
    }
  
}

void TCRComplex::LoadTCR(string filename)
{
  LoadChainFromFile(filename, my_tcra_chain, &my_tcra);
  LoadChainFromFile(filename, my_tcrb_chain, &my_tcrb);
  // cout << my_tcra.num_atoms << " atoms loaded for TCR alpha" << endl;
  // cout << my_tcrb.num_atoms << " atoms loaded for TCR beta" << endl;
  if ((my_tcra.num_atoms == 0) || (my_tcrb.num_atoms == 0)) 
    {
      cerr << "No atoms loaded for TCR a or b!" << endl;
      exit(1);
    }
}

void TCRComplex::LoadChainFromFile(string filename, string chain_id, ProteinChain* prot)
{
  string tmpbuf;
  int num_atoms = 0, n_atm = 0;

  ifstream pdbfile(filename.c_str());
  if (! pdbfile.is_open())
    {
      cerr << "Error opening pdb file: " << filename << endl;
      exit(1);
    }

  // count the number of atoms so we can allocate the array
  while (getline(pdbfile, tmpbuf))
    {
      if ((tmpbuf.substr(0, 4) != "ATOM") && (tmpbuf.substr(0, 6) != "HETATM")) continue;
      
      string ch_id = tmpbuf.substr(21, 1);
      if (ch_id == chain_id) num_atoms++;
    }
  pdbfile.close();

  prot->Atoms = new Atom[num_atoms];
  prot->chain_id = chain_id;
  prot->num_atoms = 0;
  

  // parse the atoms from the PDB file

  ifstream pdbfile2(filename.c_str());
  if (! pdbfile2.is_open())
    {
      cerr << "Error opening pdb file: " << filename << endl;
      exit(1);
    }

  // count the number of atoms so we can allocate the arrays
  while (getline(pdbfile2, tmpbuf))
    {
      if ((tmpbuf.substr(0, 4) != "ATOM") && (tmpbuf.substr(0, 6) != "HETATM")) continue;
      
      string ch_id = tmpbuf.substr(21, 1);
      if (ch_id == chain_id)
	{
	  double xcoord = atof(tmpbuf.substr(30, 8).c_str());
	  double ycoord = atof(tmpbuf.substr(38, 8).c_str());
	  double zcoord = atof(tmpbuf.substr(46, 8).c_str());
	  string res = tmpbuf.substr(17, 3);
	  int res_num = atoi(tmpbuf.substr(22, 4).c_str());
	  string ins_code = tmpbuf.substr(26, 1);
	  string atom = tmpbuf.substr(12, 4);
	  string res_name = tmpbuf.substr(17, 3);
	  
	  if (IsAlreadyLoaded(prot, res_num, ins_code, atom))
	    {
	      cout << "Warning: duplicate entries found in chain " << ch_id << " res " << res_num << ins_code << " atom " << atom << endl;
	      // exit(1);
	    }
	  prot->Atoms[prot->num_atoms].x = xcoord;
	  prot->Atoms[prot->num_atoms].y = ycoord;
	  prot->Atoms[prot->num_atoms].z = zcoord;
	  prot->Atoms[prot->num_atoms].res = res_num;
	  prot->Atoms[prot->num_atoms].ins_code = ins_code;
	  prot->Atoms[prot->num_atoms].atom = atom;
	  prot->Atoms[prot->num_atoms].res_name = res_name;
	  prot->Atoms[prot->num_atoms].line = tmpbuf;
	  prot->num_atoms++;
	}
    }
  pdbfile2.close();
  
}

// determine whether the protein atom has already been loaded into the array (hopefully not!)
bool TCRComplex::IsAlreadyLoaded(ProteinChain *prot, int res_num, string ins_code, string atom)
{
  for (int i = 0; i < prot->num_atoms; i++)
    {
      if ((prot->Atoms[i].res == res_num) && (prot->Atoms[i].ins_code == ins_code) && (prot->Atoms[i].atom == atom)) return true;
    }
  return false;
}

// determine whether the protein atom has already been loaded into the array (hopefully not!)
bool TCRComplex::GetCoords(ProteinChain *prot, int res_num, string ins_code, string atom, double &x, double &y, double &z)
{
  for (int i = 0; i < prot->num_atoms; i++)
    {
      if ((prot->Atoms[i].res == res_num) && (prot->Atoms[i].ins_code == ins_code) && (prot->Atoms[i].atom == atom))
	{
	  x = prot->Atoms[i].x;
	  y = prot->Atoms[i].y;
	  z = prot->Atoms[i].z;
	  return true;
	}
    }
  return false;
}

// return the coordinates of the specified residue number, attempting to match the residue names as well
bool TCRComplex::GetCoordsRes(ProteinChain *prot, int res_num, string ins_code, string atom, string res_name, double &x, double &y, double &z)
{
  for (int i = 0; i < prot->num_atoms; i++)
    {
      if ((prot->Atoms[i].res == res_num) && (prot->Atoms[i].ins_code == ins_code) && (prot->Atoms[i].atom == atom) && (prot->Atoms[i].res_name == res_name))
      {
	  x = prot->Atoms[i].x;
	  y = prot->Atoms[i].y;
	  z = prot->Atoms[i].z;
	  return true;
	}
    }
  return false;
}

double TCRComplex::CalcDockingAngle()
{
  double tcrx = 0.0, tcry = 0.0, tcrz = 0.0, tcrcentx = 0.0, tcrcenty = 0.0, tcrcentz = 0.0, tcrrotx = 0.0, tcrroty = 0.0, tcrrotz = 0.0;
  double mhcx = 0.0, mhcy = 0.0, mhcz = 0.0, mhccentx = 0.0, mhccenty = 0.0, mhccentz = 0.0, mhcnormx = 0.0, mhcnormy = 0.0, mhcnormz = 0.0;
  FILE *outfile, *outfile2, *outfile3, *outfile4, *outfile5, *outfile6, *dummyfile;

  CalcTCRSGVector(tcrx, tcry, tcrz, tcrcentx, tcrcenty, tcrcentz);
  CalcTCRRotVector(tcrrotx, tcrroty, tcrrotz);
  CalcMHCVectors(mhcx, mhcy, mhcz, mhccentx, mhccenty, mhccentz, mhcnormx, mhcnormy, mhcnormz);
  AlignTCRRotVector(tcrrotx, tcrroty, tcrrotz, tcrcentx, tcrcenty, tcrcentz); // negate the vector if necessary so that it's a normal vector
  
  double angle = acos(tcrx*mhcx + tcry*mhcy + tcrz*mhcz)*180.0/PI;
  //if (angle > 90.0) angle = 180 - angle;

  // get the x and y distance from the TCR rot vector intersection with MHC norm plane and z axis (assumes pre-aligned MHC
  // and TCR to the z axis)
  // 7.5 A from MHC cent seems to be a good approximation of the TCR binding site
  double line_dist = (my_z_offset - 7.5 - tcrcentz)/tcrrotz;
  double x_off = tcrrotx*line_dist + tcrcentx;
  double y_off = tcrroty*line_dist + tcrcenty;

  // cout << "Docking angle = " << angle << endl;
  
  // calculate distance between MHC normal vector and the CM point of the TCR
  double tmpx, tmpy, tmpz;
  CalcCrossProd(tcrcentx - mhccentx, tcrcenty - mhccenty, tcrcentz - mhccentz, 
		tcrcentx - mhccentx - mhcnormx, tcrcenty - mhccenty - mhcnormy, tcrcentz - mhccentz - mhcnormz,
		tmpx, tmpy, tmpz);

   // no need to divide by mag of the vector between points because that vector is normalized
  //cout << "distance = " << mag1 << endl;
  double normal_angle = acos(-tcrrotx*mhcnormx + -tcrroty*mhcnormy + -tcrrotz*mhcnormz)*180.0/PI;
  cout << tcrrotx << "\t" << tcrroty << "\t" << tcrrotz << "\t" << mhcnormx << "\t" << mhcnormy << "\t" << mhcnormz << endl;
  cout << "offsets: " << x_off << "\t" << y_off << endl;
	// get the vector to rotate between normal angles
  CalcCrossProd(-tcrrotx, -tcrroty, -tcrrotz, 
		mhcnormx, mhcnormy, mhcnormz,
		tmpx, tmpy, tmpz);
  
  double mag = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
  tmpx /= mag;
  tmpy /= mag;
  tmpz /= mag;

  // output the MHC normal axis to a PDB file  
  outfile3 = fopen("axis.pdb", "w");
  //OutputAxis(outfile3, mhccentx, mhccenty, mhccentz, mhcx, mhcy, mhcz);
  //OutputAxis(outfile3, mhccentx, mhccenty, mhccentz, mhcnormx, mhcnormy, mhcnormz);
  OutputAxis(outfile3, tcrcentx, tcrcenty, tcrcentz, tcrx, tcry, tcrz);
  OutputAxis(outfile3, tcrcentx, tcrcenty, tcrcentz, tcrrotx, tcrroty, tcrrotz, "H");
  fclose(outfile3);
 
  // align TCR to the x and z axes 
  outfile = fopen("rottcrmhc.tcrfixed.pdb", "w");
  outfile2 = fopen("rottcrmhc.pdb", "w");
  outfile4 = fopen("rottcrmhc.mhcfixed.pdb", "w");
  outfile6 = fopen("rottcrmhc2.pdb", "w");
  dummyfile = fopen("dummy.pdb", "w");

  double newtcrcentx = tcrcentx;
  double newtcrcenty = tcrcenty;
  double newtcrcentz = tcrcentz;
  double newmhccentx = mhccentx;
  double newmhccenty = mhccenty;
  double newmhccentz = mhccentz;
  double temp1, temp2, temp3;

  AlignTCRToAxis(outfile, outfile2, tcrrotx, tcrroty, tcrrotz, tcrx, tcry, tcrz, tcrcentx, tcrcenty, tcrcentz, newmhccentx, newmhccenty, newmhccentz, false);  
  // shift the TCR to a second starting point
  AlignTCRToAxis(dummyfile, outfile6, tcrrotx, tcrroty, tcrrotz, tcrx, tcry, tcrz, tcrcentx, tcrcenty, tcrcentz, temp1, temp2, temp3, true);
  fprintf(outfile2, "%s", "TER\n");
  fprintf(outfile6, "%s", "TER\n");
  AlignMHCToAxis(outfile2, outfile4, outfile6,  mhcnormx, mhcnormy, mhcnormz, mhcx, mhcy, mhcz, mhccentx, mhccenty, mhccentz, newtcrcentx, newtcrcenty, newtcrcentz);
  fclose(outfile4); 

  fclose(outfile);
  fclose(outfile2);
  fclose(outfile6);
  fclose(dummyfile);
  
  double z_diff = -newtcrcentz + my_z_offset;
  double tot_diff = sqrt(newtcrcentx*newtcrcentx + newtcrcenty*newtcrcenty + (newtcrcentz - my_z_offset)*(newtcrcentz - my_z_offset));

  // output the docking and normal angles
  cout << "ANGLES:" << "\t" << angle << "\t" << normal_angle << "\t" << z_diff << endl;
 
  // output the new tcr center as a HETATM in a PDB file
  outfile5 = fopen("tcrcent.pdb", "w");
  string tmpbuf = "HETATM    1  CA  CEN G   1       9.419  10.702  55.946  1.00 48.05           C  ";
  fprintf(outfile5, "%s%8.3f%8.3f%8.3f%s%s", tmpbuf.substr(0, 30).c_str(), newtcrcentx, newtcrcenty, newtcrcentz, tmpbuf.substr(54).c_str(), "\n");
  fclose(outfile5);
  
  return angle;
}

void TCRComplex::CalcTCRSGVector(double & tcrx, double & tcry, double & tcrz, double & tcrcentx, double & tcrcenty, double & tcrcentz)
{
  // calculate the vector between the TCR disulfides, the CM of that vector, and the angle of rotation between the Vb domains
  double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
  double tmp_x, tmp_y, tmp_z;
  if (! GetCoords(&my_tcra, 22, " ", " SG ", tmp_x, tmp_y, tmp_z))
    {
      if (verbose) cout << "error getting coordinate for SG of TCRa residue 22!" << endl;
      bool atom_found = false;
      for (int i = -4; i <= 4; i++)
	{
	  if (GetCoords(&my_tcra, 22 + i, " ", " SG ", tmp_x, tmp_y, tmp_z))
	    {
	      if (verbose) cout << "using residue " << 22 + i << " instead" << endl;
	      atom_found = true;
	      break;
	    }
	}
      if (! atom_found) exit(0);
    }
  x1 += tmp_x;
  y1 += tmp_y;
  z1 += tmp_z;
  
  if (! GetCoords(&my_tcra, 90, " ", " SG ", tmp_x, tmp_y, tmp_z))
    {
      if (verbose) cout << "error getting coordinate for SG of TCRa residue 90!" << endl;
      bool atom_found = false;
      for (int i = -6; i <= 15; i++) // going up to residue 104 for 3DXA
	{
	  if (GetCoords(&my_tcra, 90 + i, " ", " SG ", tmp_x, tmp_y, tmp_z))
	    {
	      if (verbose) cout << "using residue " << 90 + i << " instead" << endl;
	      atom_found = true;
	      tcra_cys_shift = i;
	      break;
	    }
	}
      if (! atom_found) exit(0);
    }
  x1 += tmp_x;
  y1 += tmp_y;
  z1 += tmp_z;

  if (! GetCoords(&my_tcrb, 23, " ", " SG ", tmp_x, tmp_y, tmp_z))
    {
      if (verbose) cout << "error getting coordinate for SG of TCRb residue 23!" << endl;
      bool atom_found = false;
      for (int i = -4; i <= 4; i++)
	{
	  if (GetCoords(&my_tcrb, 23 + i, " ", " SG ", tmp_x, tmp_y, tmp_z))
	    {
	      if (verbose) cout << "using residue " << 23 + i << " instead" << endl;
	      atom_found = true;
	      break;
	    }
	}
      if (! atom_found) exit(0);
    }
  x2 += tmp_x;
  y2 += tmp_y;
  z2 += tmp_z;
  
  if (! GetCoords(&my_tcrb, 92, " ", " SG ", tmp_x, tmp_y, tmp_z))
    {
      if (verbose) cout << "error getting coordinate for SG of TCRb residue 92!" << endl;
      bool atom_found = false;
      for (int i = -6; i <= 15; i++)  // going up to residue 104 for 3DXA
	{
	  if (GetCoords(&my_tcrb, 92 + i, " ", " SG ", tmp_x, tmp_y, tmp_z))
	    {
	      if (verbose) cout << "using residue " << 92 + i << " instead" << endl;
	      atom_found = true;
	      tcrb_cys_shift = i;
	      break;
	    }
	}
      if (! atom_found) exit(0);
    }
  x2 += tmp_x;
  y2 += tmp_y;
  z2 += tmp_z;

  // get the average 3D coordinate of each pair of disulfide SG atoms
  x1 /= 2.0;
  y1 /= 2.0;
  z1 /= 2.0;
  x2 /= 2.0;
  y2 /= 2.0;
  z2 /= 2.0;

  // calculate the vector from tcra to tcrb disulfides
  tcrx = x2 - x1;
  tcry = y2 - y1;
  tcrz = z2 - z1;
  
  tcrcentx = (x1 + x2)/2.0;
  tcrcenty = (y1 + y2)/2.0;
  tcrcentz = (z1 + z2)/2.0;
  
  // normalize the vector
  double mag = sqrt(tcrx*tcrx + tcry*tcry + tcrz*tcrz);
  tcrx /= mag;
  tcry /= mag;
  tcrz /= mag;

  // cout << "TCR vector: (" << tcrx << "," << tcry << "," << tcrz << ")" << endl;

  return;
}

void TCRComplex::CalcTCRRotVector(double & tcrrotx, double & tcrroty, double & tcrrotz)
{
  jz_protein protein1, *p1;
  jz_protein protein2, *p2;
  jz_pair_int *pairs;
  
  double rot_angle;
  
  float ori[3], rot[9];

  Alignment alignment;
  
  alignment_create(&alignment);
  jz_protein_null(&protein1);
  jz_protein_null(&protein2);

  // will just use the entire input chains, truncate externally
  bp_protein_load(&protein1,my_tcr_file.c_str(),255,'D',-1,200);
  bp_protein_load(&protein2,my_tcr_file.c_str(),255,'E',-1,200);
  p1 = &protein1;
  p2 = &protein2;

  // compute the alignment
  vote_pairwise(p1,p2,&alignment);
  
  int aln_len = alignment.length;

  if (aln_len == 0) 
    {
      cerr << "Error: zero length alignment for TCR chains!" << endl;
      exit(1);
    }
  
  JZ_ARRAY_INIT(pairs, aln_len);
  for (int i = 0; i < aln_len; i++)
    {
      pairs[i][0] = alignment.pairs[i][0];
      pairs[i][1] = alignment.pairs[i][1];
    }
  
  // get the translation and rotation matrix for the alignment
  vote_rms_fit(pairs, aln_len, &protein1, &protein2, ori, rot);

  rot_angle = GetRotationAxisAngle(rot, tcrrotx, tcrroty, tcrrotz);

  if (verbose) cout << "TCR alignment length = " << aln_len << " RMSD = " << alignment.rmsd << " angle = " << rot_angle << endl;

  jz_protein_clear(&protein1);
  jz_protein_clear(&protein2);
  
  return;
}

void TCRComplex::AlignTCRRotVector(double & tcrrotx, double & tcrroty, double & tcrrotz, double & tcrcentx, double & tcrcenty, double & tcrcentz)
{
  double tmp_x, tmp_y, tmp_z;
  // make the vector point outward from the CDR surface
  // trying the GLN residues around residue 37 on each chain
  double cterm_cent_x = 0.0;
  double cterm_cent_y = 0.0;
  double cterm_cent_z = 0.0;
  double gln_cent_x = 0.0;
  double gln_cent_y = 0.0;
  double gln_cent_z = 0.0;

  int cterm_posa = 88;
  int cterm_posb = 90;
  int gln_pos = 37; // this is 44 in 4G8E, 4G8F

  // try GLN 37 first
  bool gln_found = true;
  if (! GetCoordsRes(&my_tcra, gln_pos, " ", " CA ", "GLN", tmp_x, tmp_y, tmp_z))
    {
      gln_found = false;
      if (verbose) cout << "Unable to find residue " << gln_pos << " on TCRa chain" << endl;
    }
  else
    {
      gln_cent_x += tmp_x;
      gln_cent_y += tmp_y;
      gln_cent_z += tmp_z;
    }

  if (gln_found && (! GetCoordsRes(&my_tcrb, gln_pos, " ", " CA ", "GLN", tmp_x, tmp_y, tmp_z)))
    {
      gln_found = false;
      if (verbose) cout << "Unable to find residue " << gln_pos << " on TCRb chain" << endl;
    }
  else if (gln_found)
    {
      gln_cent_x += tmp_x;
      gln_cent_y += tmp_y;
      gln_cent_z += tmp_z;
    }
  
  // find the first c_term direction CA coord
  bool pos_found = true;
  if (! gln_found)
    {
      if (! GetCoordsRes(&my_tcra, cterm_posa, " ", " CA ", "TYR", tmp_x, tmp_y, tmp_z))
	{
	  if (GetCoordsRes(&my_tcra, cterm_posa + tcra_cys_shift, " ", " CA ", "TYR", tmp_x, tmp_y, tmp_z))
	    {
	      if (verbose) cout << "Using residue " << cterm_posa + tcra_cys_shift << " instead of " << cterm_posa << " for TCRa" << endl;
	    }
	  else
	    {
	      cerr << "Unable to find residue " << cterm_posa << " on TCRa" << endl;
	      pos_found = false;
	    }
	}
      cterm_cent_x += tmp_x;
      cterm_cent_y += tmp_y;
      cterm_cent_z += tmp_z;
      
      // find the c term direction CA coord
      if (! GetCoordsRes(&my_tcrb, cterm_posb, " ", " CA ", "TYR", tmp_x, tmp_y, tmp_z))
	{
	  if (GetCoordsRes(&my_tcrb, cterm_posb + tcrb_cys_shift, " ", " CA ", "TYR", tmp_x, tmp_y, tmp_z))
	    {
	      if (verbose) cout << "Using residue " << cterm_posb + tcrb_cys_shift << " instead of " << cterm_posb << " for TCRb" << endl;
	    }
	  else
	    {
	      cerr << "Unable to find residue " << cterm_posb << " on TCRb" << endl;
	      pos_found = false;
	    }
	}
      cterm_cent_x += tmp_x;
      cterm_cent_y += tmp_y;
      cterm_cent_z += tmp_z;
    }
  
  double pos_cent_x, pos_cent_y, pos_cent_z;
  if (gln_found)
    {
      pos_cent_x = gln_cent_x/2.0;
      pos_cent_y = gln_cent_y/2.0;
      pos_cent_z = gln_cent_z/2.0;
    }
  else if (pos_found)
    {
      pos_cent_x = cterm_cent_x/2.0;
      pos_cent_y = cterm_cent_y/2.0;
      pos_cent_z = cterm_cent_z/2.0;
    }
  else
    {
      cout << "Error: unable to find GLN or TYR positions for TCR vector orientation!" << endl;
      exit(1);
    }
  double tmp_vect_x = pos_cent_x - tcrcentx;
  double tmp_vect_y = pos_cent_y - tcrcenty;
  double tmp_vect_z = pos_cent_z - tcrcentz;
  
  // get the angle between the axis vector and the vector to those GLN residues
  double tmp_mag = sqrt(tmp_vect_x*tmp_vect_x + tmp_vect_y*tmp_vect_y + tmp_vect_z*tmp_vect_z);
  double axis_angle = acos((tmp_vect_x*tcrrotx + tmp_vect_y*tcrroty + tmp_vect_z*tcrrotz)/tmp_mag)*180.0/PI;

  double angle_cutoff = 25.0;
  if ((axis_angle > angle_cutoff) && (axis_angle < 180.0 - angle_cutoff))
    {
      cerr << "Error: tcr axis angle check too large: " << axis_angle << endl;
      exit(1);
    }
  if (verbose) cout << "Axis angle: " << axis_angle << endl;
  if (axis_angle < 90.0) 
    {
      tcrrotx = -tcrrotx;
      tcrroty = -tcrroty;
      tcrrotz = -tcrrotz;
    }
  return;
}

void TCRComplex::CalcMHCVectors(double & mhcx, double & mhcy, double & mhcz, double & mhccentx, double & mhccenty, double & mhccentz, double & mhcnormx, double & mhcnormy, double & mhcnormz)
{
  // calculate the vector in the direction of the MHC helices, the CM of those CA atoms, and the normal vector for the plane
  // from the Rudolph et al 2006 review, we need Ca atoms from A50-A86, A140-A176 for Class I
  // and A46-78, B54-64, B67-91 for Class II
  static const int MAX_ATMS = 100;
  static const int DIM = 3; // we are getting two weights, for x and y
  
  double tmp_x, tmp_y, tmp_z, chisq;
  double mhcorthx, mhcorthy, mhcorthz;
  double cent_x = 0.0, cent_y = 0.0, cent_z = 0.0;
  double mhc_x_coords[MAX_ATMS];
  double mhc_y_coords[MAX_ATMS];
  double mhc_z_coords[MAX_ATMS];

  // gsl variables
  gsl_matrix *A, *V;
  gsl_vector *S, *work;

  int num_coords = 0; // number of coords for all helices
  int num_helix_coords = 0; // number of coords for the first helix
  
  if (IsClassII)
    {
      for (int i = 46; i <= 78; i++)
	{
	  if (GetCoords(&my_mhc, i, " ", " CA ", tmp_x, tmp_y, tmp_z))
	    {
	      mhc_x_coords[num_coords] = tmp_x;
	      mhc_y_coords[num_coords] = tmp_y;
	      mhc_z_coords[num_coords] = tmp_z;
	      cent_x += tmp_x;
	      cent_y += tmp_y;
	      cent_z += tmp_z;
	      num_coords++;
	    }
	  else cout << "error: unable to find MHC atom: " << i << " for helix vector" << endl;
	}
      for (int i = 54; i <= 64; i++)
	{
	  if (GetCoords(&my_mhcb, i, " ", " CA ", tmp_x, tmp_y, tmp_z))
	    {
	      mhc_x_coords[num_coords] = tmp_x;
	      mhc_y_coords[num_coords] = tmp_y;
	      mhc_z_coords[num_coords] = tmp_z;
	      cent_x += tmp_x;
	      cent_y += tmp_y;
	      cent_z += tmp_z;
	      num_coords++;
	    }
	  else cout << "error: unable to find MHC atom: " << i << " for helix vector" << endl;
	}
      for (int i = 67; i <= 91; i++)
	{
	  if (GetCoords(&my_mhcb, i, " ", " CA ", tmp_x, tmp_y, tmp_z))
	    {
	      mhc_x_coords[num_coords] = tmp_x;
	      mhc_y_coords[num_coords] = tmp_y;
	      mhc_z_coords[num_coords] = tmp_z;
	      cent_x += tmp_x;
	      cent_y += tmp_y;
	      cent_z += tmp_z;
	      num_coords++;
	    }
	  else if (verbose) cout << "error: unable to find MHC atom: " << i << " for helix vector" << endl;
	}
    }
  else
    {
      int helix1_start = 50, helix1_end = 86, helix2_start = 140, helix2_end = 176;
      if (my_mhc_type == 2) // CD1d 3HUJ
	{
	  helix1_start = 52;
	  helix1_end = 89;
	  helix2_start = 141; 
	  helix2_end = 178;
	}
      else if (my_mhc_type == 3) // MR1
	{
	  helix1_start = 48;
	  helix1_end = 85;
	  helix2_start = 136;
	  helix2_end = 172;
	}
      else if (my_mhc_type == 4) // CD1d 4EI5
	{
	  helix1_start = 52;
	  helix1_end = 89;
	  helix2_start = 143;
	  helix2_end = 180; 
	}
      else if (my_mhc_type == 5) // CD1c
	{
	  helix1_start = 52;
	  helix1_end = 89;
	  helix2_start = 141;
	  helix2_end = 179;
	}
      for (int i = helix1_start; i <= helix1_end; i++)
	{
	  if (GetCoords(&my_mhc, i, " ", " CA ", tmp_x, tmp_y, tmp_z))
	    {
	      mhc_x_coords[num_coords] = tmp_x;
	      mhc_y_coords[num_coords] = tmp_y;
	      mhc_z_coords[num_coords] = tmp_z;
	      cent_x += tmp_x;
	      cent_y += tmp_y;
	      cent_z += tmp_z;
	      num_coords++;
	    }
	  else if (verbose) cout << "error: unable to find MHC atom: " << i << " for helix vector" << endl;
	}
      for (int i = helix2_start; i <= helix2_end; i++)
	{
	  if (GetCoords(&my_mhc, i, " ", " CA ", tmp_x, tmp_y, tmp_z))
	    {
	      mhc_x_coords[num_coords] = tmp_x;
	      mhc_y_coords[num_coords] = tmp_y;
	      mhc_z_coords[num_coords] = tmp_z;
	      cent_x += tmp_x;
	      cent_y += tmp_y;
	      cent_z += tmp_z;
	      num_coords++;
	    }
	  else if (verbose) cout << "error: unable to find MHC atom: " << i << " for helix vector" << endl;
	}
    }
  // calculate the centroids for all helices
  cent_x /= (double)num_coords;
  cent_y /= (double)num_coords;
  cent_z /= (double)num_coords;
 
  mhccentx = cent_x;
  mhccenty = cent_y;
  mhccentz = cent_z;

  // allocate the GSL matrices and stuff
  A = gsl_matrix_alloc(num_coords, DIM);
  V = gsl_matrix_alloc(DIM, DIM);
  S = gsl_vector_alloc(DIM);
  work = gsl_vector_alloc(DIM);
  
  // assign the atom coords to the GSL variables, centered at the centroid
  for (int i = 0; i < num_coords; i++)
    {
      gsl_matrix_set(A, i, 0, mhc_x_coords[i] - cent_x);
      gsl_matrix_set(A, i, 1, mhc_y_coords[i] - cent_y);
      gsl_matrix_set(A, i, 2, mhc_z_coords[i] - cent_z);
    }

  // perform the SV decomp
  gsl_linalg_SV_decomp(A, V, S, work);

  // cout << "singular values: " << gsl_vector_get(S, 0) << " " << gsl_vector_get(S, 1) << " " << gsl_vector_get(S, 2) << endl;
  
  mhcx = gsl_matrix_get(V, 0, 0);
  mhcy = gsl_matrix_get(V, 1, 0);
  mhcz = gsl_matrix_get(V, 2, 0);
  
  // normalize the vector
  double mag = sqrt(mhcx*mhcx + mhcy*mhcy + mhcz*mhcz);
  mhcx /= mag;
  mhcy /= mag;
  mhcz /= mag;

  // cout << "MHC vector: (" << mhcx << "," << mhcy << "," << mhcz << ")" << endl;

  // get the MHC normal vector
  mhcnormx = gsl_matrix_get(V, 0, 2);
  mhcnormy = gsl_matrix_get(V, 1, 2);
  mhcnormz = gsl_matrix_get(V, 2, 2);

   // normalize the vector
  mag = sqrt(mhcnormx*mhcnormx + mhcnormy*mhcnormy + mhcnormz*mhcnormz);
  mhcnormx /= mag;
  mhcnormy /= mag;
  mhcnormz /= mag;

  // determine whether the MHC vectors are pointing in the appropriate direction, decided to be 
  // from MHC center to residue 50
  int vect_res_num = 84;
  if (IsClassII) vect_res_num = 76;
  if (! GetCoords(&my_mhc, vect_res_num, " ", " CA ", tmp_x, tmp_y, tmp_z))
    {
      cerr << "Error: unable to obtain coordinates of MHC residue " << vect_res_num << " for vector calculations!" << endl;
      exit(1);
    }
  
  // compute the angle between the mhc helix vector and the vector to atom 84 or 76
  double vect_par_x = tmp_x - mhccentx;
  double vect_par_y = tmp_y - mhccenty;
  double vect_par_z = tmp_z - mhccentz;
  double mag_parallel = sqrt(vect_par_x*vect_par_x + vect_par_y*vect_par_y + vect_par_z*vect_par_z);
  
  double parallel_angle = acos((vect_par_x*mhcx + vect_par_y*mhcy + vect_par_z*mhcz)/mag_parallel)*180.0/PI;
  
   // compute the angle between the mhc norm vector and the vector to the base of the MHC helix domain
  if (IsClassII) 
    {
      if (! GetCoords(&my_mhcb, 13, " ", " CA ", tmp_x, tmp_y, tmp_z))
	{
	  cerr << "Error: unable to obtain coordinates of MHC residue for vector calculations!" << endl;
	  exit(1);
	}
    }
  else
    {
      if (! GetCoords(&my_mhc, 98, " ", " CA ", tmp_x, tmp_y, tmp_z))
	{
	  cerr << "Error: unable to obtain coordinates of MHC residue for vector calculations!" << endl;
	  exit(1);
	}
    }
  
  double vect_norm_x = tmp_x - mhccentx;
  double vect_norm_y = tmp_y - mhccenty;
  double vect_norm_z = tmp_z - mhccentz;
  double mag_norm = sqrt(vect_norm_x*vect_norm_x + vect_norm_y*vect_norm_y + vect_norm_z*vect_norm_z);
    
  double norm_angle = acos((vect_norm_x*mhcnormx + vect_norm_y*mhcnormy + vect_norm_z*mhcnormz)/mag_norm)*180.0/PI;
  
  // sanity check!
  if (verbose) cout << "angle = " << parallel_angle << " " << norm_angle << endl;
  double angle_cutoff = 28.0;
  if ((parallel_angle > angle_cutoff) && (parallel_angle < (180.0 - angle_cutoff)))
    {
      cerr << "error: MHC parallel vector angle out of range, too suspicious!" << endl;
      exit(1);
    }
  if (parallel_angle > 90.0) // switch the helix vector if necessary
    {
      mhcx = -mhcx;
      mhcy = -mhcy;
      mhcz = -mhcz;
    }
  // switching the norm angle may not be necessary!!!
  if ((norm_angle > angle_cutoff) && (norm_angle < (180.0 - angle_cutoff)))
    {
      if (verbose) cout << "warning: MHC norm vector angle out of range, this is suspicious!" << endl;
    }
  if (norm_angle < 90.0) // switch the norm vector if necessary
    {
      mhcnormx = -mhcnormx;
      mhcnormy = -mhcnormy;
      mhcnormz = -mhcnormz;
    }
   
  return;
}

void TCRComplex::CalcCrossProd(double x1, double y1, double z1, double x2, double y2, double z2, double & x3, double & y3, double & z3)
{
  x3 = y1*z2 - z1*y2;
  y3 = z1*x2 - x1*z2;
  z3 = x1*y2 - y1*x2;
  
  return;
}


void TCRComplex::OutputAxis(FILE* outfile, double x, double y, double z, double vectx, double vecty, double vectz, string chain)
{
  double step_length = 4.5; // 3 Angstroms per point along axis
  string tmpbuf = "ATOM      1  CA  GLY " + chain + "   1       9.419  10.702  55.946  1.00 48.05           C  ";
  string tmpbuf2 = "ATOM      2  C   GLY " + chain + "   1       9.419  10.702  55.946  1.00 48.05           C  ";

  int dummy_res_num = 1;
  
  for (int i = -10; i <= 10; i++)
    {
      if ((i == 0) && (chain == "H")) continue;
      double newx = x + vectx*step_length*(double)i;
      double newy = y + vecty*step_length*(double)i;
      double newz = z + vectz*step_length*(double)i;
      
      fprintf(outfile, "%s%4d%s%8.3f%8.3f%8.3f%s%s", tmpbuf.substr(0, 22).c_str(), dummy_res_num, tmpbuf.substr(26, 4).c_str(), newx, newy, newz,
	      tmpbuf.substr(54).c_str(), "\n");
      newx += vectx*step_length*0.3;
      newy += vecty*step_length*0.3;
      newz += vectz*step_length*0.3;
      fprintf(outfile, "%s%4d%s%8.3f%8.3f%8.3f%s%s", tmpbuf2.substr(0, 22).c_str(), dummy_res_num, tmpbuf2.substr(26, 4).c_str(), newx, newy, newz,
	      tmpbuf2.substr(54).c_str(), "\n");
      
      dummy_res_num++;
    }

  return;
}

// get the rotation axis and angle from a rotation matrix
// based on Martin Baker Java code for rotation matrix to axis/angle
double TCRComplex::GetRotationAxisAngle(float *rot, double & x_vec, double & y_vec, double & z_vec)
{
  double eps1 = 0.01; // rounding errors
  double eps2 = 0.1; // determine 0 or 180 degrees
  double angle = 0.0;

  double m00, m01, m02, m10, m11, m12, m20, m21, m22;
  m00 = rot[0];
  m01 = rot[1];
  m02 = rot[2];
  m10 = rot[3];
  m11 = rot[4];
  m12 = rot[5];
  m20 = rot[6];
  m21 = rot[7];
  m22 = rot[8];

  if ((fabs(m01 - m10) < eps1) && (fabs(m02 - m20) < eps1) && (fabs(m12 - m21) < eps1))
    {
      // singularity
      if ((fabs(m01 + m10) < eps2) && (fabs(m02 + m20) < eps2) && (fabs(m12 + m21) < eps2) && (fabs(m00 + m11 + m22 - 3.0) < eps2))
	{
	  // zero angle!
	  return angle;
	}
      angle = 180.0; // 180 degree rotation, need to find the axis now
      double xx = (m00 + 1.0)/2.0;
      double yy = (m11 + 1.0)/2.0;
      double zz = (m22 + 1.0)/2.0;
      double xy = (m01 + m10)/4.0;
      double xz = (m02 + m20)/4.0;
      double yz = (m12 + m21)/4.0;
      if ((xx > yy) && (xx > zz)) // m00 is the largest diagonal
	{
	  if (xx < eps1)
	    {
	      x_vec = 0.0;
	      y_vec = 0.7071;
	      z_vec = 0.7071;
	    }
	  else
	    {
	      x_vec = sqrt(xx);
	      y_vec = xy/x_vec;
	      z_vec = xz/x_vec;
	    }
	}
      else if (yy > zz) // m11 is the largest diagonal
	{
	  if (yy < eps1)
	    {
	      x_vec = 0.7071;
	      y_vec = 0.0;
	      z_vec = 0.7071;
	    }
	  else 
	    {
	      y_vec = sqrt(yy);
	      x_vec = xy/y_vec;
	      z_vec = yz/y_vec;
	    }
	}
      else // m22 is the largest diagonal
	{
	  if (zz < eps1)
	    {
	      x_vec = 0.7071;
	      y_vec = 0.7071;
	      z_vec = 0.0;
	    }
	  else
	    {
	      z_vec = sqrt(zz);
	      x_vec = xz/z_vec;
	      y_vec = yz/z_vec;
	    }
	}
      return angle;
    }
  double s = sqrt((m21 - m12)*(m21 - m12) + (m02 - m20)*(m02 - m20) + (m10 - m01)*(m10 - m01));
  if (s < 0.001) s = 1.0;
  angle = acos((m00 + m11 + m22 - 1.0)/2.0);
  x_vec = (m21 - m12)/s;
  y_vec = (m02 - m20)/s;
  z_vec = (m10 - m01)/s;
  return angle*180.0/PI;
}

// align the TCR rotation axis and disulf axes to the input axes and center
void TCRComplex::AlignTCRToAxis(FILE* outfile1, FILE* outfile2, double tcrrotx, double tcrroty, double tcrrotz, double tcrx, double tcry, double tcrz, double tcrcentx, double tcrcenty, double tcrcentz, double & newcentx, double & newcenty, double & newcentz, bool shift)
{
  // first, translate the TCR to the center
  int tot_atoms = my_tcra.num_atoms + my_tcrb.num_atoms;
  
  tot_atoms += my_mhc.num_atoms + my_pep.num_atoms; // also moving MHC+pep
  if (IsClassII)  tot_atoms += my_mhcb.num_atoms; // also moving MHC+pep
  tot_atoms++; // extra position to store centroid

  double* x_coords = new double[tot_atoms];
  double* y_coords = new double[tot_atoms];
  double* z_coords = new double[tot_atoms];

  int ind = 0;
  
  for (int i = 0; i < my_mhc.num_atoms; i++) 
    {
      x_coords[ind] = my_mhc.Atoms[i].x;
      y_coords[ind] = my_mhc.Atoms[i].y;
      z_coords[ind] = my_mhc.Atoms[i].z;
      ind++;
    }
  for (int i = 0; i < my_pep.num_atoms; i++) 
    {
      x_coords[ind] = my_pep.Atoms[i].x;
      y_coords[ind] = my_pep.Atoms[i].y;
      z_coords[ind] = my_pep.Atoms[i].z;
      ind++;
    }
  if (IsClassII)
    {
      for (int i = 0; i < my_mhcb.num_atoms; i++) 
	{
	  x_coords[ind] = my_mhcb.Atoms[i].x;
	  y_coords[ind] = my_mhcb.Atoms[i].y;
	  z_coords[ind] = my_mhcb.Atoms[i].z;
	  ind++;
	}
    }
  
  for (int i = 0; i < my_tcra.num_atoms; i++) 
    {
      x_coords[ind] = my_tcra.Atoms[i].x;
      y_coords[ind] = my_tcra.Atoms[i].y;
      z_coords[ind] = my_tcra.Atoms[i].z;
      ind++;
    }
  for (int i = 0; i < my_tcrb.num_atoms; i++) 
    {
      x_coords[ind] = my_tcrb.Atoms[i].x;
      y_coords[ind] = my_tcrb.Atoms[i].y;
      z_coords[ind] = my_tcrb.Atoms[i].z;
      ind++;
    }

  // add the centroid coordinates at the end
  x_coords[ind] = newcentx;
  y_coords[ind] = newcenty;
  z_coords[ind] = newcentz; 

  // translate to the center
  for (int i = 0; i < tot_atoms; i++)
    {
      x_coords[i] -= tcrcentx;
      y_coords[i] -= tcrcenty;
      z_coords[i] -= tcrcentz;
    }

  
  // rotate to align with the z axis
  AlignMolToVector(x_coords, y_coords, z_coords, tcrrotx, tcrroty, tcrrotz, tcrx, tcry, tcrz, 0, 0, 1, tot_atoms);

  // get the angle between the x of interchain vector and x-axis in the x-y plane
  double mag = sqrt(tcrx*tcrx + tcry*tcry);
  
  double rot_angle = -acos(tcrx/mag);
  if (tcry < 0.0) rot_angle = -rot_angle; // negate angle if vector is below the x-axis
  rot_angle -= PI/4.0; // rotate to average docking angle 45 degrees
  if (shift) rot_angle += -PI*my_tcr_cross_shift/180;

  // new rotation matrix around z_axis
  double a11 = cos(rot_angle);
  double a12 = -sin(rot_angle);
  double a13 = 0;
  double a21 = sin(rot_angle);
  double a22 = cos(rot_angle);
  double a23 = 0;
  double a31 = 0;
  double a32 = 0;
  double a33 = 1;

   // rotate around the axis
   for (int i = 0; i < tot_atoms; i++)
    {
      double newx = a11*x_coords[i] + a12*y_coords[i] + a13*z_coords[i];
      double newy = a21*x_coords[i] + a22*y_coords[i] + a23*z_coords[i];
      double newz = a31*x_coords[i] + a32*y_coords[i] + a33*z_coords[i];
      x_coords[i] = newx;
      y_coords[i] = newy;
      z_coords[i] = newz;
    }
   
   // tilt as necessary
   if (shift && (my_tcr_tilt_shift != 0))
   { 
       rot_angle = -(PI*my_tcr_tilt_shift/180);
       a11 = cos(rot_angle);
       a12 = 0;
       a13 = sin(rot_angle);
       a21 = 0;
       a22 = 1;
       a23 = 0;
       a31 = -sin(rot_angle);
       a32 = 0;
       a33 = cos(rot_angle);
       // rotate around y axis
       for (int i = 0; i < tot_atoms; i++)
	 {
	   double newx = a11*x_coords[i] + a12*y_coords[i] + a13*z_coords[i];
	   double newy = a21*x_coords[i] + a22*y_coords[i] + a23*z_coords[i];
	   double newz = a31*x_coords[i] + a32*y_coords[i] + a33*z_coords[i];
	   x_coords[i] = newx;
	   y_coords[i] = newy;
	   z_coords[i] = newz;
	 }
     }
   // output the new coordinates!
   ind = 0;
   
   for (int i = 0; i < my_mhc.num_atoms; i++) 
    {
      
      fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_mhc.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhc.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }
   for (int i = 0; i < my_pep.num_atoms; i++) 
    {
      
      fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_pep.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_pep.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }
   if (IsClassII)
     {
       for (int i = 0; i < my_mhcb.num_atoms; i++) 
	 {
	   
	   fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_mhcb.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhcb.Atoms[i].line.substr(54).c_str(), "\n");
	   ind++;
	 }
     }
   fprintf(outfile1, "TER\n");

   double x_shift_dist = 0.0;
   double z_shift_dist = 0.0;
   if (shift) 
     {
       x_shift_dist = my_tcrx_shift;
       z_shift_dist = my_tcrz_shift;
     } 
   for (int i = 0; i < my_tcra.num_atoms; i++) 
    {
      
      fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_tcra.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_tcra.Atoms[i].line.substr(54).c_str(), "\n");
      fprintf(outfile2, "%s%8.3f%8.3f%8.3f%s%s", my_tcra.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind]+x_shift_dist, y_coords[ind], z_coords[ind]+z_shift_dist, my_tcra.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }
   for (int i = 0; i < my_tcrb.num_atoms; i++) 
    {
      fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_tcrb.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_tcrb.Atoms[i].line.substr(54).c_str(), "\n");
      fprintf(outfile2, "%s%8.3f%8.3f%8.3f%s%s", my_tcrb.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind]+x_shift_dist, y_coords[ind], z_coords[ind]+z_shift_dist, my_tcrb.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }
   
   newcentx = x_coords[ind];
   newcenty = y_coords[ind];
   newcentz = z_coords[ind];
   
   delete[] x_coords;
   delete[] y_coords;
   delete[] z_coords;
}

// align the TCR rotation axis and disulf axes to the input axes and center
void TCRComplex::AlignMHCToAxis(FILE* outfile1, FILE* outfile2, FILE* outfile3, double mhcrotx, double mhcroty, double mhcrotz, double mhcx, double mhcy, double mhcz, double mhccentx, double mhccenty, double mhccentz, double & newcentx, double & newcenty, double & newcentz)
{
  double z_offset = my_z_offset; // distance to translate the MHC
  
  // first, translate the TCR to the center
  int tot_atoms = my_mhc.num_atoms + my_pep.num_atoms;
  if (IsClassII) tot_atoms += my_mhcb.num_atoms;
  tot_atoms += my_tcra.num_atoms + my_tcrb.num_atoms; // do the tcr chains too!
  tot_atoms++; // add the centroid coords at the end
  
  double *x_coords, *y_coords, *z_coords; 
  x_coords = new double[tot_atoms];
  y_coords = new double[tot_atoms];
  z_coords = new double[tot_atoms];

  int ind = 0;
  for (int i = 0; i < my_mhc.num_atoms; i++) 
    {
      x_coords[ind] = my_mhc.Atoms[i].x;
      y_coords[ind] = my_mhc.Atoms[i].y;
      z_coords[ind] = my_mhc.Atoms[i].z;
      ind++;
    }
  for (int i = 0; i < my_pep.num_atoms; i++) 
    {
      x_coords[ind] = my_pep.Atoms[i].x;
      y_coords[ind] = my_pep.Atoms[i].y;
      z_coords[ind] = my_pep.Atoms[i].z;
      ind++;
    }
  if (IsClassII)
    {
      for (int i = 0; i < my_mhcb.num_atoms; i++) 
	{
	  x_coords[ind] = my_mhcb.Atoms[i].x;
	  y_coords[ind] = my_mhcb.Atoms[i].y;
	  z_coords[ind] = my_mhcb.Atoms[i].z;
	  ind++;
	}
    }
    for (int i = 0; i < my_tcra.num_atoms; i++) 
    {
      x_coords[ind] = my_tcra.Atoms[i].x;
      y_coords[ind] = my_tcra.Atoms[i].y;
      z_coords[ind] = my_tcra.Atoms[i].z;
      ind++;
    }
  for (int i = 0; i < my_tcrb.num_atoms; i++) 
    {
      x_coords[ind] = my_tcrb.Atoms[i].x;
      y_coords[ind] = my_tcrb.Atoms[i].y;
      z_coords[ind] = my_tcrb.Atoms[i].z;
      ind++;
    }

  x_coords[ind] = newcentx;
  y_coords[ind] = newcenty;
  z_coords[ind] = newcentz;
  
  // translate to the center
  for (int i = 0; i < tot_atoms; i++)
    {
      x_coords[i] -= mhccentx;
      y_coords[i] -= mhccenty;
      z_coords[i] -= mhccentz;
    }
  
  // rotate to align with the z axis
  AlignMolToVector(x_coords, y_coords, z_coords, mhcrotx, mhcroty, mhcrotz, mhcx, mhcy, mhcz, 0, 0, -1, tot_atoms);
  
  // get the angle between the x of interchain vector and x-axis in the x-y plane
  // in this case the z component should be zero anyway as it is orthogonal to the mhc rotation vector
  double mag = sqrt(mhcx*mhcx + mhcy*mhcy);
  
  double rot_angle = -acos(mhcx/mag);
  if (mhcy < 0.0) rot_angle = -rot_angle; // negate angle if vector is below the x-axis
  
  // new rotation matrix around z-axis
  double a11 = cos(rot_angle);
  double a12 = -sin(rot_angle);
  double a13 = 0;
  double a21 = sin(rot_angle);
  double a22 = cos(rot_angle);
  double a23 = 0;
  double a31 = 0;
  double a32 = 0;
  double a33 = 1;

   // rotate around the axis
   for (int i = 0; i < tot_atoms; i++)
    {
      double newx = a11*x_coords[i] + a12*y_coords[i] + a13*z_coords[i];
      double newy = a21*x_coords[i] + a22*y_coords[i] + a23*z_coords[i];
      double newz = a31*x_coords[i] + a32*y_coords[i] + a33*z_coords[i];
      x_coords[i] = newx;
      y_coords[i] = newy;
      z_coords[i] = newz + z_offset;
    }
   
   // output the new coordinates!
   ind = 0;
 
   // output the pMHC
   for (int i = 0; i < my_mhc.num_atoms; i++) 
    {
      
      fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_mhc.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhc.Atoms[i].line.substr(54).c_str(), "\n");
      fprintf(outfile2, "%s%8.3f%8.3f%8.3f%s%s", my_mhc.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhc.Atoms[i].line.substr(54).c_str(), "\n");
      fprintf(outfile3, "%s%8.3f%8.3f%8.3f%s%s", my_mhc.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhc.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }
   for (int i = 0; i < my_pep.num_atoms; i++) 
    {
      
      fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_pep.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_pep.Atoms[i].line.substr(54).c_str(), "\n");
      fprintf(outfile2, "%s%8.3f%8.3f%8.3f%s%s", my_pep.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_pep.Atoms[i].line.substr(54).c_str(), "\n");
      fprintf(outfile3, "%s%8.3f%8.3f%8.3f%s%s", my_pep.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_pep.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }
   if (IsClassII)
     {
       for (int i = 0; i < my_mhcb.num_atoms; i++) 
	 {
	   
	   fprintf(outfile1, "%s%8.3f%8.3f%8.3f%s%s", my_mhcb.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhcb.Atoms[i].line.substr(54).c_str(), "\n");
	   fprintf(outfile2, "%s%8.3f%8.3f%8.3f%s%s", my_mhcb.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhcb.Atoms[i].line.substr(54).c_str(), "\n");
	   fprintf(outfile3, "%s%8.3f%8.3f%8.3f%s%s", my_mhcb.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_mhcb.Atoms[i].line.substr(54).c_str(), "\n");
	   ind++;
	 }
     }

   // print the TCR chains too!
   for (int i = 0; i < my_tcra.num_atoms; i++) 
    {
      
      fprintf(outfile2, "%s%8.3f%8.3f%8.3f%s%s", my_tcra.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_tcra.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }
   for (int i = 0; i < my_tcrb.num_atoms; i++) 
    {
      
      fprintf(outfile2, "%s%8.3f%8.3f%8.3f%s%s", my_tcrb.Atoms[i].line.substr(0, 30).c_str(),
	      x_coords[ind], y_coords[ind], z_coords[ind], my_tcrb.Atoms[i].line.substr(54).c_str(), "\n");
      ind++;
    }

   fprintf(outfile2, "TER\n");
  
   newcentx = x_coords[ind];
   newcenty = y_coords[ind];
   newcentz = z_coords[ind];
   
   delete[] x_coords;
   delete[] y_coords;
   delete[] z_coords;
}

// align set of coords to a given vector representing an axis
// also apply the same transform to an input vector, orig_axis
void TCRComplex::AlignMolToVector(double* x_coords, double* y_coords, double* z_coords, double orig_rotx, double orig_roty, double orig_rotz, double & orig_axisx, double & orig_axisy, double & orig_axisz, double axisx, double axisy, double axisz, int num_coords)
{
  double newrotx, newroty, newrotz;
  
  double mag1 = sqrt(axisx*axisx + axisy*axisy + axisz*axisz);
  double mag2 = sqrt(orig_rotx*orig_rotx + orig_roty*orig_roty + orig_rotz*orig_rotz);
    
  // determine rotation axis to align to the input axis
  CalcCrossProd(orig_rotx, orig_roty, orig_rotz, axisx, axisy, axisz, newrotx, newroty, newrotz);
  double rot_angle = acos((axisx*orig_rotx + axisy*orig_roty + axisz*orig_rotz)/(mag1*mag2));
  
  // need to normalize the cross product!
  double mag = sqrt(newrotx*newrotx + newroty*newroty + newrotz*newrotz);
  newrotx /= mag;
  newroty /= mag;
  newrotz /= mag;

  // let's try quaternions
  double q0 = cos(rot_angle/2.0);
  double q1 = sin(rot_angle/2.0)*newrotx;
  double q2 = sin(rot_angle/2.0)*newroty;
  double q3 = sin(rot_angle/2.0)*newrotz;
  
  double a11 = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  double a12 = 2.0*(q1*q2 - q0*q3);
  double a13 = 2.0*(q1*q3 + q0*q2);
  double a21 = 2.0*(q2*q1 + q0*q3);
  double a22 = q0*q0 - q1*q1 + q2*q2 - q3*q3;
  double a23 = 2.0*(q2*q3 - q0*q1);
  double a31 = 2.0*(q3*q1 - q0*q2);
  double a32 = 2.0*(q3*q2 + q0*q1);
  double a33 = q0*q0 - q1*q1 - q2*q2 + q3*q3;

  // rotate around the axis
   for (int i = 0; i < num_coords; i++)
    {
      double newx = a11*x_coords[i] + a12*y_coords[i] + a13*z_coords[i];
      double newy = a21*x_coords[i] + a22*y_coords[i] + a23*z_coords[i];
      double newz = a31*x_coords[i] + a32*y_coords[i] + a33*z_coords[i];
      x_coords[i] = newx;
      y_coords[i] = newy;
      z_coords[i] = newz;
    }

   double newx = a11*orig_axisx + a12*orig_axisy + a13*orig_axisz;
   double newy = a21*orig_axisx + a22*orig_axisy + a23*orig_axisz;
   double newz = a31*orig_axisx + a32*orig_axisy + a33*orig_axisz;
   orig_axisx = newx;
   orig_axisy = newy;
   orig_axisz = newz;
}
