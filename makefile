# Brian Pierce
# makefile for protein complex monte carlo minimization program

CXX		=	g++
GSL_LIBS	= 	-lgsl -lgslcblas
FAST_LIBS	= 	-Lfast -lfast
INCLUDES	=	 

tcr_docking_angle: tcr_docking_angle.cc tcr_complex.cc tcr_complex.h makefile
	$(CXX) $(INCLUDES) -o $@ tcr_docking_angle.cc tcr_complex.cc $(FAST_INCLUDES) $(GSL_LIBS) $(FAST_LIBS)

