/***********************************************************************
CreateProtein - Functions to create protein structures from amino acid
sequences and secondary structure prediction sequences.

Copyright (c) 2003 The Regents of the University of California, through
Lawrence Berkeley National Laboratory, University of California at
Davis, and Lawrence Livermore National Laboratory, subject to any
required approvals from the U.S. Department of Energy.

This source code is part of the ProteinShop software.

ProteinShop is copyrighted and your use is under license, subject to
any required approvals from the U.S. Department of Energy.  For details
or questions, you may contact Berkeley Lab's Technology Transfer
Department at TTD@lbl.gov (Re:  ProteinShop; CR-1877)

NOTICE OF U.S. GOVERNMENT RIGHTS.  ProteinShop was developed under
funding from the U.S. Government which consequently retains certain
rights as follows: the U.S. Government has been granted for itself and
others acting on its behalf a paid-up, nonexclusive, irrevocable,
worldwide license in ProteinShop to reproduce, prepare derivative
works, and perform publicly and display publicly.  Beginning five (5)
years after the date permission to assert copyright is obtained from the
U.S. Department of Energy, and subject to any subsequent five (5) year
renewals, the U.S. Government is granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in
ProteinShop to reproduce, prepare derivative works, distribute copies
to the public, perform publicly and display publicly, and to permit
others to do so.

Written by Nelson Max, Oliver Kreylos and James Lu.
***********************************************************************/

#ifndef CREATEPROTEIN_INCLUDED
#define CREATEPROTEIN_INCLUDED

namespace MD {

/* Forward declarations: */
class Protein;

void ReadStandards(const char* standardsDirectory); // Reads standard configuration of amino acid residues
int  checkStandardsAtomNo(const char* standardsDirectory, const char* aminoName);
int  matchStandardsAtoms(const char* standardsDirectory, const char* aminoName, const char*atomName, int index);
void setAlphaHelixAngles(double phi,double psi); // Sets default dihedral angles for alpha helices in degrees
void setBetaStrandAngles(double phi,double psi); // Sets default dihedral angles for beta strands in degrees
void setCoilRegionAngles(double phi,double psi); // Sets default dihedral angles for coil regions in degrees
Protein* ReadPredictionFile(const char* predictionFilename, int setprotein); // Reads prediction file and returns protein structure
Protein* SetDihedrals(int numResidues,const int type[],const char pred[],const double phis[],const double psis[]); // Creates protein structure from given amino acid sequence and dihedral angles
Protein* loadProtein(const char* inputFileName); // Creates protein from input file; file type is determined by extension
}

#endif
