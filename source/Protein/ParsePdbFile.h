/***********************************************************************
ParsePdbFile - Function to read a protein structure from a protein
database file.

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

Written by Oliver Kreylos and James Lu.
***********************************************************************/

#ifndef PARSEPDBFILE_INCLUDED
#define PARSEPDBFILE_INCLUDED

namespace MD {

/* Forward declarations: */
class Protein;

Protein* parsePdbFile (const char* filename,bool valid=true, const char* chain=0, const int modelId =0);
char* parseChains (const char *filename);
int   parseModels (const char *filename);
bool writePdbFile(const Protein& protein,const char* filename,bool writeStructure =true);
bool writePdbFile(const Protein& protein,int fd,bool writeStructure =true);
bool checkPdbFile(const char* filename);
bool checkHydrogen (const char *filename);
bool checkBackBones (const char *filename);
bool checkBackBone (const char *atomName);
}

#endif
