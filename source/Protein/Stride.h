/***********************************************************************
Stride - Class to invoke the Stride  

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

Written by James Lu.
***********************************************************************/

#ifndef _STRIDE_H
#define _STRIDE_H
 
class Stride {
 
public:
 
    Stride(void);
    ~Stride(void);
        
    int 	stridePrediction( const char *infilename );
    const char* 	parseStrideFile(const char *stridefile, const char* chainID); 
    int 	writePredFile(const char *stridefile); 
	void	enableChainSelection(void) { selectChain = true;}
	void 	clean(void);

private:
	char*	namebuffer;
	//char*	aaSeq;
	int 	count;
	bool   selectChain;
	char 	chainId[1];
  	char *outfilename;
};

#endif
