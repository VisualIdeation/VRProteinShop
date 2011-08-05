/***********************************************************************
Stride - Class to invoke the Stride 

Copyright (c) 2003 The Regents of the University of California, through
Lawrence Berkeley National Laboratory, University of California at
Davis, and Lawrence Livermore National Laboratory, subject to any
required approvals from the U.S. Department of Energy.
Adaptations to g++ 4.3.x copyright (c) 2009 Oliver Kreylos

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

#define MAX_RESIDUE 4000

#include <stdexcept>
#include <stdio.h>      
#include <string.h>
#include <stdlib.h>
#include "Stride.h"
#include <typeinfo>
#include "Globals.h"

Stride::Stride(void)
	:count(0), selectChain(false)
{
     namebuffer= new char[MAX_RESIDUE];
	 namebuffer[0] ='\0';
     chainId[0] ='\0';
  	outfilename = const_cast<char*>(".stridePredFile");
}
 
Stride::~Stride(void)
{
	delete namebuffer;
}

int Stride::stridePrediction( const char *infilename ) {

	#ifdef __APPLE__
  const char *stridebin   = "./bin/stride.mac";
	#endif
	#ifdef __LINUX__
  const char *stridebin   = "./bin/stride.linux";
	#endif
	#ifdef __SGI_IRIX__
  const char *stridebin   = "./bin/stride.irix";
	#endif

  /* Check the availability of Stride binary */
  if (!stridebin) {
	msg.Error(STRIDEID,"No STRIDE binary found!", stridebin, "in ./bin" );
    return 0;
  }

  char *stridecall = new char[strlen(stridebin)+strlen(infilename)
                               +strlen(outfilename)+8];
  
#ifdef CHSIN
  /* Set the specified chain for selection */
  if(chainID!=NULL) {
  	strncpy(chainId, chainID,1);
  	chainId[1] ='\0';
	msg.Debug(STRIDEID, chainId, " will be selected from the result of Stride.");
  }
  else
	msg.Debug(STRIDEID, "No chain selection for stride ", chainId);

#endif
  /* Generate command calls */
  sprintf(stridecall,"\"%s\" %s -f%s", stridebin, infilename, outfilename);

  /* Validate the Stide prediction */
  if (system(stridecall) <0 ){
  	msg.Fatal(STRIDEID,"Stride call failed");
    return 0;
  }
  
  /* Free memory */
  delete stridecall;
  
  msg.Debug(STRIDEID,"Stride Second Structure prediction\n", namebuffer);

  return 1;
} 

 
//int Stride::parseStrideFile(const char *filename)
const char* Stride::parseStrideFile(const char *filename, const char* chainID )
{ 
	filename = ".stridePredFile";
  /* Set the specified chain for selection */
  if(chainID!=NULL) {
  	strncpy(chainId, chainID,1);
  	chainId[1] ='\0';
	msg.Debug(STRIDEID, chainId, " will be selected from the result of Stride.");
  }
  else
	msg.Debug(STRIDEID, "No chain selection for stride ", chainId);

	/* open the the input stride file and check for validity: */
	FILE *file = fopen(filename, "rt");
	if (!file) {
        msg.Error(STRIDEID,"Unable to open the prediction file	");
    	return NULL;
	}
  	
	char line[100], keyword[4], residuename[4], chain[1], origResIDS[5], ss[1], aa[50];
  	int origResID, resID[1], currentID = 0;
	if(selectChain) 
  			msg.Debug(STRIDEID, chainId, " was specified for selection!");

	while (fgets(line, sizeof(line), file))
	{
         /* Parse the line just read: */
		sscanf(&line[ 0],"%3s",keyword);

		if (strcmp(keyword,"ASG")==0)
		{
			sscanf(line,"ASG  %3s %c %4s %4d    %c",
      			residuename, chain, origResIDS, resID, ss);
			#ifdef ORIGRESID
			origResID = atoi(origResIDS);  

  			if (origResID < 0)
			{
        		msg.Error(STRIDEID, "Invalid resid found in output file!",filename);
        		msg.Error(STRIDEID, ">> ",line);
     			
				fclose(file);
      			return 0;
    		}
			#endif
			
			/* Generate a sequence of Second Structure prediction from Stride */
			ss[1]='\0';
			if(!selectChain) 
			{
				strcat(namebuffer, ss);
				count++;
			}
			else
			{
				if(chainId != NULL) {
				/* Extract only the specified chain */
				if (strncmp(chain, chainId, 1)==0) 
				{
					strcat(namebuffer, ss);
					count++;
				}
				}
			}
		}
	}
    	
  	fclose(file);
	return namebuffer;
}  

void Stride::clean()
{
//  char *outfilename = ".stridePredFile";
  char *deletecall = new char[strlen("rm")+strlen(outfilename)+8];
  sprintf(deletecall,"\"%s\" %s", "rm", outfilename);
  /* Remove the output file from Stride */
  if (getenv("STRIDE")==NULL)
  system(deletecall);

  delete deletecall;

}

int Stride::writePredFile(const char *filename)
{ 
	/* open the the input stride file and check for validity: */
	FILE *file = fopen(filename, "w");
	if (!file) {
       // throw msg->Error(STRIDEID,"Stride: Unable to open input file", filename);
    	return 1;
	}
	fprintf(file,"Conf: ");
	for(int i = 0; i< count+1; i++) {
		if (i < count) {
		fprintf(file,"%c", '9');
		}
		else {
		fprintf(file,"\n");
		}
	}

  	fclose(file);
  	return 0;
}  

