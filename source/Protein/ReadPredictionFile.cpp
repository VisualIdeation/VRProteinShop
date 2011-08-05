/***********************************************************************
ReadPredictionFile - Functions to read prediction files and create their
protein structures.

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

Written by Nelson Max.
***********************************************************************/

#include <stdio.h>
#include <stdexcept>

#include <Misc/StandardValueCoders.h>

#include "Protein.h"
#include "CreateProtein.h"
#include "Globals.h"

#define PI 3.14159265358979323

static int strand = 0, nstrand = 0;
static int coil = 0, ncoil = 0;
static int strandstart[100], strandend[100];
static int coilstart[100], coilend[100];

#define MAXRES 1400
int residue_type[MAXRES] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
char pred[MAXRES] = {'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'} ;
static int numResidues = 10;
static double phi[MAXRES];
static double psi[MAXRES];
//Nelson' change
int conf[MAXRES] = { 
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5} ;
int is_core[MAXRES], has_core = 0;
//static int conf[MAXRES];
static double Cphi = -182.*PI/180;
static double Cpsi = -182.*PI/180;
static double raphi = -120.*PI/180.;
static double rapsi = 140.*PI/180.;
static double flphi = -131.85*PI/180.;
static double flpsi = 128.15*PI/180.;
static double Hphi = -60.*PI/180.;
static double Hpsi = -45.*PI/180.;
static double ProlinePhi = -60.*PI/180.;
extern char inputfilename[];

namespace MD {

Protein* ReadPredictionFile(const char* predictionFilename, int setprotein)
{
    int c,d,i,j, flatten;
    double Ephi, Epsi, flatfrac = 0;
	

    FILE* fp = fopen(predictionFilename, "r");
    strcpy (inputfilename, predictionFilename);
    if (!fp)
    {
		throw std::runtime_error(std::string("Unable to open ")+std::string(predictionFilename));
        return 0;
    }
	else
		msg.Info(PREDID, "Loading ", predictionFilename);
    
	numResidues = 0;
    flatten = configFile->retrieveValue<int,Misc::ValueCoder<int> >("/BuildBeta/flatten",0 );
    flatfrac = configFile->retrieveValue<float,Misc::ValueCoder<float> >("/BuildBeta/flatfrac",0.);
    if (!flatten) flatfrac = 0.;
    Ephi = flatfrac*flphi + (1. - flatfrac)*raphi;
    Epsi = flatfrac*flpsi + (1. - flatfrac)*rapsi;
    // printf ("flatfrac %f Ephi %f Epsi %f\n", flatfrac, Ephi, Epsi);
    while((c = fgetc(fp)) != ':')
        ;
    while((c = fgetc(fp)) == ' ')
        ;
    while(c >= '0' && c <= '9' && numResidues < MAXRES)
    {
        conf[numResidues] = c - '0';
        c = fgetc(fp);
        ++numResidues;
    }
    i = 0;
    while((c = fgetc(fp)) != ':')
        ;
    while((c = fgetc(fp)) == ' ')
        ;
    while((c == 'C' || c == 'E' || c == 'H') && i < MAXRES)
    {
        pred[i] = c;
        if(c == 'C')
        {
            phi[i] = Cphi;
            psi[i] = Cpsi;
            strand = 0;
            if(!coil)
            {
                coil = 1;
                ++ncoil;
                coilstart[ncoil] = i;
            }
            else
                coilend[ncoil] = i;
        }
        else if(c == 'E')
        {
            phi[i] = Ephi;
            psi[i] = Epsi;
            coil = 0;
            if(!strand)
            {
                strand = 1;
                ++nstrand;
                strandstart[nstrand] = i;
            }
            else
                strandend[nstrand] = i;
        }
        else if(c == 'H')
        {
            phi[i] = Hphi;
            psi[i] = Hpsi;
            coil = 0;
            strand = 0;
        }
        c = fgetc(fp);
        ++i;
    }
    if (i != numResidues)
    {
        printf("numbers of confidences %d and predictions %d do not agree\n", numResidues, i);
        return 0;
    }
    i = 0;
    while((c = fgetc(fp)) != ':')
        ;
    while((c = fgetc(fp)) == ' ')
        ;
    while(c > 64 && i < MAXRES)
    {
        for( j = 0; j < 25; ++j)
        {
            d = MD::Protein::Residue::singleLetterNames[j];
            if (c == d || c + ('a' - 'A') == d)
            {
                residue_type[i] = j;
                //if(j == 18) phi[i] = ProlinePhi;
                break;
            }
        }
        if(j == 25)
        {
            printf("unknown single letter residue name %c\n", c);
            return 0;
        }
        c = fgetc(fp);
        ++i;
    }
    if (i != numResidues)
    {
        printf("numbers of confidences %d and single letter names %d do not agree\n", numResidues, i);
        return 0;
    }
    printf("setprotein in ReadPredictionFile is %d\n", setprotein);
    if(setprotein)
	return SetDihedrals( numResidues, residue_type, pred, phi, psi);
    else
        return 0;
}

}
