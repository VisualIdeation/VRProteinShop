/***********************************************************************
SetDihedrals - Reads angles and residues from arrays and builds a PDB
protein structure.

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

Written by Nelson Max and James Lu.
***********************************************************************/

/* Change this to !=0 to write output file "AlphaBeta.pdb": */
#define WRITEPDBFILE 0
#define WRITEHELIXFILE 0

#include <ctype.h>
#include <stdio.h>
#include <iostream>
#include <string.h>

#include "Protein.h"
#include "CreateProtein.h"
#define PI 3.14159265358979323

extern double residueAtomPos[25][25][3];
extern char residueAtomName[25][25][5];
extern double residued2[25];
extern double residued3[25];
extern double StandardPhi[25];
extern double StandardPsi[25];
extern double StandardAlpha[25];
extern double StandardBeta[25];
extern double StandardGamma[25];
extern int numbAtoms[25];

extern void normalize( double a[3]);
extern void cross( double a[3], double b[3], double c[3]);
extern double dot(double a[3], double b[3]);
double saxis[3], stail[3], ptail[3], dtail, daxis, dptail;

void identity4 (double a[4][4]) {
   int i, j;
   for (i = 0; i < 4; ++i) {
     for (j = 0; j < 4; ++j) a[i][j] = 0.;
     a[i][i] = 1.;
     }
   }

void translate4(double a[3], double b[4][4]) {
   int i;
   identity4(b);
   for (i = 0; i < 3; ++i) b[i][3] = a[i];
   }

void translate4d(double x, double y, double z, double b[4][4]) {
   int i;
   identity4(b);
   b[0][3] = x;
   b[1][3] = y;
   b[2][3] = z;
   }

void matmult4 (double a[4][4], double b[4][4], double c[4][4]) {

   int i, j, k;

   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 4; ++k)
            c[i][j] += a[i][k] * b[k][j];
         }
   }

void matmult_transp4 (double a[4][4], double b[4][4], double c[4][4]) {

   int i, j, k;

   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 4; ++k)
            c[i][j] += a[i][k] * b[j][k];
         }
   }

void matrix_vector4 (double a[4][4], double b[4], double c[4]) {

   int i, k;

   for (i = 0; i < 4; ++i) {
      c[i] = 0.;
      for (k = 0; k < 4; ++k)
         c[i] += a[i][k] * b[k];
      }
   }

void rotX4(double angle, double a[4][4]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
         a[i][j] = 0;
   a[1][1] = a[2][2] = cosa;
   a[1][2] = -sina;
   a[2][1] = sina;
   a[0][0] = 1.;
   a[3][3] = 1.;
   }

void rotY4(double angle, double a[4][4]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
         a[i][j] = 0;
   a[0][0] = a[2][2] = cosa;
   a[2][0] = -sina;
   a[0][2] = sina;
   a[1][1] = 1.;
   a[3][3] = 1.;
   }

void rotZ4(double angle, double a[4][4]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
         a[i][j] = 0;
   a[1][1] = a[0][0] = cosa;
   a[0][1] = -sina;
   a[1][0] = sina;
   a[2][2] = 1.;
   a[3][3] = 1.;
   }

namespace MD {

void SetStandardHelix() {
   int numResidues = 9, type[9] = {8, 8, 8, 8, 8, 8, 8, 8, 8};
   double phi[9], psi[9], head[3], tail[3], diff[3];
   int i, j, k, l, m, n;
   double v0[4], v1[4], v2[4], v3[4], v4[4], pos[4];
   double mainpos[4]={0.0,0.0,0.0,0.0};
   double a[4][4], b[4][4], c[4][4], d[4][4], e[4][4], f[4][4], g[4][4];
   double chainN[9][3];
   double chainCA[9][3];
   double chainC[9][3];
   double center[9][3];
   double axis[3], plane[3], dd, radius;
   #if WRITEHELIXFILE
   FILE *fp;
   const char outfile[34] = "StandardHelix.pdb";
   #endif
   char chainId[2], pc;
   char elementName[2];
   int residueIndex;
   char* atomNamePtr;


   for (i = 0; i < numResidues; ++i) {
      phi[i] = -60.*PI/180.;
      psi[i] = -45.*PI/180.;
      }
   #if WRITEHELIXFILE
   fp = fopen(outfile, "w");
   if (fp == 0) {
      printf("unable to open file %s\n", outfile);
      }
   #endif
   n = 1;
   identity4(a);

//  Standard amino acids have CA at the origin and N at -residue_d2,
//  so translate to put N at the origin. Standard amino acids also
//  have C in the xy plane, but not necessarily carbonyl O.

   for (i = 0; i < numResidues; ++i)
         {
         j = type[i];

         translate4d(residued2[j], 0., 0., b);
         matmult4(a, b, c);
         rotX4(PI-StandardPhi[j], b);
         //   printf("StandardPhi[%d] = %f\n", j, StandardPhi[j]);
         matmult4(c, b, f);              // f is for the amide H

         //  First do NH

         l = 2;
         for (k = 0; k < l; ++k)
         {
             for (m = 0; m < 3; ++m)
                 v0[m] = residueAtomPos[j][k][m];
             v0[3] = 1.;
             if(k == 0) {
                matrix_vector4(c, v0, mainpos);
                for (m = 0; m < 3; ++m) 
                   chainN[i][m] = mainpos[m];
                }
             else  matrix_vector4(f, v0, mainpos);
             #if WRITEHELIXFILE
             fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                     n, residueAtomName[j][k], Protein::Residue::abbreviatedNames[j],
                     i+1, mainpos[0], mainpos[1], mainpos[2], 0., 0.);
             #endif
             ++n;

         }
         rotX4(phi[i], b);
         matmult4(c, b, a);              // a continuing backbone
         rotZ4(StandardAlpha[j], b);
         matmult4(a, b, c);
         rotX4(psi[i], b);
         matmult4(c, b, d);              // d is to help continue backbone
         rotX4(StandardPsi[j], b);      // to bring O to plane of next N and CA
         matmult4(d, b, g);
         rotZ4(-StandardAlpha[j], b);
         matmult4(g, b, e);              // e is for carbonyl oxygen

         //  Now do rest of atoms

         for (k = l; k <  numbAtoms[j]; ++k)
             {
             for (m = 0; m < 3; ++m)
                 v0[m] = residueAtomPos[j][k][m];
             v0[3] = 1.;
             if (k == l+3)
                 matrix_vector4(e, v0, pos);
             else
                 matrix_vector4(a, v0, pos);
	     if(k == 2)
	        for (m = 0; m < 3; ++m)
		   chainCA[i][m] = pos[m];
	     if(k == 4)
		for (m = 0; m < 3; ++m)
		   chainC[i][m] = pos[m];
             #if WRITEHELIXFILE
             fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                     n, residueAtomName[j][k], Protein::Residue::abbreviatedNames[j],
                     i+1, pos[0], pos[1], pos[2], 0., 0.);
             #endif
             ++n;
             }
         // Now build up rest of main chain effect on matrix a

         translate4d(residued3[j], 0., 0., b);
         matmult4(d, b, c);
         rotZ4(StandardBeta[j], b);
         matmult4(c, b, d);
         translate4d(1.325, 0., 0., b);
         matmult4(d, b, c);
         rotX4(PI, b);
         matmult4(c, b, d);
         rotZ4(StandardGamma[j], b);
         matmult4(d, b, a);
         //   printf("Alpha %f  Beta %f  Gamma %f\n", 180*StandardAlpha[j]/PI,
         //      180*StandardBeta[j]/PI, 180*StandardGamma[j]/PI);
         }  
   for (i = 1; i < numResidues-1; ++i) {
      for(m = 0; m < 3; ++m) {
         v0[m] = chainN[i-1][m] - chainN[i][m];
	 v1[m] = chainN[i+1][m] - chainN[i][m];
	 }
      normalize(v0);
      normalize(v1);
      for(m = 0; m < 3; ++m)
         center[i][m] = v0[m] + v1[m];
      normalize(center[i]);
      }
   cross(center[1], center[2], axis);
   normalize(axis);
   cross(axis, center[2], plane);
   dd = dot(plane, chainN[2]) - dot(plane, chainN[1]) ;
   radius = dd/dot(plane,center[1]);
   for(m = 0; m < 3; ++m) {
      tail[m] = chainN[1][m] + radius*center[1][m];
      head[m] = chainN[7][m] + radius*center[7][m];
      diff[m] = (head[m] - tail[m])/6.;
      tail[m] -= diff[m];
      }
   for(m = 0; m < 3; ++m) {
     v0[m] = chainN[0][m] - chainCA[0][m];
     v1[m] = chainC[0][m] - chainCA[0][m];
     }
   normalize(v0);
   normalize(v1);
   dd = dot(v0, v1);
   for(m = 0; m < 3; ++m)
      v1[m] -= dd*v0[m];
   normalize(v1);
   cross(v0, v1, v2);
   normalize(v2);
   saxis[0] = dot(v0, diff);
   saxis[1] = dot(v1, diff);
   saxis[2] = dot(v2, diff);
   stail[0] = dot(v0, tail);
   stail[1] = dot(v1, tail);
   stail[2] = dot(v2, tail);
   cross(saxis, stail, ptail);
   normalize(ptail);
   dtail= sqrt( stail[0]*stail[0] + stail[1]*stail[1] + stail[2]*stail[2] );
   daxis= sqrt( saxis[0]*saxis[0] + saxis[1]*saxis[1] + saxis[2]*saxis[2] );
   for(m = 0; m < 3; ++m)
      ptail[m] = dtail*ptail[m];

/*
   normalize(diff);
   fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           61, "  C", Protein::Residue::abbreviatedNames[j],
           9, head[0], head[1], head[2], 0., 0.);
   fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           62, "  O", Protein::Residue::abbreviatedNames[j],
           9, tail[0], tail[1], tail[2], 0., 0.);
   cout << "radius  " << radius << "  dd " << dd << "  plane " << *plane <<
      "  center[0] " << *center[1] << "  tail " << *tail << "\n";

   fprintf(stderr, "axis %g %g %g\n", axis[0], axis[1], axis[2]);
   fprintf(stderr, "diff %g %g %g\n", diff[0], diff[1], diff[2]);
   fprintf(stderr, "head %g %g %g\n", head[0], head[1], head[2]);
   fprintf(stderr, "tail %g %g %g\n", tail[0], tail[1], tail[2]);
   fprintf(stderr, "chainN[7] %g %g %g\n", chainN[7][0], chainN[7][1], chainN[7][2]);
   fprintf(stderr, "chainN[1] %g %g %g\n", chainN[1][0], chainN[1][1], chainN[1][2]);
*/
   #if WRITEHELIXFILE
   fclose(fp);
   #endif
   }


Protein* SetDihedrals (int numResidues, const int type[],
      const char pred[], const double phi[], const double psi[]) {
   int i, j, k, l, m, n;
   double v0[4], v1[4], v2[4], v3[4], v4[4], pos[4];
     double mainpos[4]={0.0,0.0,0.0,0.0};
   double a[4][4], b[4][4], c[4][4], d[4][4], e[4][4], f[4][4], g[4][4];
   #if WRITEPDBFILE
   FILE *fp;
   const char outfile[14] = "StandardHelix.pdb";
   #endif
   Protein* result=new Protein;
   int currentResidueIndex=-1;
   Protein::ResidueCreator proteinCreator(result);
     Protein::SecondaryStructure::StructureType currentStructureType=Protein::SecondaryStructure::NONE;
   char chainId[2], pc;
   char elementName[2];
   int residueIndex;
   char* atomNamePtr;


   if(numResidues <= 0) return result;
   #if WRITEPDBFILE
   fp = fopen(outfile, "w");
   if (fp == 0) {
      printf("unable to open file %s\n", outfile);
      assert (fp != 0);
      }
   #endif
  n = 1;
  identity4(a);

//  Standard amino acids have CA at the origin and N at -residue_d2,
//  so translate to put N at the origin. Standard amino acids also
//  have C in the xy plane, but not necessarily carbonyl O.

    for (i = -1; i < numResidues+1; ++i)
    {
        if(i==-1||i==numResidues)
        {
            const char* residuePdbName;
            if (i == -1)
            {
                j = 0;
                residuePdbName="ACE";
            }
            else
            {
                j = 14;
                residuePdbName="NME";
            }
            Protein::SecondaryStructure::StructureType newStructureType;
            newStructureType=Protein::SecondaryStructure::COIL;
            if(newStructureType!=currentStructureType)
            {
                proteinCreator.newSecondaryStructure(newStructureType);
                currentStructureType=newStructureType;
            }
            proteinCreator.newResidue(residuePdbName,i+1);
            for (k = 0; k <  numbAtoms[j]; ++k)
            {
                for (m = 0; m < 3; ++m)
                    pos[m] = residueAtomPos[j][k][m] + mainpos[m];
                pos[1] += 1;
                pos[0] -= 0.8;
                #if WRITEPDBFILE
                fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                        n, residueAtomName[j][k], residuePdbName,
                        i+1, pos[0], pos[1], pos[2], 0., 0.);
                #endif
                ++n;
                atomNamePtr=residueAtomName[j][k];
                while(isspace(*atomNamePtr))
                    ++atomNamePtr;
                elementName[0]=*atomNamePtr;
                ++atomNamePtr;
                elementName[1]='\0';
                proteinCreator.addAtom(elementName, n-1 ,Position(pos), atomNamePtr);
            }
        }
        else
        {
            j = type[i];
            pc = pred[i];
            Protein::SecondaryStructure::StructureType newStructureType;
            if (pc == 'C')
                newStructureType=Protein::SecondaryStructure::COIL;
            else if (pc == 'H')
                newStructureType=Protein::SecondaryStructure::ALPHA_HELIX;
            else if (pc == 'E')
                newStructureType=Protein::SecondaryStructure::BETA_STRAND;
            else
            {
                printf("Unknown secondary structure type.\n");
                exit(-1);
            }
            if(newStructureType!=currentStructureType)
            {
                proteinCreator.newSecondaryStructure(newStructureType);
                currentStructureType=newStructureType;
            }
            proteinCreator.newResidue(Protein::Residue::abbreviatedNames[j],i+1);

            translate4d(residued2[j], 0., 0., b);
            matmult4(a, b, c);
            rotX4(PI-StandardPhi[j], b);
            //   printf("StandardPhi[%d] = %f\n", j, StandardPhi[j]);
            matmult4(c, b, f);              // f is for the amide H

            //  First do NH

            if(j == 16)
                l = 1; 
            else
                l = 2;
            for (k = 0; k < l; ++k)
            {
                for (m = 0; m < 3; ++m)
                    v0[m] = residueAtomPos[j][k][m];
                v0[3] = 1.;
                if(k == 0) matrix_vector4(c, v0, mainpos);
                else  matrix_vector4(f, v0, mainpos);
                #if WRITEPDBFILE
                fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                        n, residueAtomName[j][k], Protein::Residue::abbreviatedNames[j],
                        i+1, mainpos[0], mainpos[1], mainpos[2], 0., 0.);
                #endif
                ++n;

                atomNamePtr=residueAtomName[j][k];
                while(isspace(*atomNamePtr))
                    ++atomNamePtr;
                elementName[0]=*atomNamePtr;
                ++atomNamePtr;
                elementName[1]='\0';
                proteinCreator.addAtom(elementName, n-1 ,Position(mainpos), atomNamePtr);
                //  printf("%d %4s %4s %4s\n",
                //            n-1, elementName, atomNamePtr, residueAtomName[j][k]);
            }
            rotX4(phi[i], b);
            matmult4(c, b, a);              // a continuing backbone
            rotZ4(StandardAlpha[j], b);
            matmult4(a, b, c);
            rotX4(psi[i], b);
            matmult4(c, b, d);              // d is to help continue backbone
            rotX4(StandardPsi[j], b);      // to bring O to plane of next N and CA
            matmult4(d, b, g);
            rotZ4(-StandardAlpha[j], b);
            matmult4(g, b, e);              // e is for carbonyl oxygen

            //  Now do rest of atoms

            for (k = l; k <  numbAtoms[j]; ++k)
            {
                for (m = 0; m < 3; ++m)
                    v0[m] = residueAtomPos[j][k][m];
                v0[3] = 1.;
                if (k == l+3)
                    matrix_vector4(e, v0, pos);
                else
                    matrix_vector4(a, v0, pos);
                #if WRITEPDBFILE
                fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                        n, residueAtomName[j][k], Protein::Residue::abbreviatedNames[j],
                        i+1, pos[0], pos[1], pos[2], 0., 0.);
                #endif
                ++n;
                atomNamePtr=residueAtomName[j][k];
                while(isspace(*atomNamePtr))
                    ++atomNamePtr;
                elementName[0]=*atomNamePtr;
                ++atomNamePtr;
                elementName[1]='\0';
                proteinCreator.addAtom(elementName, n-1 ,Position(pos), atomNamePtr);
            }
            // Now build up rest of main chain effect on matrix a

            translate4d(residued3[j], 0., 0., b);
            matmult4(d, b, c);
            rotZ4(StandardBeta[j], b);
            matmult4(c, b, d);
            translate4d(1.325, 0., 0., b);
            matmult4(d, b, c);
            rotX4(PI, b);
            matmult4(c, b, d);
            rotZ4(StandardGamma[j], b);
            matmult4(d, b, a);
            //   printf("Alpha %f  Beta %f  Gamma %f\n", 180*StandardAlpha[j]/PI,
            //      180*StandardBeta[j]/PI, 180*StandardGamma[j]/PI);
        }
    }

  #if WRITEPDBFILE
  fclose(fp);
  #endif
  proteinCreator.finishProtein();
  SetStandardHelix();
  return result;
  }

}
