/***********************************************************************
Globals - Global variables for interactive protein manipulation program.

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

Written by Oliver Kreylos.
Modified by Clark Crawford and James Lu.
***********************************************************************/


class ProteinState;


#ifndef GLOBALS_INCLUDED
#define GLOBALS_INCLUDED

#ifdef __GNUC__
/***************************************************************
This must be some bug in the STL implementation coming with g++:
operator== for lists is declared as a friend function, but not
using the new(?) C++-syntax for template friend functions.
**************************************************************/

namespace std {
template <class _Tp, class _Alloc>
class list;
template <class _Tp, class _Alloc>
bool operator==(const list<_Tp,_Alloc>& __x,const list<_Tp,_Alloc>& __y);
}
#endif
#include <list>
#include <pthread.h>

#include <Geometry/Box.h>
#include <Misc/ConfigurationFile.h>

#include "DragBox.h"
#include "Protein.h"
#include "ProteinClient.h"
#include "ProteinInteractor.h"
#include "ProteinRenderer.h"
#include "EnergyRenderer.h"
#include "EnergyAPI.h"
#include "UndoBuffer.h"
#include "Message.h"

/* Forward declarations: */
class GLContextData;

// proteinMutex guards or regulates the variables in this block except ikUpdatePosted
extern pthread_mutex_t proteinMutex;
extern pthread_cond_t ikUpdateRequestedCond;
extern volatile bool ikUpdateRequested;
extern pthread_cond_t ikUpdateDoneCond;
extern volatile bool ikUpdateDone;
extern volatile bool ikUpdatePosted;
extern pthread_t ikUpdateThreadId;

extern Misc::ConfigurationFile * configFile; // Configuration file containing program parameters
extern UndoBuffer undoBuffer;
extern EnergyLibrary* energyLibrary; // Object representing a dynamically linked energy calculation library
extern PShopMessage msg;

extern void updateDihedralAngles(void);
extern void recalculateEnergy(void);
extern void updateAtomEnergies(void);
extern double getEnergyOutput(void);
extern void updateGui(void);
extern int zipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern int zipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern int altZipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern int altZipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern void redrawRenderWindows(void);
extern void updateProteinNow(void);


struct ProteinState
{
    MD::Protein * protein; // The visualized protein
    char *name; // name to use when displaying a list of proteins

    ProteinClient* client; // A protein client talking to a global optimization server
    uint proteinId; // Unique ID of current protein in server's optimization tree

    MD::Protein::Residue* selectedResidue; // The selected residue
    ProteinInteractor* interactor; // Interaction object

    MD::ProteinRenderer* proteinRenderer; // Renderer attached to the protein
    EnergyRenderer *energyRenderer; // energy renderer attached to the protein

    EnergyCalculator* energyCalculator; // Energy calculator attached to the protein
    bool visualizeEnergy; // Flag to toggle energy visualization
    double visualizeEnergyMinRange,visualizeEnergyMaxRange, engUpdateRate; // Value range of energy visualization color map

    ProteinState (const char *filename);
    ~ProteinState();
    MD::Protein * getProtein(void) {
    	return protein;
    }
};

extern void checkDrawToggles(ProteinState* state);
extern void checkShowToggles(ProteinState* state);

// ProteinState managment
extern ProteinState *createProtein (const char *filename);
extern ProteinState *curProtein();
extern void deleteAllProteins();
extern bool deleteProtein (ProteinState *state);
extern ProteinState *getProtein (uint index);
extern uint numProteins();
extern bool setCurProtein (ProteinState *state);

/// Translate the errno variable into a string.
extern const char *errnoString();

/// Get the next larger power of two.
extern uint ceilingPowerOfTwo (uint number);

/// Clamp a scalar variable to an interval.
template <typename Scalar>
void clamp (Scalar &toClamp, const Scalar &lowerBound, const Scalar &upperBound)
{
    if ( toClamp < lowerBound ) toClamp = lowerBound;
    if ( toClamp > upperBound ) toClamp = upperBound;
}


#endif
