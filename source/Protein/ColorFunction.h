/***********************************************************************
ColorFunction - Map the unit interval to RGBA space.

Copyright (c) 2005 The Regents of the University of California, through
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

Written by Clark Crawford.
***********************************************************************/

#include <GL/GLColor.h>

class ColorFunction;

#ifndef COLOR_FUNCTION_INCLUDED
#define COLOR_FUNCTION_INCLUDED


typedef unsigned int uint;

/** Function object for arbitrary mappings from the unit interval to RGBA
    space.  Some default implementations are provided, all but one of which set
    full opacity (A = 1).  Singleton instances of these can be obtained from the
    class factory. */
class ColorFunction
{
public:

    /** Mapping range type, component vector in RGBA space. */
    typedef GLColor<double,4> Color;

    /** Remove this function from the factory's repertoire, if applicable. */
    virtual ~ColorFunction();

    /** The color function.
        @param scalar
            Clamped to the interval [0,1] before mapping.
        @param result
            The color corresponding to @a scalar is stored here. */
    virtual void map (double scalar, Color &result) = 0;

    /** Provide a name for this function suitable for utilization by a user
        interface component.
        @return
            User-friendly name in a null-terminated string. */
    virtual const char *name() = 0;

    /** Add a new function to the factory's repertoire.
        @param function
            Address of the object to add.  It will not be added if it is null or
            has already been added.
        @return
            The number that is assigned to this function, UINT_MAX if it is a
            duplicate. */
    static uint add (ColorFunction *function);

    /** Factory method for color functions.
        @param number
            Number of the function requested, in [0, numFunctions()).
        @return
            Pointer to the requested function, or null if @a number is out of
            range. */
    static ColorFunction *get (uint number);

    /** Get the number of functions available in the factory.
        @return
            The number of default functions, plus the number of functions that
            have been added, minus the number of functions that have been
            removed. */
    static uint numFunctions();

    /** Remove a function from the factory's repertoire.
        @param number
            Number of the function requested, in [0, numFunctions()).
        @return
            Pointer to the function that was removed, if @a number was in range;
            null otherwise.  The indices of all higher-numbered functions will
            decrease by 1 after removal. */
    static ColorFunction *remove (uint number);
};


#endif
