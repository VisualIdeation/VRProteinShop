/***********************************************************************
DragBox - Class for 6-DOF rigid body movement interaction.

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
***********************************************************************/

#ifndef DRAGBOX_INCLUDED
#define DRAGBOX_INCLUDED

#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/ProjectiveTransformation.h>

class DragBox
{
    /* Embedded classes: */
    public:
    typedef Geometry::Vector<double,3> Vector;
    typedef Geometry::Point<double,3> Point;
    typedef Geometry::OrthonormalTransformation<double,3> OTransformation;
    typedef Geometry::ProjectiveTransformation<double,3> PTransformation;
    
    enum PickMode // Enumerated type for picking modes
    {
        TRANSPARENT,OPAQUE
    };
    
    enum DragMode // Enumerated type for dragging modes
    {
        NONE,TRANSLATING,ROTATING_AXIS,ROTATING_VERTEX,SIXDOF
    };
    
    /* Elements: */
    private:
    Point center; // The drag box's center
    Vector axis[3]; // The drag box's three (normalized) axis vectors
    double size[3]; // The drag box's sizes along the axis vectors
    double edgeRadius; // The radius of an edge cylinder
    double vertexRadius; // The radius of a vertex sphere
    
    PickMode currentPickMode; // The current picking mode (determines how box is picked with 2D devices)
    DragMode currentDragMode; // The current dragging mode (gets set after a successful pick operation)
    bool edgeHighlightFlags[12]; // Highlight flags for each edge for rendering
    Point lastIntersection; // The last intersection point for dragging
    double rotateDepth; // Z distance between picked point and rotation center in clip coordinates
    Point lastMouse; // Last picked/dragged mouse position in clip coordinates
    Vector translateFaceNormal; // Normal vector of translating plane
    double translateFaceOffset; // Offset of translating plane
    Vector rotateAxis; // Axis of rotation
    Point rotateCenter; // Center of rotation
    double rotateCylinderRadius; // Radius of rotation cylinder
    OTransformation initialTransformation; // Transformation matrix to pre-apply when dragging with a 6-DOF dragger
    OTransformation dragTransformation; // Complete transformation matrix during dragging
    
    /* Private methods */
    double intersectVertex(int vertexIndex,const Point& start,const Vector& direction) const;
    double intersectEdge(int axisIndex,int vertexMask,const Point& start,const Vector& direction) const;
    double intersectFace(int axisIndex,int faceSign,const Point& start,const Vector& direction) const;
    
    /* Constructors and destructors: */
    public:
    DragBox(void); // Constructs a default drag box
    
    /* Methods: */
    const Point& getCenter(void) // Returns box's current center
    {
        return center;
    };
    void setCenter(const Point& newCenter); // Sets a new center point
    void setAxis(int axisIndex,const Vector& newAxis); // Sets one axis
    void setSize(int axisIndex,double newSize); // Sets one size
    void setEdgeRadius(double newEdgeRadius); // Sets a new radius for ray/edge intersection tests
    const Point& getRotateCenter(void) const // Returns current center of rotation
    {
        return rotateCenter;
    };
    void setRotateCenter(const Point& newRotateCenter); // Sets a new center of rotation
    PickMode getPickMode(void) const // Returns current picking mode
    {
        return currentPickMode;
    };
    void setPickMode(PickMode newPickMode) // Sets current picking mode
    {
        currentPickMode=newPickMode;
    };
    bool pickIncremental(void); // Prepares box for subsequent incremental 6-DOF dragging
    bool pick(const OTransformation& transformation); // Returns true if the transformation's origin is inside the box; prepares for subsequent 6-DOF box dragging
    bool pick(const Point& start,const Vector& direction); // Returns true if ray intersects box; sets up internal state for subsequent dragging
    bool pick(const PTransformation& modelView,const PTransformation& projection,const Point& mouseClip); // Returns true if mouse picks box; sets up internal state for dragging
    DragMode getDragMode(void) const // Returns current dragging mode
    {
        return currentDragMode;
    };
    void dragIncremental(const OTransformation& transformation); // Drags the box with an incremental 6-DOF dragger
    void drag(const OTransformation& transformation); // Drags the box with a 6-DOF dragger
    void drag(const Point& start,const Vector& direction); // Drags the box via a ray intersecting the 3D widget
    void drag(const PTransformation& modelView,const PTransformation& projection,const Point& mouseClip); // Drags the box via screen-space position
    const OTransformation& getDragTransformation(void) const // Returns the current dragging transformation
    {
        return dragTransformation;
    };
    void release(void); // Resets the dragging mode to NONE
    void draw(void) const; // Draws the box
};

#endif
