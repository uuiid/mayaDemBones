#include "doodleWeight.h"
#include <maya/MSyntax.h>
DoodleWeight::DoodleWeight( )
{
}

DoodleWeight::~DoodleWeight( )
{
}

void* DoodleWeight::creator( )
{
    return new DoodleWeight( );
}

MStatus DoodleWeight::doIt(const MArgList& arg)
{
    MStatus status = MS::kSuccess;




    return status;
}

bool DoodleWeight::isUndoable( ) const
{
    return false;
}

bool DoodleWeight::hasSyntax( ) const
{
    return true;
}

MSyntax DoodleWeight::newSyntax( )
{
    MSyntax syntax;
    return syntax;
}
