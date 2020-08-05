#include "doodleWeight.h"
#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MSyntax.h>
#include <maya/MSelectionList.h>
#include <maya/MArgDatabase.h>
#include <maya/MDagPath.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnMesh.h>
DoodleWeight::DoodleWeight( )
{
}

DoodleWeight::~DoodleWeight( )
{
}

void* DoodleWeight::creator( )
{
    return new DoodleWeight;
}

MStatus DoodleWeight::doIt(const MArgList& arg)
{
    MStatus status = MS::kSuccess;
    
    CHECK_MSTATUS_AND_RETURN_IT(this->getsyntaxFlag(arg));

    MSelectionList selectList;
    CHECK_MSTATUS_AND_RETURN_IT(selectList.add(commNodeString));
    CHECK_MSTATUS_AND_RETURN_IT(selectList.add(skinNodeString));
    MObject commNodeObj;
    MObject skinNodeObj;
    CHECK_MSTATUS_AND_RETURN_IT(selectList.getDependNode(0, commNodeObj));
    CHECK_MSTATUS_AND_RETURN_IT(selectList.getDependNode(1, skinNodeObj));

    MFnDependencyNode commNode;
    MFnSkinCluster skinNode;
    CHECK_MSTATUS_AND_RETURN_IT(commNode.setObject(commNodeObj));
    CHECK_MSTATUS_AND_RETURN_IT(skinNode.setObject(skinNodeObj));
    
    MFnMesh bindMesh;

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
    syntax.addFlag("-dc", "-doodleConvertNode", MSyntax::kString);
    syntax.addFlag("-sc", "-doodleSkinCluster", MSyntax::kString);
    syntax.enableEdit(false);
    syntax.enableQuery(false);
    return syntax;
}

MStatus DoodleWeight::getsyntaxFlag(const MArgList& arg)
{
    MStatus status = MS::kSuccess;
    MArgDatabase argData(syntax(), arg, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    if (argData.isFlagSet("-doodleConvertNode"))
    {
        CHECK_MSTATUS_AND_RETURN_IT(argData.getFlagArgument("-doodleConvertNode", 0, commNodeString));
    }
    else
    {
        displayError("请设置转换节点");
        return MS::kFailure;
    }
    if (argData.isFlagSet("-doodleSkinCluster"))
    {
        CHECK_MSTATUS_AND_RETURN_IT(argData.getFlagArgument("-doodleSkinCluster", 0, skinNodeString));
    }
    else
    {
        displayError("请设置皮肤簇节点");
        return MS::kFailure;
    }
    return MS::kSuccess;
}
