#include "doodleWeight.h"
#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MSyntax.h>
#include <maya/MSelectionList.h>
#include <maya/MArgDatabase.h>
#include <maya/MDagPath.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnMesh.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MItMeshVertex.h>
#include <maya/MComputation.h>

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
    //if (!(commNode.className( ) == "DoodleConvertBone"))
    //{
    //    displayError("不是doodle转换节点");
    //    return MS::kFailure;
    //}
    MItDependencyGraph iterGetMesh(skinNodeObj,
                                   MFn::kMesh,
                                   MItDependencyGraph::kDownstream,
                                   MItDependencyGraph::kDepthFirst,
                                   MItDependencyGraph::kPlugLevel,
                                   &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    MObject meshobj;//获得网格物体
    //for (iterGetMesh.reset();iterGetMesh.isDone();iterGetMesh.next())
    //{
    //    meshobj = iterGetMesh.currentItem(&status);
    //    CHECK_MSTATUS_AND_RETURN_IT(status);
    //    if (meshobj.hasFn(MFn::kMesh)) { break; };
    //}
    if (!iterGetMesh.isDone( ))
    {
        meshobj = iterGetMesh.currentItem( );
    }
    if (!meshobj.hasFn(MFn::kMesh))
    {
        displayError("找不到绑定网格");
        return MS::kFailure;
    }
    //CHECK_MSTATUS_AND_RETURN_IT(selectList.add(meshobj));
    //CHECK_MSTATUS_AND_RETURN_IT(selectList.getDagPath(2, meshPath));//获得网格dagpath

    MPlug plugWeightList = commNode.findPlug("bindWeightsList", &status);//获得绑定权重列表
    CHECK_MSTATUS_AND_RETURN_IT(status);

    MFnMesh fnmesh(meshobj);//获得网格功能集
    MDagPath meshPath;
    CHECK_MSTATUS_AND_RETURN_IT(fnmesh.getPath(meshPath));
    MItMeshVertex iterMeshVertex(meshPath);
    CHECK_MSTATUS_AND_RETURN_IT(status);


    //测试权重列表是否符合顶点列表
    if (!(plugWeightList.isArray( ) && (plugWeightList.evaluateNumElements( ) == fnmesh.numVertices( ))))
    {
        displayError("找到的网格顶点和权重列表中的顶点数量不相同");
        return MS::kFailure;
    }
    // 循环网格点复制权重

    MComputation comm;
    comm.beginComputation(true);
    comm.setProgressRange(0, 100);

    for (int index = 0; index < fnmesh.numVertices( ); index++)
    {
        if (comm.isInterruptRequested( )) return MS::kFailure;
        MIntArray jointIndex;
        MDoubleArray weightList;
        MPlug plugWeightListItem = plugWeightList.elementByLogicalIndex(index, &status);//获得权重
        CHECK_MSTATUS_AND_RETURN_IT(status);
        MPlug plugWeight = plugWeightListItem.child(0, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        //#pragma omp parallel for
        for (int i = 0; i < plugWeight.evaluateNumElements( ); i++)
        {
            double w;
            plugWeight.elementByLogicalIndex(i).getValue(w);
            jointIndex.append(i);
            weightList.append(w);
            //CHECK_MSTATUS_AND_RETURN_IT();
            //CHECK_MSTATUS_AND_RETURN_IT();
            //CHECK_MSTATUS_AND_RETURN_IT();
        }
        comm.setProgress(int(double(index) / fnmesh.numVertices( ) * 100.0));
        // 赋予权重
        CHECK_MSTATUS_AND_RETURN_IT(skinNode.setWeights(meshPath,
                                    iterMeshVertex.currentItem( ),
                                    jointIndex, weightList, false));
        iterMeshVertex.next( );
    }
    displayInfo("是否结束 ==" + MString( ) + iterMeshVertex.isDone( ));
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
    MArgDatabase argData(syntax( ), arg, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    if (argData.isFlagSet("-doodleConvertNode"))
    {
        CHECK_MSTATUS_AND_RETURN_IT(argData.getFlagArgument("-doodleConvertNode", 0, commNodeString));
    }
    else
    {
        displayError("请传入转换节点");
        return MS::kFailure;
    }
    if (argData.isFlagSet("-doodleSkinCluster"))
    {
        CHECK_MSTATUS_AND_RETURN_IT(argData.getFlagArgument("-doodleSkinCluster", 0, skinNodeString));
    }
    else
    {
        displayError("请传入皮肤簇节点");
        return MS::kFailure;
    }
    return MS::kSuccess;
}
