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

#include <maya/MFnIkJoint.h>
#include <maya/MDataHandle.h>
#include <maya/MDagModifier.h>
#include <maya/MFnAnimCurve.h>
#include <maya/MQuaternion.h>
#include <maya/MEulerRotation.h>
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

    CHECK_MSTATUS_AND_RETURN_IT(selectList.getDependNode(0, commNodeObj));
    CHECK_MSTATUS_AND_RETURN_IT(selectList.getDependNode(1, skinNodeObj));
    // 创建动画曲线
    CHECK_MSTATUS_AND_RETURN_IT(connectAimCurve( ));
    // 复制权重
    CHECK_MSTATUS_AND_RETURN_IT(copyWeight())
    return status;
}

MStatus DoodleWeight::copyWeight()
{
    MFnDependencyNode commNode;
    MFnSkinCluster skinNode;
    CHECK_MSTATUS_AND_RETURN_IT(commNode.setObject(commNodeObj));
    CHECK_MSTATUS_AND_RETURN_IT(skinNode.setObject(skinNodeObj));

    MStatus status;
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
    //for (iterGetMesh.reset();!iterGetMesh.isDone();iterGetMesh.next())
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

MStatus DoodleWeight::createAimCurve(MPlug& plug, MFnIkJoint& jointNode)
{
    MDagModifier dagModifier;
    MFnDependencyNode commNode;
    MStatus status;
    //获得解算节点
    CHECK_MSTATUS_AND_RETURN_IT(commNode.setObject(commNodeObj));
    int statr = commNode.findPlug("startFrame", &status).asInt( );
    CHECK_MSTATUS_AND_RETURN_IT(status);
    int end = commNode.findPlug("endFrame", &status).asInt( );
    CHECK_MSTATUS_AND_RETURN_IT(status);
    //创建动画节点
    //MObject aimNodeObj
    //CHECK_MSTATUS_AND_RETURN_IT(status);
    MFnAnimCurve aim;
    MTimeArray timeArray;
    MDoubleArray TX;
    MDoubleArray TY;
    MDoubleArray TZ;
    MDoubleArray RX;
    MDoubleArray RY;
    MDoubleArray RZ;
    MDoubleArray SX;
    MDoubleArray SY;
    MDoubleArray SZ;

    for (int i = 0; i < (end - statr); i++)
    {
        CHECK_MSTATUS_AND_RETURN_IT(timeArray.append(MTime((i + statr), MTime::uiUnit( ))));
        MDataHandle HandleMatrix;
        //MMatrix jointNodeFn;
        CHECK_MSTATUS_AND_RETURN_IT(
            plug.elementByLogicalIndex(i).getValue(HandleMatrix));

        MTransformationMatrix matrix(HandleMatrix.asMatrix( ));
        MVector tVector = matrix.getTranslation(MSpace::kWorld, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        double size[3];
        CHECK_MSTATUS_AND_RETURN_IT(matrix.getScale(size, MSpace::kWorld));
        MQuaternion rQuater = matrix.rotation();
        
        CHECK_MSTATUS_AND_RETURN_IT(TX.append(tVector.x));
        CHECK_MSTATUS_AND_RETURN_IT(TY.append(tVector.y));
        CHECK_MSTATUS_AND_RETURN_IT(TZ.append(tVector.z));
        
        CHECK_MSTATUS_AND_RETURN_IT(RX.append(rQuater.asEulerRotation( ).asVector( ).x));
        CHECK_MSTATUS_AND_RETURN_IT(RY.append(rQuater.asEulerRotation( ).asVector( ).y));
        CHECK_MSTATUS_AND_RETURN_IT(RZ.append(rQuater.asEulerRotation( ).asVector( ).z));
        
        CHECK_MSTATUS_AND_RETURN_IT(SX.append(size[0]));
        CHECK_MSTATUS_AND_RETURN_IT(SY.append(size[1]));
        CHECK_MSTATUS_AND_RETURN_IT(SZ.append(size[2]));
    }
    MObject aimTx = aim.create(jointNode.findPlug("tx"), MFnAnimCurve::kAnimCurveTL, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &TX));
    
    MObject aimTy = aim.create(jointNode.findPlug("ty"), MFnAnimCurve::kAnimCurveTL, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &TY));
    
    MObject aimTz = aim.create(jointNode.findPlug("tz"), MFnAnimCurve::kAnimCurveTL, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &TZ));
    
    MObject aimRx = aim.create(jointNode.findPlug("rx"), MFnAnimCurve::kAnimCurveTA, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &RX));
    
    MObject aimRy = aim.create(jointNode.findPlug("ry"), MFnAnimCurve::kAnimCurveTA, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &RY));
    
    MObject aimRz = aim.create(jointNode.findPlug("rz"), MFnAnimCurve::kAnimCurveTA, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &RZ));
    
    MObject aimSx = aim.create(jointNode.findPlug("sx"), MFnAnimCurve::kAnimCurveTL, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &SX));
    
    MObject aimSy = aim.create(jointNode.findPlug("sy"), MFnAnimCurve::kAnimCurveTL, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &SY));
    
    MObject aimSz = aim.create(jointNode.findPlug("sz"), MFnAnimCurve::kAnimCurveTL, &dagModifier);
    CHECK_MSTATUS_AND_RETURN_IT(aim.addKeys(&MTimeArray(timeArray), &SZ));

    dagModifier.renameNode(aimTx, jointNode.name( ) + "_tx");
    dagModifier.renameNode(aimTy, jointNode.name( ) + "_ty");
    dagModifier.renameNode(aimTz, jointNode.name( ) + "_tz");
    dagModifier.renameNode(aimRx, jointNode.name( ) + "_rx");
    dagModifier.renameNode(aimRy, jointNode.name( ) + "_ry");
    dagModifier.renameNode(aimRz, jointNode.name( ) + "_rz");
    dagModifier.renameNode(aimSx, jointNode.name( ) + "_sx");
    dagModifier.renameNode(aimSy, jointNode.name( ) + "_sy");
    dagModifier.renameNode(aimSz, jointNode.name( ) + "_sz");
    dagModifier.doIt( );
    return MS::kSuccess;
}

MStatus DoodleWeight::connectAimCurve()
{
    MFnDependencyNode commNode;
    MFnSkinCluster skinNode;
    MStatus status;

    CHECK_MSTATUS_AND_RETURN_IT(commNode.setObject(commNodeObj));
    CHECK_MSTATUS_AND_RETURN_IT(skinNode.setObject(skinNodeObj));
    MItDependencyGraph iterGetMesh(skinNodeObj,
                                   MFn::kJoint,
                                   MItDependencyGraph::kUpstream,
                                   MItDependencyGraph::kDepthFirst,
                                   MItDependencyGraph::kNodeLevel,
                                   &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugWeightList = commNode.findPlug("localAnimList", &status);//获得动画数据
    MObject JointObj;
    int index = 0;
    for (iterGetMesh.reset( ); !iterGetMesh.isDone( ); iterGetMesh.next( ))
    {

        JointObj = iterGetMesh.currentItem(&status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        MFnIkJoint jointNodeFn(JointObj);
        createAimCurve(plugWeightList.elementByLogicalIndex(index).child(0), jointNodeFn);

        index++;
    }
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
