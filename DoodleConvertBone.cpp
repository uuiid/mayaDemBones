//-
// ==========================================================================
// Copyright 1995,2006,2008 Autodesk, Inc. All rights reserved.
//
// Use of this software is subject to the terms of the Autodesk
// license agreement provided at the time of installation or download,
// or which otherwise accompanies this software in either electronic
// or hard copy form.
// ==========================================================================
//+

#include <math.h>
#include <vector>


#include <maya/MIOStream.h>
#include <maya/MArgList.h>
#include <maya/MPxCommand.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MPoint.h>
#include <maya/MNurbsIntersector.h>
#include <maya/MDagPath.h>
#include <maya/MMatrix.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnTransform.h>
#include <maya/MVector.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnNurbsSurface.h>

#include <maya/MArgDatabase.h>
#include <maya/MAnimControl.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MComputation.h>
#include <maya/MTime.h>
#include <maya/MFnMesh.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItMeshVertex.h>
#include <maya/MDGContextGuard.h>
#include <maya/MFnComponentListData.h>
#include <maya/MFnArrayAttrsData.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnIkJoint.h>

#include <Eigen/Dense>
#include <DemBones/DemBonesExt.h>

#include "DoodleConvertBone.h"

#define DIVISION(a) if(a) cout << "====================================================================================================" <<endl;

MObject DoodleConvertBone::startFrame;
MObject DoodleConvertBone::endFrame;
MObject DoodleConvertBone::bindFrame;
MObject DoodleConvertBone::inputMesh;
MObject DoodleConvertBone::nBones;
MObject DoodleConvertBone::nInitIters;
MObject DoodleConvertBone::nIters;
MObject DoodleConvertBone::nTransIters;
MObject DoodleConvertBone::isBindUpdate;
MObject DoodleConvertBone::transAffine;
MObject DoodleConvertBone::transAffineNorm;
MObject DoodleConvertBone::nWeightsIters;
MObject DoodleConvertBone::nonZeroWeightsNum;
MObject DoodleConvertBone::weightsSmooth;
MObject DoodleConvertBone::weightsSmoothStep;

MObject DoodleConvertBone::bindWeights;
MObject DoodleConvertBone::bindWeightsList;

//MObject DoodleConvertBone::subjectIndex;
//MObject DoodleConvertBone::subjectIndex;

MObject DoodleConvertBone::localRotation;
MObject DoodleConvertBone::localRotationList;

MObject DoodleConvertBone::localTranslation;
MObject DoodleConvertBone::localTranslationList;

MObject DoodleConvertBone::globalBindMatrices;
MObject DoodleConvertBone::globalBindMatricesList;

MObject DoodleConvertBone::localBindPoseRotation;
MObject DoodleConvertBone::localBindPoseRotationList;

MObject DoodleConvertBone::localBindPoseTranslation;
MObject DoodleConvertBone::localBindPoseTranslationList;

const MTypeId DoodleConvertBone::id(2333);

// CONSTRUCTOR:
DoodleConvertBone::DoodleConvertBone()
{
    this->__bindFrame__ = 0;
    this->IsGetFrame = false;
}

// DESTRUCTOR:
DoodleConvertBone::~DoodleConvertBone()
{
}

// FOR CREATING AN INSTANCE OF THIS COMMAND:
void* DoodleConvertBone::creator()
{
    return new DoodleConvertBone;
}

MStatus DoodleConvertBone::AnalysisCommand(MDataBlock & datablock)
{
    this->_startFrame_ = datablock.inputValue(startFrame).asInt( );
    this->_endFrame_ = datablock.inputValue(endFrame).asInt( );
    this->__bindFrame__ = datablock.inputValue(bindFrame).asInt( );
    // 解算设置
    this->_nBones_ = datablock.inputValue(nBones).asInt( );
    this->_nInitIters_ = datablock.inputValue(nInitIters).asInt( );
    this->_nIters_ = datablock.inputValue(nIters).asInt( );
    this->_nTransIters_ = datablock.inputValue(nTransIters).asInt( );
    this->_isBindUpdate_ = datablock.inputValue(isBindUpdate).asInt( );
    this->_transAffine_ = datablock.inputValue(transAffine).asDouble( );
    this->_transAffineNorm_ = datablock.inputValue(transAffineNorm).asDouble( );
    this->_nWeightsIters_ = datablock.inputValue(nWeightsIters).asInt( );
    this->_nonZeroWeightsNum_ = datablock.inputValue(nonZeroWeightsNum).asInt64( );
    this->_weightsSmooth_ = datablock.inputValue(weightsSmooth).asDouble( );
    this->_weightsSmoothStep_ = datablock.inputValue(weightsSmoothStep).asDouble( );

    this->_inputmesh_ = datablock.inputValue(inputMesh).asMesh( );

    //测试开始帧和结束帧
    if (this->_startFrame_ > this->_endFrame_) {
        MGlobal::displayError("开始帧不可大于结束帧");
        return MS::kFailure;
    }
    //测试骨骼数量
    if (this->_nBones_ <= 0) {
        MGlobal::displayError("骨骼的个数不可小于零");
        return MS::kFailure;
    }
    if (datablock.inputValue(inputMesh).type() != MFnData::Type::kMesh) {
        MGlobal::displayError("不是网格类型");
        return MS::kFailure;
    }
    // 设置解算体
    this->DoodleConvert.nB = this->_nBones_;
    this->DoodleConvert.nIters = this->_nIters_;
    this->DoodleConvert.nTransIters = this->_nTransIters_;
    this->DoodleConvert.nWeightsIters = this->_nWeightsIters_;
    this->DoodleConvert.bindUpdate = this->_isBindUpdate_;
    this->DoodleConvert.transAffine = this->_transAffine_;
    this->DoodleConvert.transAffineNorm = this->_transAffineNorm_;
    this->DoodleConvert.nnz = this->_nonZeroWeightsNum_;
    this->DoodleConvert.weightsSmooth = this->_weightsSmooth_;
    this->DoodleConvert.weightsSmoothStep = this->_weightsSmoothStep_;

    this->DoodleConvert.nS = 1;
    //设置总帧数
    this->DoodleConvert.nF = this->_endFrame_ - this->_startFrame_;
    return MS::kSuccess;
}

MStatus DoodleConvertBone::InitAndCimpute( )
{
    MComputation computtation;
    computtation.beginComputation( );

    MGlobal::displayInfo("开始计算分布骨骼......请等待");
    this->DoodleConvert.init( );
    if (computtation.isInterruptRequested( )) {
        return MS::kFailure;
    }

    MGlobal::displayInfo("开始计算权重......请等待");
    this->DoodleConvert.compute( );
    computtation.endComputation( );
    return MS::kSuccess;
}


MStatus DoodleConvertBone::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    bool debug = true;
    bool treeBased = true;
    MStatus Doolstatus = MStatus::kSuccess;
    //帮助文档暂时不写
    
    CHECK_MSTATUS(this->AnalysisCommand(dataBlock));
    //获得网格体
    MFnMesh inputMesh_(this->_inputmesh_);
    this->initConvertAttr(inputMesh_);
    // 获得序列网格数据
    this->GetMeshData();

    if (this->IsGetFrame && plug == bindWeights) {
        // 设置计算需要的索引
        this->SetSubObjectIndex( );
        CHECK_MSTATUS(this->InitAndCimpute( ));
        // 获得计算输出
        this->DoodleConvert.computeRTB(0, 
                                       _localRotation_,
                                       _localTranslation_,
                                       _globalBindMatrices_,
                                       _localBindPoseRotation_,
                                       _localBindPoseTranslation_);
        // 设置输出
        MStatus status;
        MArrayDataHandle bindWeightList_ = dataBlock.outputArrayValue(bindWeightsList,&status);
        CHECK_MSTATUS(status);
        MArrayDataBuilder buildBindWeightList_ = bindWeightList_.builder(&status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        CHECK_MSTATUS_AND_RETURN_IT(buildBindWeightList_.growArray(this->DoodleConvert.w.cols( )));
        CHECK_MSTATUS_AND_RETURN_IT(bindWeightList_.set(buildBindWeightList_));
        cout << bindWeightList_.elementCount( ) << endl;
        // cout << this->DoodleConvert.w << endl;
        for (int i = 0; i < this->DoodleConvert.w.cols( ); i++)
        {
            CHECK_MSTATUS_AND_RETURN_IT(bindWeightList_.jumpToArrayElement(i));
            MArrayDataHandle hamdleBindWeight_ = bindWeightList_.outputArrayValue();
            MArrayDataBuilder bindWeight_ = hamdleBindWeight_.builder(&status);
            CHECK_MSTATUS_AND_RETURN_IT(status);
            CHECK_MSTATUS_AND_RETURN_IT(bindWeight_.growArray(this->DoodleConvert.w.rows( )));
            //MArrayDataHandle bindWeight_ = .outputArrayValue( );
            for (int index = 0; index < this->DoodleConvert.w.rows( ); index++)
            {
                MDataHandle weight = bindWeight_.addElement(index, &status);
                CHECK_MSTATUS_AND_RETURN_IT(status);
                double __w = this->DoodleConvert.w.coeff(index, i);
                // weight.setDouble(__w);
            }
            hamdleBindWeight_.set(bindWeight_);
        }
        bindWeightList_.setAllClean( );
    }
    else
    {
        if (!this->IsGetFrame)MGlobal::displayError("没有获得绑定帧, 请播放动画, 将收集动画和绑定帧");
        return MS::kFailure;
    }
    return MS::kSuccess;
}

void DoodleConvertBone::SetSubObjectIndex( )
{
    this->DoodleConvert.fStart(1) = this->DoodleConvert.fStart(0) + this->DoodleConvert.nF;
    // 添加子物体的索引
    for (int s = 0; s < this->DoodleConvert.nS; s++) {
        for (int k = this->DoodleConvert.fStart(s); k < this->DoodleConvert.fStart(s + 1); k++)
        {
            this->DoodleConvert.subjectID(k) = s;
        }
    }
}

void DoodleConvertBone::initConvertAttr(MFnMesh& inputMesh_)
{
    //设置静止的顶点数
    this->DoodleConvert.nV = inputMesh_.numVertices( );
    //设置一些全局值
    this->DoodleConvert.v.resize(3 * this->DoodleConvert.nF, this->DoodleConvert.nV);
    this->DoodleConvert.fTime.resize(this->DoodleConvert.nF);
    this->DoodleConvert.fStart.resize(this->DoodleConvert.nS + 1);
    this->DoodleConvert.fStart(0) = 0;
    this->DoodleConvert.fv.resize(inputMesh_.numFaceVertices( ));
    this->DoodleConvert.subjectID.resize(this->DoodleConvert.nF);
}

void DoodleConvertBone::GetFrameMeshData(int frame,MObject& MobjMesh)
{
    /// <summary>
    /// 获得传入帧的网格顶点数据
    /// </summary>
    /// <param name="frame"></param>
    /// <param name="MobjMesh"></param>
    this->DoodleConvert.fTime(frame) = frame;
    MItMeshVertex vexpoint(MobjMesh);
    // 使用点迭代器获得网格位置
    for (vexpoint.reset( ); !vexpoint.isDone( ); vexpoint.next( ))
    {
        int index = vexpoint.index( );
        MPoint pos = vexpoint.position(MSpace::kWorld);
        this->DoodleConvert.v.col(index).segment<3>(3 * frame) << pos.x, pos.y, pos.z;
    }
}

void DoodleConvertBone::GetBindFrame(MObject& MobjMesh)
{
    /// <summary>
    /// 获得网格中绑定帧的数据
    /// </summary>
    /// <param name="MobjMesh"></param>
    MItMeshVertex vexpointBindFrame(MobjMesh);
    MGlobal::displayInfo("已获得绑定帧");
    this->IsGetFrame = true;
    this->DoodleConvert.u.resize(this->DoodleConvert.nS * 3, this->DoodleConvert.nV);
    // 设置多边形拓扑网格;
    for (vexpointBindFrame.reset( ); !vexpointBindFrame.isDone( ); vexpointBindFrame.next( ))
    {
        int index = vexpointBindFrame.index( );
        MPoint pos = vexpointBindFrame.position(MSpace::kWorld);
        this->DoodleConvert.u.col(index).segment(0, 3) << pos.x, pos.y, pos.z;
    }
    // 获得相对于polygon obj的顶点
    MItMeshPolygon vexIter(MobjMesh);
    for (vexIter.reset( ); !vexIter.isDone( ); vexIter.next( )) {
        int index = vexIter.index( );
        MIntArray vexIndexArray;
        vexIter.getVertices(vexIndexArray);

        std::vector<int> mindex;
        for (unsigned int vexindex = 0; vexindex < vexIndexArray.length( ); vexindex++)
        {
            mindex.push_back(vexIndexArray[vexindex]);
        }
        this->DoodleConvert.fv[index] = mindex;
    }
}

void DoodleConvertBone::GetMeshData()
{
    /// <summary>
    /// 获得序列网格数据
    /// </summary>
    /// <param name="inputMeshPath"></param>
    MAnimControl aimControl;
    int currentFrame = aimControl.currentTime( ).asUnits(MTime::uiUnit( ));// 获得当前时间
    int frame = currentFrame - this->_startFrame_;// 获得内部帧
    if ((frame >= 0) && (frame < (this->_endFrame_ - this->_startFrame_))) { 
        //当帧大于零时开始获得数据 小于最小帧时取消获取
        MGlobal::displayInfo("已获得" + MString( ) + (currentFrame)+"帧数据");
        this->GetFrameMeshData(frame, this->_inputmesh_);//设置解算帧
        if (currentFrame == this->__bindFrame__) {
            GetBindFrame(this->_inputmesh_);//获得绑定帧
        }
    }
}

//void DoodleConvertBone::createJoins(const std::vector<MString>& name)
//{
//    /// <summary>
//    /// 创建骨骼
//    /// </summary>
//    /// <param name="name"></param>
//    for (std::vector<MString>::const_iterator iter = name.begin( ); iter != name.end(); ++iter)
//    {
//        MObject jointObject = this->dgModidier.createNode("joint");
//        this->dgModidier.renameNode(jointObject, *iter);
//        this->dgModidier.doIt( );
//        MFnIkJoint joint(jointObject);
//        joint.setRotationOrder(MTransformationMatrix::RotationOrder::kXYZ, true);
//        doolJoint.push_back(jointObject);
//    }
//    
//}

//void DoodleConvertBone::addCurve()
//{
//    MFnAnimCurve aim;
//    for (std::vector<MFnIkJoint>::const_iterator i = this->doolJoint.begin();
//         i != this->doolJoint.end() ; ++i)
//    {
//        MPlug plugtx = i->findPlug("tx");
//        MPlug plugty = i->findPlug("ty");
//        MPlug plugtz = i->findPlug("tz");
//        MPlug plugrx = i->findPlug("rx");
//        MPlug plugry = i->findPlug("ry");
//        MPlug plugrz = i->findPlug("rz");
//        //平移曲线
//        MObject aimTX = aim.create(plugtx, MFnAnimCurve::AnimCurveType::kAnimCurveTL, &this->dgModidier);
//        MTimeArray timeArray;
//        this->DoodleConvert.fTime.size( );
//
//        //aim.addKeys()
//        MObject aimTY = aim.create(plugty, MFnAnimCurve::AnimCurveType::kAnimCurveTL, &this->dgModidier);
//        MObject aimTZ = aim.create(plugtz, MFnAnimCurve::AnimCurveType::kAnimCurveTL, &this->dgModidier);
//        //旋转曲线
//        MObject aimRX = aim.create(plugrx, MFnAnimCurve::AnimCurveType::kAnimCurveTA, &this->dgModidier);
//        MObject aimRY = aim.create(plugry, MFnAnimCurve::AnimCurveType::kAnimCurveTA, &this->dgModidier);
//        MObject aimRZ = aim.create(plugrz, MFnAnimCurve::AnimCurveType::kAnimCurveTA, &this->dgModidier);
//        //重命名曲线
//        this->dgModidier.renameNode(aimTX, i->name( ) + "acTX");
//        this->dgModidier.renameNode(aimTX, i->name( ) + "acTY");
//        this->dgModidier.renameNode(aimTX, i->name( ) + "acTZ");
//        this->dgModidier.renameNode(aimTX, i->name( ) + "acRX");
//        this->dgModidier.renameNode(aimTX, i->name( ) + "acRY");
//        this->dgModidier.renameNode(aimTX, i->name( ) + "acRZ");
//    }
//}

MStatus DoodleConvertBone::initialize( ) {
    MStatus status;
    MFnTypedAttribute tattr;
    MFnNumericAttribute numAttr;
    // 开始帧设置
    startFrame = numAttr.create("startFrame", "sf", MFnNumericData::Type::kInt, 20, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(addAttribute(startFrame));
    // 结束帧设置
    endFrame = numAttr.create("endFrame", "ef", MFnNumericData::Type::kInt, 80, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(addAttribute(endFrame));
    //绑定帧
    bindFrame = numAttr.create("bindFrame", "bf", MFnNumericData::Type::kInt, 35, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(addAttribute(bindFrame));
    // 输入网格设置
    inputMesh = tattr.create("inputMesh", "inm", MFnData::kMesh, MObject::kNullObj, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(tattr.setWritable(true));
    CHECK_MSTATUS(tattr.setStorable(true));
    CHECK_MSTATUS(addAttribute(inputMesh));
    // 骨骼数量
    nBones = numAttr.create("nBones", "nB", MFnNumericData::Type::kInt, 60, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(addAttribute(nBones));
    // 初始集群拆分次数
    nInitIters = numAttr.create("nInitIters", "nIt", MFnNumericData::Type::kInt, 10, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(addAttribute(nInitIters));
    // 总体迭代次数
    nIters = numAttr.create("nIters", "nIters", MFnNumericData::Type::kInt64, 30, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("总体迭代次数"));
    CHECK_MSTATUS(addAttribute(nIters));
    // 次级迭代次数
    nTransIters = numAttr.create("nTransIters", "nti", MFnNumericData::Type::kInt64, 5, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("次级迭代次数"));
    CHECK_MSTATUS(addAttribute(nTransIters));
    // 更新绑定
    isBindUpdate = numAttr.create("isBindUpdate", "isB", MFnNumericData::Type::kBoolean, 0, &status);
    numAttr.setWritable(true);
    numAttr.setStorable(true);
    numAttr.setChannelBox(true);
    numAttr.setHidden(true);
    numAttr.setNiceNameOverride("更新绑定");
    CHECK_MSTATUS(addAttribute(isBindUpdate));
    //平移约束
    transAffine = numAttr.create("transAffine", "traf", MFnNumericData::Type::kDouble, 10, &status);
    numAttr.setWritable(true);
    numAttr.setStorable(true);
    numAttr.setChannelBox(true);
    numAttr.setHidden(false);
    numAttr.setNiceNameOverride("平移约束");
    CHECK_MSTATUS(addAttribute(transAffine));
    //旋转约束
    transAffineNorm = numAttr.create("transAffineNorm", "trafn", MFnNumericData::Type::kDouble, 4, &status);
    numAttr.setWritable(true);
    numAttr.setStorable(true);
    numAttr.setChannelBox(true);
    numAttr.setHidden(false);
    numAttr.setNiceNameOverride("旋转约束");
    CHECK_MSTATUS(addAttribute(transAffineNorm));
    //权重次级迭代
    nWeightsIters = numAttr.create("nWeightsIters", "nWI", MFnNumericData::Type::kInt64, 3, &status);
    numAttr.setWritable(true);
    numAttr.setStorable(true);
    numAttr.setChannelBox(true);
    numAttr.setHidden(false);
    numAttr.setNiceNameOverride("权重次级迭代");
    CHECK_MSTATUS(addAttribute(nWeightsIters));
    //权重次级迭代
    nonZeroWeightsNum = numAttr.create("nonZeroWeightsNum", "nZWN", MFnNumericData::Type::kInt64, 8, &status);
    numAttr.setWritable(true);
    numAttr.setStorable(true);
    numAttr.setChannelBox(true);
    numAttr.setHidden(false);
    numAttr.setNiceNameOverride("非零权重骨骼数");
    CHECK_MSTATUS(addAttribute(nonZeroWeightsNum));
    //权重平滑约束
    weightsSmooth = numAttr.create("weightsSmooth", "wS", MFnNumericData::Type::kDouble, 0.0001, &status);
    numAttr.setWritable(true);
    numAttr.setStorable(true);
    numAttr.setChannelBox(true);
    numAttr.setHidden(false);
    numAttr.setNiceNameOverride("权重平滑约束");
    CHECK_MSTATUS(addAttribute(weightsSmooth));
    //权重平滑步长
    weightsSmoothStep = numAttr.create("weightsSmoothStep", "wSS", MFnNumericData::Type::kDouble, 0.5, &status);
    numAttr.setWritable(true);
    numAttr.setStorable(true);
    numAttr.setChannelBox(true);
    numAttr.setHidden(false);
    numAttr.setNiceNameOverride("权重平滑步长");
    CHECK_MSTATUS(addAttribute(weightsSmoothStep));

    /// <summary>
    /// 输出属性
    /// </summary>
    /// <returns></returns>
    
    MFnCompoundAttribute comAttr;
    MFnMatrixAttribute matrixAttr;
    // 绑定权重
    bindWeights = matrixAttr.create("bindWeights", "bw", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(matrixAttr.setArray(true));
    CHECK_MSTATUS(addAttribute(bindWeights));

    bindWeightsList = comAttr.create("bindWeightsList", "bWl", &status);
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(comAttr.addChild(bindWeights));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("绑定权重"));
    CHECK_MSTATUS(addAttribute(bindWeightsList));
    
    // 动画旋转
    localRotation = matrixAttr.create("localRotation", "lr", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(matrixAttr.setArray(true));
    CHECK_MSTATUS(addAttribute(localRotation));

    localRotationList = comAttr.create("localRotationList", "lRl", &status);
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(comAttr.addChild(localRotation));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("动画旋转"));
    CHECK_MSTATUS(addAttribute(localRotationList));
    
    // 动画平移
    localTranslation = matrixAttr.create("localTranslation", "lT", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(matrixAttr.setArray(true));
    CHECK_MSTATUS(addAttribute(localTranslation));

    localTranslationList = comAttr.create("localTranslationList", "lTl", &status);
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(comAttr.addChild(localTranslation));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("动画平移"));
    CHECK_MSTATUS(addAttribute(localTranslationList));
    
    // 全局绑定矩阵
    globalBindMatrices = matrixAttr.create("globalBindMatrices", "gbm", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(matrixAttr.setArray(true));
    CHECK_MSTATUS(addAttribute(globalBindMatrices));

    globalBindMatricesList = comAttr.create("globalBindMatricesList", "gbml", &status);
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(comAttr.addChild(globalBindMatrices));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setArray(true));
    comAttr.setNiceNameOverride("全局绑定矩阵");
    CHECK_MSTATUS(addAttribute(globalBindMatricesList));
    
    // 本地绑定pose
    localBindPoseRotation = matrixAttr.create("localBindPoseRotation", "lbpr", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(matrixAttr.setArray(true));
    CHECK_MSTATUS(addAttribute(localBindPoseRotation));

    localBindPoseRotationList = comAttr.create("localBindPoseRotationList", "lprl", &status);
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true))
    CHECK_MSTATUS(comAttr.addChild(localBindPoseRotation));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("本地绑定pose"));
    CHECK_MSTATUS(addAttribute(localBindPoseRotationList));
    // 本地频移pose
    localBindPoseTranslation = matrixAttr.create("localBindPoseTranslation", "lbpt", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(matrixAttr.setArray(true));
    CHECK_MSTATUS(addAttribute(localBindPoseTranslation));

    localBindPoseTranslationList = comAttr.create("localBindPoseTranslationList", "lptl", &status);
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true))
    CHECK_MSTATUS(comAttr.addChild(localBindPoseTranslation));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("本地平移pose"));
    CHECK_MSTATUS(addAttribute(localBindPoseTranslationList));

    // 添加影响
    CHECK_MSTATUS(attributeAffects(startFrame, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(endFrame, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(inputMesh, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(nBones, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(nInitIters, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(nIters, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(nTransIters, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(isBindUpdate, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(transAffine, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(transAffineNorm, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(nWeightsIters, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(nonZeroWeightsNum, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(weightsSmooth, bindWeightsList));
    CHECK_MSTATUS(attributeAffects(weightsSmoothStep, bindWeightsList));

    return status;
}

MStatus initializePlugin(MObject obj)
{
    MStatus   status;
    MFnPlugin plugin(obj, PLUGIN_COMPANY, "0.1", "Any");

    status = plugin.registerNode("doodleConvertBone",
                                 DoodleConvertBone::id,
                                 DoodleConvertBone::creator,
                                 DoodleConvertBone::initialize);
    if (!status) {
        status.perror("registerCommand");
        return status;
    }

    return status;
}

MStatus uninitializePlugin(MObject obj)
{
    MStatus   status;
    MFnPlugin plugin(obj);

    status = plugin.deregisterNode(DoodleConvertBone::id);
    if (!status) {
        status.perror("deregisterCommand");
        return status;
    }

    return status;
}

DoodleDemBones::DoodleDemBones()
{
}

DoodleDemBones::~DoodleDemBones()
{
}
