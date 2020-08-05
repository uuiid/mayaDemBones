#include "doodleConvert.h"
#include "doodleConvertPlug.h"

#include <maya/MGlobal.h>
#include <maya/MTypeId.h>
#include <maya/MDataBlock.h>
#include <maya/MFnMesh.h>
#include <maya/MPlug.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MAnimControl.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MPxNode.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MEulerRotation.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MMatrix.h>
//#include <maya/mf>
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

MObject DoodleConvertBone::getOutPut;
MObject DoodleConvertBone::getFrameData;

MObject DoodleConvertBone::localAnim;
MObject DoodleConvertBone::localAnimList;

MObject DoodleConvertBone::globalBindMatrices;
MObject DoodleConvertBone::globalBindMatricesList;

MObject DoodleConvertBone::localBindPose;
MObject DoodleConvertBone::localBindPoseList;

const MTypeId DoodleConvertBone::id(2333);

// CONSTRUCTOR:
DoodleConvertBone::DoodleConvertBone( )
{
    this->__bindFrame__ = 0;
    this->IsGetFrame = false;
    //this->DoodleConvert = DoodleDemBones( );
}

// DESTRUCTOR:
DoodleConvertBone::~DoodleConvertBone( )
{
}

// FOR CREATING AN INSTANCE OF THIS COMMAND:
void* DoodleConvertBone::creator( )
{
    return new DoodleConvertBone;
}

MStatus DoodleConvertBone::getInputAttr(MDataBlock& datablock)
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
    if (datablock.inputValue(inputMesh).type( ) != MFnData::Type::kMesh) {
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


    MGlobal::displayInfo("开始计算分布骨骼......请等待");
    this->DoodleConvert.dolCom.beginComputation( );
    this->DoodleConvert.init( );

    MGlobal::displayInfo("开始计算权重......请等待");
    this->DoodleConvert.compute( );
    this->DoodleConvert.dolCom.endComputation( );
    MGlobal::displayInfo("===========解算完成============");
    return MS::kSuccess;
}


MStatus DoodleConvertBone::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    // 计算序列
    if (this->IsGetFrame && plug == getOutPut) {
        // 设置计算需要的索引
        this->SetSubObjectIndex( );
        CHECK_MSTATUS(this->InitAndCimpute( ));
        // 获得计算输出
        this->DoodleConvert.computeRTB(0,
                                       _localRotation_,
                                       _localTranslation_,
                                       _globalBindMatrices_,
                                       _localBindPoseRotation_,
                                       _localBindPoseTranslation_,
                                       false);
        // 设置输出
        CHECK_MSTATUS_AND_RETURN_IT(this->setOutBindWeight(dataBlock));
        CHECK_MSTATUS_AND_RETURN_IT(this->setBindPose(dataBlock));
        CHECK_MSTATUS_AND_RETURN_IT(this->setAim(dataBlock));
        dataBlock.setClean(plug);
    }
    //获得每帧数据
    else if (plug == getFrameData)
    {
        bool debug = true;
        bool treeBased = true;
        MStatus Doolstatus = MStatus::kSuccess;
        //帮助文档暂时不写

        CHECK_MSTATUS(this->getInputAttr(dataBlock));
        //获得网格体
        MFnMesh inputMesh_(this->_inputmesh_);
        this->initConvertAttr(inputMesh_);
        // 获得序列网格数据
        this->GetMeshData( );
        MDataHandle handleOutFrame = dataBlock.outputValue(getFrameData);
        handleOutFrame.set(1.0);
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
    MStatus status;
    //设置静止的顶点数
    this->DoodleConvert.nV = inputMesh_.numVertices( );
    //设置一些全局值
    this->DoodleConvert.v.resize(3 * this->DoodleConvert.nF, this->DoodleConvert.nV);
    this->DoodleConvert.fTime.resize(this->DoodleConvert.nF);
    this->DoodleConvert.fStart.resize(this->DoodleConvert.nS + 1);
    this->DoodleConvert.fStart(0) = 0;
    this->DoodleConvert.fv.resize(inputMesh_.numPolygons(&status));
    CHECK_MSTATUS(status);
    this->DoodleConvert.subjectID.resize(this->DoodleConvert.nF);
}

void DoodleConvertBone::GetFrameMeshData(int frame, MObject& MobjMesh)
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
        this->DoodleConvert.v.col(index).segment(3 * frame, 3) << pos.x, pos.y, pos.z;
    }
}

void DoodleConvertBone::GetBindFrame(MObject& MobjMesh)
{
    /// <summary>
    /// 获得网格中绑定帧的数据
    /// </summary>
    /// <param name="MobjMesh"></param>
    //this->DoodleConvert.fv.clear( );
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
    MItMeshPolygon polyIter(MobjMesh);
    for (polyIter.reset( ); !polyIter.isDone( ); polyIter.next( )) {
        int index = polyIter.index( );
        MIntArray vexIndexArray;
        polyIter.getVertices(vexIndexArray);

        std::vector<int> mindex;
        for (unsigned int vexindex = 0; vexindex < vexIndexArray.length( ); vexindex++)
        {
            mindex.push_back(vexIndexArray[vexindex]);
        }
        this->DoodleConvert.fv[index] = mindex;
    }
}

void DoodleConvertBone::GetMeshData( )
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

MStatus DoodleConvertBone::setOutBindWeight(MDataBlock& dataBlock)
{
    MStatus status;
    MPlug plugBindWeight(thisMObject( ), bindWeights);
    for (int i = 0; i < this->DoodleConvert.w.cols( ); i++)
    {
        CHECK_MSTATUS_AND_RETURN_IT(plugBindWeight.selectAncestorLogicalIndex(i, bindWeightsList));
        MDataHandle handleBindWeight = plugBindWeight.constructHandle(dataBlock);
        MArrayDataHandle arrayHandleBindWeight(handleBindWeight, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        MArrayDataBuilder arrayBuilderBindWeight = arrayHandleBindWeight.builder(&status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        for (int index = 0; index < this->DoodleConvert.w.rows( ); index++)
        {
            MDataHandle handleWeight = arrayBuilderBindWeight.addElement(index, &status);
            CHECK_MSTATUS_AND_RETURN_IT(status);
            handleWeight.setDouble(this->DoodleConvert.w.coeff(index, i));
        }
        CHECK_MSTATUS_AND_RETURN_IT(arrayHandleBindWeight.set(arrayBuilderBindWeight));
        plugBindWeight.setValue(handleBindWeight);
        plugBindWeight.destructHandle(handleBindWeight);
    }
    return MS::kSuccess;
}

MStatus DoodleConvertBone::setBindPose(MDataBlock& dataBlock)
{
    MStatus status;
    MPlug plugBindPose(thisMObject( ), localBindPose);
    for (int i = 0; i < this->DoodleConvert.nB; i++)
    {
        //获得属性
        CHECK_MSTATUS_AND_RETURN_IT(plugBindPose.
                                    selectAncestorLogicalIndex(i, localBindPoseList));
        MDataHandle handleBindPose = plugBindPose.constructHandle(dataBlock);
        MTransformationMatrix LocalBindPose;
        Eigen::Vector3d doodleTran = this->_localBindPoseTranslation_.col(i).segment<3>(0);
        //设置矩阵平移
        CHECK_MSTATUS_AND_RETURN_IT(
            LocalBindPose.setTranslation(
            MVector(doodleTran.x( ), doodleTran.y( ), doodleTran.z( )),
            MSpace::kWorld));
        Eigen::Vector3d doodleRon = this->_localBindPoseRotation_.col(i).segment<3>(0);
        double ro[3] = { doodleRon.x( ),doodleRon.y( ),doodleRon.z( ) };
        //设置矩阵旋转
        CHECK_MSTATUS_AND_RETURN_IT(LocalBindPose.setRotation(ro, MTransformationMatrix::RotationOrder::kXYZ));
        handleBindPose.setMMatrix(LocalBindPose.asMatrix( ));
        plugBindPose.setValue(handleBindPose);
        plugBindPose.destructHandle(handleBindPose);
        CHECK_MSTATUS_AND_RETURN_IT(dataBlock.setClean(plugBindPose));
    }
    CHECK_MSTATUS_AND_RETURN_IT(dataBlock.setClean(localBindPoseList));
    return MS::kSuccess;
}

MStatus DoodleConvertBone::setAim(MDataBlock& dataBlock)
{
    MStatus status;
    MPlug plugLocalAnim(thisMObject( ), localAnim);
    for (int ibone = 0; ibone < this->DoodleConvert.nB; ibone++)
    {
        CHECK_MSTATUS_AND_RETURN_IT(plugLocalAnim.selectAncestorLogicalIndex(ibone, localAnimList));
        MDataHandle handleLocalAim = plugLocalAnim.constructHandle(dataBlock);
        MArrayDataHandle arrayHandleAim(handleLocalAim, &status);//获得矩阵数组手柄
        CHECK_MSTATUS_AND_RETURN_IT(status);
        MArrayDataBuilder builderAim = arrayHandleAim.builder(&status);//获得数组矩阵生成器
        CHECK_MSTATUS_AND_RETURN_IT(status);
        for (int iframe = 0; iframe < this->DoodleConvert.nF; iframe++)
        {
            // 设置矩阵
            MTransformationMatrix mTranLocalAim;
            Eigen::Vector3d doodleTran = this->_localTranslation_.col(ibone).segment<3>(3 * iframe);
            Eigen::Vector3d doodleRotn = this->_localRotation_.col(ibone).segment<3>(3 * iframe);
            // 设置旋转
            MEulerRotation rot;
            rot.x = doodleRotn[0];
            rot.y = doodleRotn[1];
            rot.z = doodleRotn[2];
            mTranLocalAim = mTranLocalAim.rotateBy(rot, MSpace::kTransform, &status);
            CHECK_MSTATUS_AND_RETURN_IT(status);
            CHECK_MSTATUS_AND_RETURN_IT(
                mTranLocalAim.setTranslation(MVector(doodleTran[0], doodleTran[1], doodleTran[2]),
                MSpace::kWorld));
            MDataHandle handleAim = builderAim.addElement(iframe, &status);
            CHECK_MSTATUS_AND_RETURN_IT(status);
            handleAim.setMMatrix(mTranLocalAim.asMatrix( ));
        }
        CHECK_MSTATUS_AND_RETURN_IT(arrayHandleAim.set(builderAim));
        plugLocalAnim.setValue(handleLocalAim);
        plugLocalAnim.destructHandle(handleLocalAim);
    }
    return MS::kSuccess;
}

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
    CHECK_MSTATUS(numAttr.setNiceNameOverride("开始帧"))
    CHECK_MSTATUS(addAttribute(startFrame));
    
    // 结束帧设置
    endFrame = numAttr.create("endFrame", "ef", MFnNumericData::Type::kInt, 80, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("结束帧"))
    CHECK_MSTATUS(addAttribute(endFrame));
    //绑定帧
    bindFrame = numAttr.create("bindFrame", "bf", MFnNumericData::Type::kInt, 35, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("绑定帧"))
    CHECK_MSTATUS(addAttribute(bindFrame));
    // 输入网格设置
    inputMesh = tattr.create("inputMesh", "inm", MFnData::kMesh, MObject::kNullObj, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(tattr.setWritable(true));
    CHECK_MSTATUS(tattr.setStorable(true));
    CHECK_MSTATUS(tattr.setNiceNameOverride("输入网格"))
    CHECK_MSTATUS(addAttribute(inputMesh));
    // 骨骼数量
    nBones = numAttr.create("nBones", "nB", MFnNumericData::Type::kInt, 50, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("骨骼目标数量"))
    CHECK_MSTATUS(addAttribute(nBones));
    // 初始集群拆分次数
    nInitIters = numAttr.create("nInitIters", "nIt", MFnNumericData::Type::kInt, 10, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("初始拆分迭代数"))
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
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("更新绑定"));
    CHECK_MSTATUS(addAttribute(isBindUpdate));
    //平移约束
    transAffine = numAttr.create("transAffine", "traf", MFnNumericData::Type::kDouble, 10, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("平移约束"));
    CHECK_MSTATUS(addAttribute(transAffine));
    //旋转约束
    transAffineNorm = numAttr.create("transAffineNorm", "trafn", MFnNumericData::Type::kDouble, 4, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("旋转约束"));
    CHECK_MSTATUS(addAttribute(transAffineNorm));
    //权重次级迭代
    nWeightsIters = numAttr.create("nWeightsIters", "nWI", MFnNumericData::Type::kInt64, 3, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("权重次级迭代"));
    CHECK_MSTATUS(addAttribute(nWeightsIters));
    //权重次级迭代
    nonZeroWeightsNum = numAttr.create("nonZeroWeightsNum", "nZWN", MFnNumericData::Type::kInt64, 8, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("非零权重骨骼数"));
    CHECK_MSTATUS(addAttribute(nonZeroWeightsNum));
    //权重平滑约束
    weightsSmooth = numAttr.create("weightsSmooth", "wS", MFnNumericData::Type::kDouble, 0.0035, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("权重平滑约束"));
    CHECK_MSTATUS(addAttribute(weightsSmooth));
    //权重平滑步长
    weightsSmoothStep = numAttr.create("weightsSmoothStep", "wSS", MFnNumericData::Type::kDouble, 0.5, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("权重平滑步长"));
    CHECK_MSTATUS(addAttribute(weightsSmoothStep));

    /// <summary>
    /// 输出属性
    /// </summary>
    /// <returns></returns>
    getOutPut = numAttr.create("getOutPut", "out", MFnNumericData::Type::kDouble, 0.0, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(true));
    CHECK_MSTATUS(addAttribute(getOutPut));

    getFrameData = numAttr.create("getFrameData", "gfd", MFnNumericData::Type::kDouble, 0.0, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setHidden(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("获得帧序列"))
    CHECK_MSTATUS(addAttribute(getFrameData));

    MFnCompoundAttribute comAttr;
    MFnMatrixAttribute matrixAttr;
    // 绑定权重
    bindWeights = numAttr.create("bindWeights", "bw", MFnNumericData::Type::kDouble, 0.0, &status);
    CHECK_MSTATUS(numAttr.setArray(true));
    CHECK_MSTATUS(numAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(true));
    CHECK_MSTATUS(addAttribute(bindWeights));

    bindWeightsList = comAttr.create("bindWeightsList", "bWl", &status);
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.addChild(bindWeights));
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("绑定权重"));
    CHECK_MSTATUS(addAttribute(bindWeightsList));

    // 动画矩阵
    localAnim = matrixAttr.create("localAnim", "la", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setArray(true));
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setWritable(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(addAttribute(localAnim));

    localAnimList = comAttr.create("localAnimList", "lal", &status);
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.addChild(localAnim));
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("动画矩阵"));
    CHECK_MSTATUS(addAttribute(localAnimList));

    // 全局绑定矩阵
    globalBindMatrices = matrixAttr.create("globalBindMatrices", "gbm", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setWritable(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(addAttribute(globalBindMatrices));

    globalBindMatricesList = comAttr.create("globalBindMatricesList", "gbml", &status);
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.addChild(globalBindMatrices));
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(comAttr.setHidden(true));
    comAttr.setNiceNameOverride("全局绑定矩阵");
    CHECK_MSTATUS(addAttribute(globalBindMatricesList));

    // 本地绑定pose
    localBindPose = matrixAttr.create("localBindPose", "lbp", MFnMatrixAttribute::Type::kDouble, &status);
    CHECK_MSTATUS(matrixAttr.setUsesArrayDataBuilder(true));
    CHECK_MSTATUS(matrixAttr.setWritable(true));
    CHECK_MSTATUS(matrixAttr.setHidden(true));
    CHECK_MSTATUS(addAttribute(localBindPose));

    localBindPoseList = comAttr.create("localBindPoseList", "lpr", &status);
    CHECK_MSTATUS(comAttr.setArray(true));
    CHECK_MSTATUS(comAttr.addChild(localBindPose));
    CHECK_MSTATUS(comAttr.setUsesArrayDataBuilder(true))
        CHECK_MSTATUS(comAttr.setHidden(true));
    CHECK_MSTATUS(comAttr.setNiceNameOverride("本地绑定矩阵"));
    CHECK_MSTATUS(addAttribute(localBindPoseList));


    // 添加影响
    CHECK_MSTATUS(attributeAffects(startFrame, getOutPut));
    CHECK_MSTATUS(attributeAffects(endFrame, getOutPut));
    CHECK_MSTATUS(attributeAffects(inputMesh, getOutPut));
    CHECK_MSTATUS(attributeAffects(nBones, getOutPut));
    CHECK_MSTATUS(attributeAffects(nInitIters, getOutPut));
    CHECK_MSTATUS(attributeAffects(nIters, getOutPut));
    CHECK_MSTATUS(attributeAffects(nTransIters, getOutPut));
    CHECK_MSTATUS(attributeAffects(isBindUpdate, getOutPut));
    CHECK_MSTATUS(attributeAffects(transAffine, getOutPut));
    CHECK_MSTATUS(attributeAffects(transAffineNorm, getOutPut));
    CHECK_MSTATUS(attributeAffects(nWeightsIters, getOutPut));
    CHECK_MSTATUS(attributeAffects(nonZeroWeightsNum, getOutPut));
    CHECK_MSTATUS(attributeAffects(weightsSmooth, getOutPut));
    CHECK_MSTATUS(attributeAffects(weightsSmoothStep, getOutPut));

    // 添加影响
    CHECK_MSTATUS(attributeAffects(startFrame, getFrameData));
    CHECK_MSTATUS(attributeAffects(endFrame, getFrameData));
    CHECK_MSTATUS(attributeAffects(inputMesh, getFrameData));
    CHECK_MSTATUS(attributeAffects(nBones, getFrameData));
    CHECK_MSTATUS(attributeAffects(nInitIters, getFrameData));
    CHECK_MSTATUS(attributeAffects(nIters, getFrameData));
    CHECK_MSTATUS(attributeAffects(nTransIters, getFrameData));
    CHECK_MSTATUS(attributeAffects(isBindUpdate, getFrameData));
    CHECK_MSTATUS(attributeAffects(transAffine, getFrameData));
    CHECK_MSTATUS(attributeAffects(transAffineNorm, getFrameData));
    CHECK_MSTATUS(attributeAffects(nWeightsIters, getFrameData));
    CHECK_MSTATUS(attributeAffects(nonZeroWeightsNum, getFrameData));
    CHECK_MSTATUS(attributeAffects(weightsSmooth, getOutPut));
    CHECK_MSTATUS(attributeAffects(weightsSmoothStep, getOutPut));

    return status;
}