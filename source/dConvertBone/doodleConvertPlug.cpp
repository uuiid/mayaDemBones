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
    // ��������
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

    //���Կ�ʼ֡�ͽ���֡
    if (this->_startFrame_ > this->_endFrame_) {
        MGlobal::displayError("��ʼ֡���ɴ��ڽ���֡");
        return MS::kFailure;
    }
    //���Թ�������
    if (this->_nBones_ <= 0) {
        MGlobal::displayError("�����ĸ�������С����");
        return MS::kFailure;
    }
    if (datablock.inputValue(inputMesh).type( ) != MFnData::Type::kMesh) {
        MGlobal::displayError("������������");
        return MS::kFailure;
    }
    // ���ý�����
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
    //������֡��
    this->DoodleConvert.nF = this->_endFrame_ - this->_startFrame_;
    return MS::kSuccess;
}

MStatus DoodleConvertBone::InitAndCimpute( )
{
    this->DoodleConvert.w.resize(0, 0);
    this->DoodleConvert.m.resize(0, 0);

    MGlobal::displayInfo("��ʼ����ֲ�����......��ȴ�");
    this->DoodleConvert.dolCom.beginComputation( );
    this->DoodleConvert.init( );

    MGlobal::displayInfo("��ʼ����Ȩ��......��ȴ�");
    this->DoodleConvert.compute( );
    this->DoodleConvert.dolCom.endComputation( );
    MGlobal::displayInfo("===========�������============");
    return MS::kSuccess;
}


MStatus DoodleConvertBone::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    // ��������
    if (this->IsGetFrame && plug == getOutPut) {
        // ���ü�����Ҫ������
        this->SetSubObjectIndex( );
        CHECK_MSTATUS(this->InitAndCimpute( ));
        // ��ü������
        this->DoodleConvert.computeRTB(0,
                                       _localRotation_,
                                       _localTranslation_,
                                       _globalBindMatrices_,
                                       _localBindPoseRotation_,
                                       _localBindPoseTranslation_,
                                       false);
        // �������
        CHECK_MSTATUS_AND_RETURN_IT(this->setOutBindWeight(dataBlock));
        CHECK_MSTATUS_AND_RETURN_IT(this->setBindPose(dataBlock));
        CHECK_MSTATUS_AND_RETURN_IT(this->setAim(dataBlock));
        dataBlock.setClean(plug);
    }
    //���ÿ֡����
    else if (plug == getFrameData)
    {
        bool debug = true;
        bool treeBased = true;
        MStatus Doolstatus = MStatus::kSuccess;

        CHECK_MSTATUS(this->getInputAttr(dataBlock));
        //���������
        MFnMesh inputMesh_(this->_inputmesh_);
        this->initConvertAttr(inputMesh_);
        // ���������������
        this->GetMeshData( );
        MDataHandle handleOutFrame = dataBlock.outputValue(getFrameData);
        handleOutFrame.set(1.0);
    }
    else
    {
        if (!this->IsGetFrame)MGlobal::displayError("û�л�ð�֡, �벥�Ŷ���, ���ռ������Ͱ�֡");
        return MS::kFailure;
    }
    return MS::kSuccess;
}

void DoodleConvertBone::SetSubObjectIndex( )
{
    this->DoodleConvert.fStart(1) = this->DoodleConvert.fStart(0) + this->DoodleConvert.nF;
    // ��������������
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
    //���þ�ֹ�Ķ�����
    this->DoodleConvert.nV = inputMesh_.numVertices( );
    //����һЩȫ��ֵ
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
    /// ��ô���֡�����񶥵�����
    /// </summary>
    /// <param name="frame"></param>
    /// <param name="MobjMesh"></param>
    this->DoodleConvert.fTime(frame) = frame;
    MItMeshVertex vexpoint(MobjMesh);
    // ʹ�õ�������������λ��
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
    /// ��������а�֡������
    /// </summary>
    /// <param name="MobjMesh"></param>
    //this->DoodleConvert.fv.clear( );
    MItMeshVertex vexpointBindFrame(MobjMesh);
    MGlobal::displayInfo("�ѻ�ð�֡");
    this->IsGetFrame = true;
    this->DoodleConvert.u.resize(this->DoodleConvert.nS * 3, this->DoodleConvert.nV);
    // ���ö������������;
    for (vexpointBindFrame.reset( ); !vexpointBindFrame.isDone( ); vexpointBindFrame.next( ))
    {
        int index = vexpointBindFrame.index( );
        MPoint pos = vexpointBindFrame.position(MSpace::kWorld);
        this->DoodleConvert.u.col(index).segment(0, 3) << pos.x, pos.y, pos.z;
    }
    // ��������polygon obj�Ķ���
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
    /// ���������������
    /// </summary>
    /// <param name="inputMeshPath"></param>
    MAnimControl aimControl;
    int currentFrame = aimControl.currentTime( ).asUnits(MTime::uiUnit( ));// ��õ�ǰʱ��
    int frame = currentFrame - this->_startFrame_;// ����ڲ�֡
    if ((frame >= 0) && (frame < (this->_endFrame_ - this->_startFrame_))) {
        //��֡������ʱ��ʼ������� С����С֡ʱȡ����ȡ
        MGlobal::displayInfo("�ѻ��" + MString( ) + (currentFrame)+"֡����");
        this->GetFrameMeshData(frame, this->_inputmesh_);//���ý���֡
        if (currentFrame == this->__bindFrame__) {
            GetBindFrame(this->_inputmesh_);//��ð�֡
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
        //�������
        CHECK_MSTATUS_AND_RETURN_IT(plugBindPose.
                                    selectAncestorLogicalIndex(i, localBindPoseList));
        MDataHandle handleBindPose = plugBindPose.constructHandle(dataBlock);
        MTransformationMatrix LocalBindPose;
        Eigen::Vector3d doodleTran = this->_localBindPoseTranslation_.col(i).segment<3>(0);
        //���þ���ƽ��
        CHECK_MSTATUS_AND_RETURN_IT(
            LocalBindPose.setTranslation(
            MVector(doodleTran.x( ), doodleTran.y( ), doodleTran.z( )),
            MSpace::kWorld));
        Eigen::Vector3d doodleRon = this->_localBindPoseRotation_.col(i).segment<3>(0);
        double ro[3] = { doodleRon.x( ),doodleRon.y( ),doodleRon.z( ) };
        //���þ�����ת
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
        MArrayDataHandle arrayHandleAim(handleLocalAim, &status);//��þ��������ֱ�
        CHECK_MSTATUS_AND_RETURN_IT(status);
        MArrayDataBuilder builderAim = arrayHandleAim.builder(&status);//����������������
        CHECK_MSTATUS_AND_RETURN_IT(status);
        for (int iframe = 0; iframe < this->DoodleConvert.nF; iframe++)
        {
            // ���þ���
            MTransformationMatrix mTranLocalAim;
            Eigen::Vector3d doodleTran = this->_localTranslation_.col(ibone).segment<3>(3 * iframe);
            Eigen::Vector3d doodleRotn = this->_localRotation_.col(ibone).segment<3>(3 * iframe);
            // ������ת
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
    // ��ʼ֡����
    startFrame = numAttr.create("startFrame", "sf", MFnNumericData::Type::kInt, 20, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("��ʼ֡"))
    CHECK_MSTATUS(addAttribute(startFrame));
    
    // ����֡����
    endFrame = numAttr.create("endFrame", "ef", MFnNumericData::Type::kInt, 80, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("����֡"))
    CHECK_MSTATUS(addAttribute(endFrame));
    //��֡
    bindFrame = numAttr.create("bindFrame", "bf", MFnNumericData::Type::kInt, 35, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("��֡"))
    CHECK_MSTATUS(addAttribute(bindFrame));
    // ������������
    inputMesh = tattr.create("inputMesh", "inm", MFnData::kMesh, MObject::kNullObj, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(tattr.setWritable(true));
    CHECK_MSTATUS(tattr.setStorable(true));
    CHECK_MSTATUS(tattr.setNiceNameOverride("��������"))
    CHECK_MSTATUS(addAttribute(inputMesh));
    // ��������
    nBones = numAttr.create("nBones", "nB", MFnNumericData::Type::kInt, 50, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("����Ŀ������"))
    CHECK_MSTATUS(addAttribute(nBones));
    // ��ʼ��Ⱥ��ִ���
    nInitIters = numAttr.create("nInitIters", "nIt", MFnNumericData::Type::kInt, 10, &status);
    CHECK_MSTATUS(status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("��ʼ��ֵ�����"))
    CHECK_MSTATUS(addAttribute(nInitIters));
    // �����������
    nIters = numAttr.create("nIters", "nIters", MFnNumericData::Type::kInt64, 30, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("�����������"));
    CHECK_MSTATUS(addAttribute(nIters));
    // �μ���������
    nTransIters = numAttr.create("nTransIters", "nti", MFnNumericData::Type::kInt64, 5, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("�μ���������"));
    CHECK_MSTATUS(addAttribute(nTransIters));
    // ���°�
    isBindUpdate = numAttr.create("isBindUpdate", "isB", MFnNumericData::Type::kBoolean, 0, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(true));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("���°�"));
    CHECK_MSTATUS(addAttribute(isBindUpdate));
    //ƽ��Լ��
    transAffine = numAttr.create("transAffine", "traf", MFnNumericData::Type::kDouble, 10, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("ƽ��Լ��"));
    CHECK_MSTATUS(addAttribute(transAffine));
    //��תԼ��
    transAffineNorm = numAttr.create("transAffineNorm", "trafn", MFnNumericData::Type::kDouble, 4, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("��תԼ��"));
    CHECK_MSTATUS(addAttribute(transAffineNorm));
    //Ȩ�شμ�����
    nWeightsIters = numAttr.create("nWeightsIters", "nWI", MFnNumericData::Type::kInt64, 3, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("Ȩ�شμ�����"));
    CHECK_MSTATUS(addAttribute(nWeightsIters));
    //Ȩ�شμ�����
    nonZeroWeightsNum = numAttr.create("nonZeroWeightsNum", "nZWN", MFnNumericData::Type::kInt64, 8, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("����Ȩ�ع�����"));
    CHECK_MSTATUS(addAttribute(nonZeroWeightsNum));
    //Ȩ��ƽ��Լ��
    weightsSmooth = numAttr.create("weightsSmooth", "wS", MFnNumericData::Type::kDouble, 0.0035, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("Ȩ��ƽ��Լ��"));
    CHECK_MSTATUS(addAttribute(weightsSmooth));
    //Ȩ��ƽ������
    weightsSmoothStep = numAttr.create("weightsSmoothStep", "wSS", MFnNumericData::Type::kDouble, 0.5, &status);
    CHECK_MSTATUS(numAttr.setWritable(true));
    CHECK_MSTATUS(numAttr.setStorable(true));
    CHECK_MSTATUS(numAttr.setChannelBox(true));
    CHECK_MSTATUS(numAttr.setHidden(false));
    CHECK_MSTATUS(numAttr.setNiceNameOverride("Ȩ��ƽ������"));
    CHECK_MSTATUS(addAttribute(weightsSmoothStep));

    /// <summary>
    /// �������
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
    CHECK_MSTATUS(numAttr.setNiceNameOverride("���֡����"))
    CHECK_MSTATUS(addAttribute(getFrameData));

    MFnCompoundAttribute comAttr;
    MFnMatrixAttribute matrixAttr;
    // ��Ȩ��
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
    CHECK_MSTATUS(comAttr.setNiceNameOverride("��Ȩ��"));
    CHECK_MSTATUS(addAttribute(bindWeightsList));

    // ��������
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
    CHECK_MSTATUS(comAttr.setNiceNameOverride("��������"));
    CHECK_MSTATUS(addAttribute(localAnimList));

    // ȫ�ְ󶨾���
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
    comAttr.setNiceNameOverride("ȫ�ְ󶨾���");
    CHECK_MSTATUS(addAttribute(globalBindMatricesList));

    // ���ذ�pose
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
    CHECK_MSTATUS(comAttr.setNiceNameOverride("���ذ󶨾���"));
    CHECK_MSTATUS(addAttribute(localBindPoseList));


    // ���Ӱ��
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

    // ���Ӱ��
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