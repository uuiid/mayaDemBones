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

// MAYA HEADERS
#include <set>

#include <maya/MArgList.h>
#include <maya/MPxCommand.h>
#include <maya/MGlobal.h>
#include <maya/MTime.h>
#include <maya/MSyntax.h>
#include <maya/MDagModifier.h>

class DoodleDemBones : public Dem::DemBonesExt<double, float> {
public:
	DoodleDemBones();
	~DoodleDemBones();
};


// MAIN CLASS 
class DoodleConvertBone : public MPxCommand
{
public:
	DoodleConvertBone();
	~DoodleConvertBone() override;
	static void* creator();

	MStatus doIt(const MArgList& arge) override;

	const bool isHistoryOn( );

	const bool isUndoable( ) ;

	void SetSubObjectIndex( );

	void initConvertAttr(MFnMesh& inputMesh_);

	void GetFrameMeshData(int i, MObject& MobjMesh);

	void GetBindFrame(MObject& MobjMesh);

	void GetMeshData(MDagPath& inputMeshPath);

	void createJoins();

	void addCurve();

	void addSkinCluster( );


	MStatus AnalysisCommand(MArgList arge, MStatus Doolstatus );
	
	MStatus InitAndCimpute( );
	static MSyntax createSyntax();
	// 核心转化类实例
	DoodleDemBones DoodleConvert;
	std::vector<MObject> doolJoint;
	
	// 比较重要的值
	int startFrame;
	int endFrame;
	MString inputMesh;
	MObject bindObj;
	MString objName;
	int nBones = 30;
	MDGModifier dgModidier;
	MDagModifier dagModifier;

	// 获得设置值
	int nInitIters = 10;
	int nIters = 30;
	int nTransIters = 5;
	int isBindUpdate = 0;
	double transAffine = 10;
	double transAffineNorm = 4;
	int nWeightsIters = 3;
	int nonZeroWeightsNum = 8;
	double weightsSmooth = 0.0001;
	double weightsSmoothStep = 0.5;

	// 蒙皮权重
	//static MObject bindWeights;
	// 子物体索引
	int subjectIndex;
	// 参考输出
	Eigen::MatrixXd localRotation;
	// 参考转换
	Eigen::MatrixXd localTranslation;
	// 全局绑定矩阵输出
	Eigen::MatrixXd globalBindMatrices;
	// 本地旋转绑定pose
	Eigen::MatrixXd localBindPoseRotation;
	// 本地输出平移pose
	Eigen::MatrixXd localBindPoseTranslation;
private:
	int __bindFrame__;
	bool IsGetFrame;
};
