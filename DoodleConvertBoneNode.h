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

	static MStatus initialize();

	MStatus doIt(const MArgList& arge) override;

	static MSyntax createSyntax();
	//// 蒙皮权重
	//static MObject bindWeights;
	//// 子物体索引
	//static MObject subjectIndex;
	//// 参考输出
	//static MObject aimRotation;
	//// 参考转换
	//static MObject aimTranslation;
	//// 全局绑定矩阵输出
	//static MObject globalBindMatrices;
	//// 本地旋转绑定pose
	//static MObject localBindPoseRotation;
	//// 本地输出平移pose
	//static MObject localBindPoseTranslation;
private:
	DoodleDemBones DoodleConvert;
	int __bindFrame__;
	std::set<double> getFrame;
	bool IsGetFrame;
};
