#pragma once

#include "doodleConvert.h"
#include <maya/MGlobal.h>
#include <maya/MComputation.h>
#include <maya/MPxNode.h>

class DoodleConvertBone : public MPxNode
{
public:
	DoodleConvertBone( );
	~DoodleConvertBone( ) override;
	static void* creator( );

	MStatus compute(const MPlug& plug, MDataBlock& dataBlock) override;

	static MStatus initialize( );

	void SetSubObjectIndex( );

	void initConvertAttr(MFnMesh& inputMesh_);

	void GetFrameMeshData(int i, MObject& MobjMesh);

	void GetBindFrame(MObject& MobjMesh);

	void GetMeshData( );

	MStatus setOutBindWeight(MDataBlock& dataBlock);

	MStatus setBindPose(MDataBlock& dataBlock);

	MStatus setAim(MDataBlock& dataBlock);

	MStatus getInputAttr(MDataBlock& datablock);

	MStatus InitAndCimpute( );
	// 核心转化类实例
	DoodleDemBones DoodleConvert;
	
	// 比较重要的值
	static MObject startFrame;
	static MObject endFrame;
	static MObject inputMesh;
	static MObject bindFrame;
	static MObject nBones;

	// 获得设置值
	static MObject nInitIters;
	static MObject nIters;
	static MObject nTransIters;
	static MObject isBindUpdate;
	static MObject transAffine;
	static MObject transAffineNorm;
	static MObject nWeightsIters;
	static MObject nonZeroWeightsNum;
	static MObject weightsSmooth;
	static MObject weightsSmoothStep;

	// 蒙皮权重
	static MObject bindWeights;
	static MObject bindWeightsList;
	// 获得计算
	static MObject getOutPut;
	// 拉取每帧动画
	static MObject getFrameData;
	// 动画输出
	static MObject localAnim;
	static MObject localAnimList;
	// 全局绑定矩阵输出
	static MObject globalBindMatrices;
	static MObject globalBindMatricesList;
	// 本地绑定pose
	static MObject localBindPose;
	static MObject localBindPoseList;

	// id
	static const MTypeId id;
private:

	int _startFrame_;
	int _endFrame_;
	int _bindFrame_;
	int _nBones_;
	int _nInitIters_;
	int _nIters_;
	int _nTransIters_;
	int _isBindUpdate_;
	double _transAffine_;
	double _transAffineNorm_;
	int _nWeightsIters_;
	int _nonZeroWeightsNum_;
	double _weightsSmooth_;
	double _weightsSmoothStep_;

	MObject _inputmesh_;


	// 蒙皮权重
	Eigen::MatrixXd _bindWeights_;
	// 参考输出
	Eigen::MatrixXd _localRotation_;
	// 参考转换
	Eigen::MatrixXd _localTranslation_;
	// 全局绑定矩阵输出
	Eigen::MatrixXd _globalBindMatrices_;
	// 本地旋转绑定pose
	Eigen::MatrixXd _localBindPoseRotation_;
	// 本地输出平移pose
	Eigen::MatrixXd _localBindPoseTranslation_;

	int __bindFrame__;
	bool IsGetFrame;
};
