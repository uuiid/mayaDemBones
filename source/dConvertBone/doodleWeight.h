#pragma once

#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
class DoodleWeight:public MPxCommand
{
public:
	DoodleWeight( );
	~DoodleWeight( ) override;

	static void* creator( );
	// 创建语法对象
	static MSyntax newSyntax( );
	//工作
	MStatus doIt(const MArgList& arg) override;
	//不可撤销
	bool isUndoable() const override;
	//有语法对象
	bool hasSyntax( ) const override;

	//获得传入参数
	MStatus getsyntaxFlag(const MArgList& arg);

	//设置曲线和节点
	//MStatus create
private:
	MStatus copyWeight();
	MStatus createAimCurve(MPlug& plug,MFnIkJoint& jointNode);
	MStatus connectAimCurve( );

private:
	MString commNodeString;
	MString skinNodeString;

	MObject skinNodeObj;
	MObject commNodeObj;
};