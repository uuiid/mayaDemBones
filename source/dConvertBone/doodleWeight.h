#pragma once

#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
class DoodleWeight:public MPxCommand
{
public:
	DoodleWeight( );
	~DoodleWeight( ) override;

	static void* creator( );
	// 设置语法解析器
	static MSyntax newSyntax( );
	//实际工作内容
	MStatus doIt(const MArgList& arg) override;
	//指示是否可以撤销
	bool isUndoable() const override;
	// 指示是否具有语法指示器
	bool hasSyntax( ) const override;

	// 解析参数函数
	MStatus getsyntaxFlag(const MArgList& arg);
private:
	MString commNodeString;
	MString skinNodeString;
};