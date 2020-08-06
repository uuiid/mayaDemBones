#pragma once

#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
class DoodleWeight:public MPxCommand
{
public:
	DoodleWeight( );
	~DoodleWeight( ) override;

	static void* creator( );
	// �����﷨������
	static MSyntax newSyntax( );
	//ʵ�ʹ�������
	MStatus doIt(const MArgList& arg) override;
	//ָʾ�Ƿ���Գ���
	bool isUndoable() const override;
	// ָʾ�Ƿ�����﷨ָʾ��
	bool hasSyntax( ) const override;

	// ������������
	MStatus getsyntaxFlag(const MArgList& arg);
private:
	MString commNodeString;
	MString skinNodeString;
};