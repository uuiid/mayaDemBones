#include <maya/MGlobal.h>
#include "doodleConvert.h"
#include <maya/MString.h>
DoodleDemBones::DoodleDemBones( )
{
}

DoodleDemBones::~DoodleDemBones( )
{
}

void DoodleDemBones::cbInitSplitBegin( )
{
    MGlobal::displayInfo("开始分割");
}

void DoodleDemBones::cbInitSplitEnd( )
{
    //MGlobal::displayInfo("分割骨骼个数 : " + nB);
}

void DoodleDemBones::cbIterBegin( )
{
    //MGlobal::displayInfo("迭代次数 # " + iter);
}

bool DoodleDemBones::cbIterEnd( )
{
    double rmseValue = rmse( );
    MGlobal::displayInfo("收敛值 = " + MString() + rmseValue);
    if (this->dolCom.isInterruptRequested( )) {
        return true;
    }
    return false;
}

void DoodleDemBones::cbWeightsBegin( )
{
    //MGlobal::displayInfo("更新权重");
}

void DoodleDemBones::cbWeightsEnd( )
{
    //MGlobal::displayInfo("更新权重结束");
}

void DoodleDemBones::cbTranformationsBegin( )
{
    //MGlobal::displayInfo("更新转换");
}

void DoodleDemBones::cbTransformationsEnd( )
{
    //MGlobal::displayInfo("更新转换结束");
}

void DoodleDemBones::cbWeightsIterBegin( )
{
    //MGlobal::displayInfo("权重迭代开始");
}

bool DoodleDemBones::cbWeightsIterEnd( )
{
    MGlobal::displayInfo("权重迭代结束");
    if (this->dolCom.isInterruptRequested( )) {
        return true;
    }
    return false;
}
