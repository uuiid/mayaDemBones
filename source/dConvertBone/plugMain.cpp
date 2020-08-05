#include "doodleConvertPlug.h"
#include "doodleWeight.h"

#include <maya/MFnPlugin.h>

MStatus initializePlugin(MObject obj)
{
    MStatus   status;
    MFnPlugin plugin(obj, PLUGIN_COMPANY, "0.1", "Any");
    CHECK_MSTATUS_AND_RETURN_IT(plugin.registerCommand("doodleWeight",
                                DoodleWeight::creator,
                                DoodleWeight::newSyntax));

    CHECK_MSTATUS_AND_RETURN_IT(plugin.registerNode("doodleConvertBone",
                                DoodleConvertBone::id,
                                DoodleConvertBone::creator,
                                DoodleConvertBone::initialize));
    return status;
}

MStatus uninitializePlugin(MObject obj)
{
    MStatus   status;
    MFnPlugin plugin(obj);

    CHECK_MSTATUS_AND_RETURN_IT(plugin.deregisterCommand("doodleWeight"));
    CHECK_MSTATUS_AND_RETURN_IT(plugin.deregisterNode(DoodleConvertBone::id));

    return status;
}