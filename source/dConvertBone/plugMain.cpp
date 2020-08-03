#include "doodleConvertPlug.h"
#include <maya/MFnPlugin.h>

MStatus initializePlugin(MObject obj)
{
    MStatus   status;
    MFnPlugin plugin(obj, PLUGIN_COMPANY, "0.1", "Any");

    status = plugin.registerNode("doodleConvertBone",
                                 DoodleConvertBone::id,
                                 DoodleConvertBone::creator,
                                 DoodleConvertBone::initialize);
    if (!status) {
        status.perror("registerCommand");
        return status;
    }

    return status;
}

MStatus uninitializePlugin(MObject obj)
{
    MStatus   status;
    MFnPlugin plugin(obj);

    status = plugin.deregisterNode(DoodleConvertBone::id);
    if (!status) {
        status.perror("deregisterCommand");
        return status;
    }

    return status;
}