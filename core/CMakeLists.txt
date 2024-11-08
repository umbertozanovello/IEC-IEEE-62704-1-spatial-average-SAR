cmake_minimum_required(VERSION 3.16)

# set the project name and version
project(SpatialAveragingSAR VERSION 1.0 LANGUAGES C )

##### USER DEFINABLE OPTIONS & VARIABLES #####
option(CMROOTS_ENABLE "Enable precise roots calculation" ON)
option(WARNINGS_ENABLE "Activate/deactivate warnings" OFF)
option(DEBUG_ENABLE "Activate/deactivate debug" OFF)
option(REPORT_ENABLE "Generate the IEC/IEEE report in ${REPORT_PATH} with voxel volumes defined in ${VOX_VOL}" OFF)

set(REPORT_PATH ${PROJECT_SOURCE_DIR}/iecReport_unsorted.txt)
set(VOX_VOL 1) # Voxel volumes in mm^3 used for filling the data in the IEC/IEEE report

if(CMROOTS_ENABLE)
    set(CM_ROOTS ON)
    # Find packages
    find_package(LAPACK REQUIRED)
endif()

if (WARNINGS_ENABLE)
    set(WARNINGS ON)
endif()

if (DEBUG_ENABLE)
    set(DEBUG ${REPORT_PATH})
endif()

if (REPORT_ENABLE)
    set(REPORT ${REPORT_PATH})
endif()

configure_file(${PROJECT_SOURCE_DIR}/include/globals.h.in ${PROJECT_SOURCE_DIR}/include/globals.h @ONLY)


add_library(auxFuncs STATIC auxFuncs.c)
add_library(AvgSARStep1 SHARED avgSARStep1.c)
add_library(AvgSARStep2 SHARED avgSARStep2.c)

if(CMROOTS_ENABLE)
    target_link_libraries(auxFuncs PRIVATE LAPACK::LAPACK)
endif()

target_link_libraries(AvgSARStep1 PUBLIC auxFuncs)
target_link_libraries(AvgSARStep2 PUBLIC auxFuncs)
