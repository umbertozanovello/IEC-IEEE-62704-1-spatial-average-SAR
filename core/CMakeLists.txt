cmake_minimum_required(VERSION 3.16)

# set the project name and version
project(SpatialAveragingSAR VERSION 0.1 LANGUAGES C)

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
    target_link_libraries(auxFuncs PUBLIC LAPACK::LAPACK)
endif()

target_link_libraries(AvgSARStep1 PUBLIC auxFuncs)
target_link_libraries(AvgSARStep2 PUBLIC auxFuncs)

# Installing #
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    include(InstallRequiredSystemLibraries)
endif()

install(FILES ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
    DESTINATION bin)
	
install(TARGETS AvgSARStep1
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include)

install(TARGETS AvgSARStep2
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/licenses/
    DESTINATION licenses)

if (WIN32 AND CMROOTS_ENABLE)
install(DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/Debug/ ${CMAKE_CURRENT_BINARY_DIR}/Release/
    DESTINATION bin
    FILE_MATCHING PATTERN "*.dll")
endif()

# Packaging #

set(CPACK_PACKAGE_NAME "SpatialAverageSAR")
set(CPACK_PACKAGE_VENDOR "Umberto Zanovello")
set(CPACK_PACKAGE_CONTACT "Umberto Zanovello")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_PACKAGE_DESCRIPTION "Python-C implementation of the spatial-average SAR alogrithm compliant to the requirements of the IEC/IEEE 62704-1 standard. ")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/umbertozanovello/IEC-IEEE-62704-1-spatial-average-SAR")
if(CMROOTS_ENABLE)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "liblapack-dev (>= 3.12.0-3)")
endif()
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/licenses/LICENSE.txt")
# set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")

if(WIN32)
    set(CPACK_GENERATOR "WIX")
	set(CPACK_WIX_VERSION 4)
elseif(UNIX)
    set(CPACK_GENERATOR DEB)
endif()

include(CPack)