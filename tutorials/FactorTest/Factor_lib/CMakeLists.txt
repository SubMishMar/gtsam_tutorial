project(Factor_lib)

set(HEADER_FILES PointToPlaneFactor.h LBLFactor.h)
set(SOURCE_FILES PointToPlaneFactor.cpp LBLFactor.cpp)

find_package(GTSAM REQUIRED)

include_directories(${GTSAM_INCLUDE_DIR})
add_library(Factor_lib STATIC ${HEADER_FILES} ${SOURCE_FILES})