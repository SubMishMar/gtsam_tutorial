project(Factor_lib)

set(HEADER_FILES PointToPlaneFactor.h LBLFactor.h PointToPlaneFactor2.h)
set(SOURCE_FILES PointToPlaneFactor.cpp LBLFactor.cpp PointToPlaneFactor2.cpp)

find_package(GTSAM REQUIRED)

include_directories(${GTSAM_INCLUDE_DIR})
add_library(Factor_lib STATIC ${HEADER_FILES} ${SOURCE_FILES})