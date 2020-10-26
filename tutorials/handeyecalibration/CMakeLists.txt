cmake_minimum_required(VERSION 2.6)
project(handeyecalibration)

find_package(GTSAM REQUIRED)
include_directories(
	${GTSAM_INCLUDE_DIR}
)

link_directories(
	${GTSAM_LIBRARY_DIRS}
)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR} ${CERES_INCLUDE_DIRS})

add_executable(handeyecalibration src/handeyecalibration.cpp)
target_link_libraries(handeyecalibration gtsam ${Boost_LIBRARIES} ${CERES_LIBRARIES})
