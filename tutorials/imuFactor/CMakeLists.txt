cmake_minimum_required(VERSION 2.6)
project(imuFactor)

find_package(GTSAM REQUIRED)
include_directories(
	${GTSAM_INCLUDE_DIR}
)

link_directories(
	${GTSAM_LIBRARY_DIRS}
)

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(imuFactor src/imuFactor.cpp)
target_link_libraries(imuFactor gtsam ${Boost_LIBRARIES})
