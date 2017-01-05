cmake_minimum_required(VERSION 2.6)
 
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
 
add_executable(test1 test1.cpp)
target_link_libraries(test1 ${GTEST_LIBRARIES} pthread)
