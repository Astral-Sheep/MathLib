cmake_minimum_required(VERSION 3.15)
project(MathLib VERSION 1.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (BUILD_SHARED_LIBS)
	list(APPEND DEFINITIONS MATH_EXPORT)
	add_library(${PROJECT_NAME} SHARED)
else()
	add_library(${PROJECT_NAME} STATIC)
endif()

add_subdirectory(src)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Release")
	list(APPEND DEFINITIONS MATH_DEBUG)
endif()

target_compile_definitions(${PROJECT_NAME} PRIVATE ${DEFINITIONS})

