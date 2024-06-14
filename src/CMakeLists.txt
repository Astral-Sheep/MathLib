target_sources(${PROJECT_NAME} PRIVATE
	Math.cpp
	Matrix2x2.cpp
	Matrix3x3.cpp
	Matrix4x4.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
)

