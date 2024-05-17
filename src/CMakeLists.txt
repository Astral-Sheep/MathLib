target_sources(${PROJECT_NAME} PRIVATE
	Math.cpp
	Vector2.cpp
	Vector3.cpp
	Vector4.cpp
	Quaternion.cpp
	QuaternionD.cpp
	Matrix2x2.cpp
	Matrix3x3.cpp
	Matrix4x4.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
)

