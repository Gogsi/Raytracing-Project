#pragma once
#include <tucano/mesh.hpp>

struct HitInfo {
	float t;
	int faceId;
	Eigen::Vector3f rayOrigin;
	Eigen::Vector3f rayDirection;
};