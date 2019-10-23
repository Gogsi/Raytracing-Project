#pragma once
#include <tucano/mesh.hpp>

struct HitInfo {
	float t;
	int faceId;
	Eigen::Vector3f normal;
	Eigen::Vector3f point;
};