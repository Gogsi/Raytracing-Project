#pragma once
#include <tucano/mesh.hpp>

struct HitInfo {
	float t;

	/* 
	NOTE: When using an accelerated structure, faceId contains the index to the face 
	relative the the box triangle list, 
	not the mesh triangle list
	*/
	int faceId;

	Eigen::Vector3f normal;
	Eigen::Vector3f point;
};