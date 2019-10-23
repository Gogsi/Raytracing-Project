#pragma once
#include <tucano/mesh.hpp>

using namespace Eigen;

// class representation of a ray
class Ray {
private:
	// origin of the ray (in world space)
	Vector3f origin;
	// direction of the ray (int world space)
	Vector3f dir;

public:
	Ray(Vector3f _origin, Vector3f _dir) {
		origin = _origin;
		dir = _dir;
	}

	Vector3f getOrigin() {
		return origin;
	}
	Vector3f getDirection() {
		return dir;
	}

	// returns: origin + t * dir
	Vector3f getPointOnRay(int t) {
		return origin + t * dir;
	}

	Ray reflectRay(Vector3f normal, Vector3f point) {
		Vector3f newDir = dir - 2 * (dir.dot(normal)) * normal;
		return Ray(point, newDir);
	}
	Ray refractRay(Vector3f normal, Vector3f point) {
		// TODO
	}
};