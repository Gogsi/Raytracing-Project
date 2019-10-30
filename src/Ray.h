#pragma once
#include <tucano/mesh.hpp>

using namespace Eigen;

// class representation of a ray
class Ray {
private:
	// origin of the ray (in world space)
	Eigen::Vector3f origin;
	// direction of the ray (int world space)
	Eigen::Vector3f dir;

public:
	Ray(Eigen::Vector3f _origin, Eigen::Vector3f _dir) {
		origin = _origin;
		dir = _dir;
	}

	Eigen::Vector3f getOrigin() {
		return origin;
	}
	Eigen::Vector3f getDirection() {
		return dir;
	}

	// returns: origin + t * dir
	Vector3f getPointOnRay(float t) {
		return origin + t * dir;
	}

	Ray reflectRay(Eigen::Vector3f normal, Eigen::Vector3f point) {
		Eigen::Vector3f newDir = dir - 2 * (dir.dot(normal)) * normal;
		return Ray(point, newDir);
	}
	Ray refractRay(Eigen::Vector3f normal, Eigen::Vector3f point, Eigen::Vector3f viewDirection, float refractionIndex) {
		float dot = viewDirection.dot(normal);
		float numerator = (1 - pow(dot, 2));
		float root = sqrt(1 - (numerator / pow(refractionIndex, 2)));
		Eigen::Vector3f t = (1 / refractionIndex) * (viewDirection - dot * normal) - normal * root;
		return Ray(point, t);
	}
	std::vector<Ray> resendRay(Eigen::Vector3f normal, Eigen::Vector3f point, int n) {
		std::vector<Ray> res;
		float theta = M_PI * 2 / n;
		Eigen::AngleAxisf rotation = Eigen::AngleAxisf(theta, normal);

		for (size_t i = 0; i < n; i++)
		{
			Eigen::Vector3f newDirection = rotation * reflectRay(normal, point).getDirection();

			res.push_back(Ray(point, newDirection));
		}
		return res;
	}
};