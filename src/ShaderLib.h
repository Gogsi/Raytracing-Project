#pragma once
#include <tucano/mesh.hpp>
#include "Flyscene.hpp"

class Shader {
public:
	static Eigen::Vector3f reflect(Eigen::Vector3f light, Eigen::Vector3f normal) {
		return  light - 2 * normal.dot(light) * normal;
	}

	static vector<Eigen::Vector3f> getNPointsOnCircle(Eigen::Vector3f center, float radius, Eigen::Vector3f normal, int n)
	{
		Eigen::Vector3f notNormal = Eigen::Vector3f(1.0, 0.0, 0.0);
		if (notNormal == normal || notNormal == -normal) {
			notNormal = Eigen::Vector3f(1.0, 1.0, 0.0).normalized();
		}

		Eigen::Vector3f radiusVector = normal.cross(notNormal).normalized() * radius;
		float theta = M_PI * 2 / n;
		Eigen::AngleAxisf rotation = Eigen::AngleAxisf(theta, normal);

		vector<Eigen::Vector3f> res;
		res.push_back(radiusVector + center);

		for (int i = 0; i < n - 1; i++) {
			radiusVector = rotation * radiusVector;
			res.push_back(radiusVector + center);
		}

		return res;
	}
};