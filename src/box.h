#pragma once
#include <tucano/mesh.hpp>

class Box {
public:
	Eigen::Vector3f tmax, tmin;

	Box(Tucano::Mesh mesh) {

		float tx_max = std::numeric_limits<float>::min();
		float ty_max = std::numeric_limits<float>::min();
		float tz_max = std::numeric_limits<float>::min();
		float tx_min = std::numeric_limits<float>::max();
		float ty_min = std::numeric_limits<float>::max();
		float tz_min = std::numeric_limits<float>::max();

		for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
			Eigen::Vector4f v = mesh.getVertex(i);

			float x = v.x();
			float y = v.y();
			float z = v.z();

			if (x > tx_max) {
				tx_max = x;
			}
			if (x < tx_min) {
				tx_min = x;
			}
			if (y > ty_max) {
				ty_max = y;
			}
			if (y < ty_min) {
				ty_min = y;
			}
			if (z > tz_max) {
				tz_max = z;
			}
			if (z < tz_min) {
				tz_min = z;
			}
		};



		tmax = mesh.getShapeModelMatrix() * Eigen::Vector3f(tx_max, ty_max, tz_max);
		tmin = mesh.getShapeModelMatrix() * Eigen::Vector3f(tx_min, ty_min, tz_min);

		/*std::cout << tmax << std::endl;
		std::cout << tmin << std::endl;*/

	};
};
