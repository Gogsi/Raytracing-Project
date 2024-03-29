#pragma once
#include <tucano/mesh.hpp>
#include <vector>

using namespace std;
using namespace Eigen;

class Box {

public:
	Eigen::Vector3f tmax, tmin;

	vector<Tucano::Face> triangles;

	vector<Box> children;

	Box() { 
		this->tmin = Eigen::Vector3f(0, 0, 0);
		this->tmax = Eigen::Vector3f(0, 0, 0);
	};
	Box(const Box& other) {
		this->children = other.children;
		this->tmin = other.tmin;
		this->tmax = other.tmax;
		this->triangles = other.triangles;
	};

	Box& operator= (const Box& other) {
		this->children = other.children;
		this->tmin = other.tmin;
		this->tmax = other.tmax;
		this->triangles = other.triangles;
		return *this;
	};

	Box(Eigen::Vector3f tmin, Eigen::Vector3f tmax ) {
		this->tmin = tmin;
		this->tmax = tmax;
	};

	Box(Tucano::Mesh mesh) {

		float tx_max = std::numeric_limits<float>::min();
		float ty_max = std::numeric_limits<float>::min();
		float tz_max = std::numeric_limits<float>::min();
		float tx_min = std::numeric_limits<float>::max();
		float ty_min = std::numeric_limits<float>::max();
		float tz_min = std::numeric_limits<float>::max();

		for (auto i = 0; i < mesh.getNumberOfFaces(); i++)
		{
			this->triangles.push_back(mesh.getFace(i));
		}

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
	};

	bool isInBox(Vector3f origin) {
		bool x = (origin.x() >= tmin.x() && origin.x() <= tmax.x());
		bool y = (origin.y() >= tmin.y() && origin.y() <= tmax.y());
		bool z = (origin.z() >= tmin.z() && origin.z() <= tmax.z());

		return (x && y && z);
	}



};
