#pragma once
#include <tucano/mesh.hpp>
#include <vector>
#include "Ray.h";
#include "hitInfo.h"
#include "box.h"
#include "Intersect.h"

struct ModelFace {
	Tucano::Face face;
	int mesh;
};

class Scene {
private :
	// list of faces (in world space)
	std::vector<ModelFace> sceneFaces;
	std::vector<Intersect::Sphere> spheres;
	std::vector<Tucano::Mesh> meshes;

public:

	Scene() {

	}

	Tucano::Mesh getMesh(int index) {
		return meshes[index];
	}

	void addSphere(Eigen::Vector3f center, float radius) {
		spheres.push_back( Intersect::Sphere(center, radius));
	}

	Intersect::Sphere getSphere(int index) {
		return spheres[index];
	}

	int getNumberOfSpheres() {
		return spheres.size();
	}

	void addObj(Tucano::Mesh mesh) {
		int faces = mesh.getNumberOfFaces();

		meshes.push_back(mesh);
		int meshIndex = meshes.size() - 1;

		for (size_t i = 0; i < faces; i++)
		{
			sceneFaces.push_back(ModelFace{ mesh.getFace(i), meshIndex });
		}
	}

	int getNumberOfFaces() {
		return sceneFaces.size();
	}

	ModelFace getSceneFace(int index) {
		return sceneFaces[index];
	}
};

