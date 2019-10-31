#pragma once
#include <tucano/effects/phongmaterialshader.hpp>
#include <tucano/mesh.hpp>

#include <vector>
#include "Ray.h";
#include "hitInfo.h"
#include "box.h"
#include "Intersect.h"

#define RADIUS_CORRECTION 3.0f

struct ModelFace {
	Tucano::Face face;
	int mesh;
};

class Scene {
private :
	// list of faces (in world space)
	std::vector<ModelFace> sceneFaces;
	std::vector<Intersect::Face *> spheres;
	std::vector<Tucano::Mesh> meshes;

public:

	Scene() {

	}

	Tucano::Mesh getMesh(int index) {
		return meshes[index];
	}

	// Material is currently hardcoded for color
	Tucano::Material::Mtl addSphere(Eigen::Vector3f center, float radius, Eigen::Vector3f color) {
		spheres.push_back( new Intersect::Sphere(center, radius / RADIUS_CORRECTION));

		Tucano::Material::Mtl sphereMat = Tucano::Material::Mtl();
		sphereMat.setAmbient(Eigen::Vector3f(color) / 3);
		sphereMat.setDiffuse(color);
		sphereMat.setSpecular(Eigen::Vector3f(0.5f, 0.5f, 0.5f));
		sphereMat.setIlluminationModel(2);
		sphereMat.setShininess(10);
		sphereMat.setDissolveFactor(0);
		sphereMat.setOpticalDensity(1);

		return sphereMat;
	}

	Intersect::Face& getSphere(int index) {
		return *spheres[index];
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

