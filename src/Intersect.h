#pragma once
#include <tucano/mesh.hpp>
#include <vector>
#include "Ray.h";
#include "hitInfo.h"
#include "box.h"

namespace Intersect {

	class Face {
	protected:
		std::vector<Eigen::Vector3f> vertices;
		Eigen::Vector3f normal;
		const float precision = 0.001f;
		int material_id;
	public:
		virtual HitInfo intersects(Ray& ray) = 0;

		Face(std::vector<Eigen::Vector3f> _vertices, Eigen::Vector3f _normal) {
			vertices = _vertices;
			/*for (size_t i = 0; i < _vertices->size(); i++)
			{
				vertices.push_back(_vertices[i]);
			}*/
			normal = _normal;
			
		}
		Face() {
			vertices = std::vector<Eigen::Vector3f>{};
			normal = Eigen::Vector3f(0, 0, 0);
		}
		Face(Eigen::Vector3f center) {
			vertices = std::vector<Eigen::Vector3f>{center};
			normal = Eigen::Vector3f(0, 0, 0);
		}
		void setMaterialID(int index) {
			material_id = index;
		}
		int getMaterialID() {
			return material_id;
		}

		~Face(){}
	};

	class Triangle : public Face {
	public:
		virtual HitInfo intersects(Ray& ray);

		Triangle(std::vector<Eigen::Vector3f> vertices, Eigen::Vector3f normal) :Face(vertices, normal) {
			//std::cout << "kaas" << std::endl;
		}
	private :
		bool isInTriangle(Eigen::Vector3f& hit, std::vector<Eigen::Vector3f>& vertices);
	};

	class Plane : public Face {
	public:
		virtual HitInfo intersects(Ray& ray);

		Plane(std::vector<Eigen::Vector3f> vertices, Eigen::Vector3f normal) :Face(vertices, normal) {
			//std::cout << "kaas" << std::endl;
		}
	};

	class Sphere : public Face {
	private:
		float radius;
	public:
		virtual HitInfo intersects(Ray& ray);

		Sphere(Eigen::Vector3f center, float _radius) : Face(center) {
			radius = _radius;
		}
		bool solveQuadratic(float& a,  float& b,  float& c, float& x0, float& x1);
	};

	class BoxObject : public Face {
	private:
		Box& box;

	public:
		virtual HitInfo intersects(Ray& ray);

		BoxObject(Box& _box) :Face(), box(_box) { // not using the '=' operator bacause it returns a pointer causing memory leaks

		}
	};
}