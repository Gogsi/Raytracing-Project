#pragma once
#include "flyscene.hpp"

namespace Intersect {

	class Face {
	protected:
		std::vector<Eigen::Vector3f> vertices;
		Eigen::Vector3f normal;
		const float precision = 0.001f;
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

	class Sphere : Face {
	public:
		HitInfo intersects(Ray& ray);
	};

	class BoxObject : public Face {
	private:
		Box box;
	public:
		virtual HitInfo intersects(Ray& ray);

		BoxObject(Box _box) :Face(vertices, normal) {
			box = _box;
		}
	};
}