#pragma once
#include "Intersect.h"

namespace Intersect {

	HitInfo Triangle::intersects(Ray& ray) {

		HitInfo res = HitInfo{ INFINITY, -1, Eigen::Vector3f(0,0,0), Eigen::Vector3f(0,0,0) };

		float D = normal.dot(vertices[0]);
		float denominator = normal.dot(ray.getDirection());

		if (denominator != 0.0f) {

			float t = (D - ray.getOrigin().dot(normal)) / denominator;
			Eigen::Vector3f hitPoint = ray.getPointOnRay(t);

			if (t >= precision && isInTriangle(hitPoint, vertices)) {
				res = HitInfo{ t, 0, normal, hitPoint }; // 0 because no i is present
			}
		}
		return res;
	}
	bool Triangle::isInTriangle(Eigen::Vector3f& hit, std::vector<Eigen::Vector3f>& vertices)
	{
		Eigen::Vector3f u = vertices[1] - vertices[0];
		Eigen::Vector3f v = vertices[2] - vertices[0];
		Eigen::Vector3f w = hit - vertices[0];

		float uDotU = u.dot(u);
		float uDotV = u.dot(v);
		float vDotV = v.dot(v);

		float wDotU = w.dot(u);
		float wDotV = w.dot(v);

		float denom = uDotU * vDotV - uDotV * uDotV;
		float s1 = vDotV * wDotU - uDotV * wDotV;
		float t1 = uDotU * wDotV - uDotV * wDotU;

		float s = s1 / denom;
		float t = t1 / denom;

		if (s < 0 || t < 0 || s + t > 1) return false;
		return true;
	}

	HitInfo Plane::intersects(Ray& ray) {
		HitInfo res = HitInfo{ INFINITY, -1, Eigen::Vector3f(0,0,0), Eigen::Vector3f(0,0,0) };

		float D = normal.dot(vertices[0]);
		float denominator = normal.dot(ray.getDirection());

		if (denominator != 0.0f) {

			float t = (D - ray.getOrigin().dot(normal)) / denominator;
			Eigen::Vector3f hitPoint = ray.getPointOnRay(t);

			if (t >= precision) {
				res = HitInfo{ t, 0, normal, hitPoint }; // 0 because no i is present
			}
		}
		return res;
	}

	HitInfo BoxObject::intersects(Ray& ray) {

		Eigen::Vector3f dir = ray.getDirection();
		Eigen::Vector3f origin = ray.getOrigin();
		Eigen::Vector3f invDir = Eigen::Vector3f(1 / dir.x(), 1 / dir.y(), 1 / dir.z());

		float tmin, tmax, tymin, tymax, tzmin, tzmax;

		if (invDir.x() >= 0) {
			tmin = (box.tmin.x() - origin.x()) * invDir.x();
			tmax = (box.tmax.x() - origin.x()) * invDir.x();
		}
		else {
			tmin = (box.tmax.x() - origin.x()) * invDir.x();
			tmax = (box.tmin.x() - origin.x()) * invDir.x();
		}

		if (invDir.y() >= 0) {
			tymin = (box.tmin.y() - origin.y()) * invDir.y();
			tymax = (box.tmax.y() - origin.y()) * invDir.y();
		}
		else {
			tymin = (box.tmax.y() - origin.y()) * invDir.y();
			tymax = (box.tmin.y() - origin.y()) * invDir.y();
		}
		if ((tmin > tymax) || (tymin > tmax)) {
			return HitInfo{ INFINITY, -1 };
		}
		if (tymin > tmin) {
			tmin = tymin;
		}
		if (tymax < tmax) {
			tmax = tymax;
		}

		if (invDir.z() >= 0) {
			tzmin = (box.tmin.z() - origin.z()) * invDir.z();
			tzmax = (box.tmax.z() - origin.z()) * invDir.z();
		}
		else {
			tzmin = (box.tmax.z() - origin.z()) * invDir.z();
			tzmax = (box.tmin.z() - origin.z()) * invDir.z();
		}

		if ((tmin > tzmax) || (tzmin > tmax)) {
			return HitInfo{ INFINITY, -1 };
		}
		if (tzmin > tmin) {
			tmin = tzmin;
		}
		if (tzmax < tmax) {
			tmax = tzmax;
		}

		return HitInfo{ tmin, -1 };
	}

	HitInfo Sphere::intersects(Ray& ray) {
		//TODO
	}
}