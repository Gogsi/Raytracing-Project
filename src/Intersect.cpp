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

		Box _box = box;

		Eigen::Vector3f dir = ray.getDirection();
		Eigen::Vector3f origin = ray.getOrigin();
		Eigen::Vector3f invDir = Eigen::Vector3f(1 / dir.x(), 1 / dir.y(), 1 / dir.z());

		float tmin, tmax, tymin, tymax, tzmin, tzmax;

		if (invDir.x() >= 0) {
			tmin = (_box.tmin.x() - origin.x()) * invDir.x();
			tmax = (_box.tmax.x() - origin.x()) * invDir.x();
		}
		else {
			tmin = (_box.tmax.x() - origin.x()) * invDir.x();
			tmax = (_box.tmin.x() - origin.x()) * invDir.x();
		}

		if (invDir.y() >= 0) {
			tymin = (_box.tmin.y() - origin.y()) * invDir.y();
			tymax = (_box.tmax.y() - origin.y()) * invDir.y();
		}
		else {
			tymin = (_box.tmax.y() - origin.y()) * invDir.y();
			tymax = (_box.tmin.y() - origin.y()) * invDir.y();
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
			tzmin = (_box.tmin.z() - origin.z()) * invDir.z();
			tzmax = (_box.tmax.z() - origin.z()) * invDir.z();
		}
		else {
			tzmin = (_box.tmax.z() - origin.z()) * invDir.z();
			tzmax = (_box.tmin.z() - origin.z()) * invDir.z();
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
		// get vector perpindicular to ray
		HitInfo res = HitInfo{ INFINITY, -1, Eigen::Vector3f(0,0,0), Eigen::Vector3f(0,0,0) };
		float t0, t1;
		float radius2 = radius * radius;

		// analytic solution
		Eigen::Vector3f L = ray.getOrigin() - vertices[0];
		float a = ray.getDirection().dot(ray.getDirection());
		float b = 2 * ray.getDirection().dot(L);
		float c = L.dot(L) - radius2;
		if (solveQuadratic(a, b, c, t0, t1)) {
			if (t0 > t1) std::swap(t0, t1);

			if (t0 < 0) {
				t0 = t1; // if t0 is negative, let's use t1 instead 
			}

			if (t0 > 0 && t0 >= precision) {

				float t = t0;

				res.t = t0;
				res.faceId = -1;
				res.point = ray.getPointOnRay(t);
				res.normal = (res.point - vertices[0]).normalized();
			}
		}
		return res;
	}

	bool Sphere::solveQuadratic(float& a,  float& b, float& c, float& x0, float& x1)
	{
		float discr = b * b - 4 * a * c;
		if (discr < 0) return false;
		else if (discr == 0) x0 = x1 = -0.5 * b / a;
		else {
			float q = (b > 0) ?
				-0.5 * (b + sqrt(discr)) :
				-0.5 * (b - sqrt(discr));
			x0 = q / a;
			x1 = c / q;
		}
		if (x0 > x1) std::swap(x0, x1);

		return true;
	}
}