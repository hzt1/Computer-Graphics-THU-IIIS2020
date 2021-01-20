#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include "plane.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object3D(m) {
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
		normal = Vector3f::cross(b - a, c - a);
		material = m;
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin, bool fl = 0) override {
		if ((normal[0] != normal[0])) return false;
		Hit h = hit;
		/* if ((normal[0] != normal[0])) {
			vertices[0].print();
			vertices[1].print();
			vertices[2].print();
			cout << Vector3f::cross(b - a, c - a) << endl;
		}
		assert(normal[0] == normal[0]); */
		if (Plane(normal, Vector3f::dot(normal, vertices[0]), material).intersect(ray, h, tmin, fl) == false) {
			return false;
		}
		Vector3f p = ray.pointAtParameter(h.getT());
		/* for (int i = 0; i < 3; i ++)
			printf("%.5f %.5f %.5f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
		puts(""); */
		for (int i = 0; i < 3; i ++) {
			Vector3f tmp = Vector3f::cross(vertices[i] - p, vertices[(i + 1) % 3] - p);
			if (Vector3f::dot(tmp, normal) < 0) {
				return false;
			}
		}
		hit = h;
		
        assert(h.getT() == h.getT());
		return true;
	}

	Vector3f getCenter() {}
	
	Vector3f normal;

protected:
	Material* material;
	Vector3f vertices[3];

};

#endif //TRIANGLE_H
