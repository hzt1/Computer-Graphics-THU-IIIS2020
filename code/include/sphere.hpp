#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement functions and add more fields as necessary

#define eps (1e-3)

class Sphere : public Object3D {
public:
    Sphere() {
        // unit ball at the center
        center = Vector3f(0, 0, 0);
        radius = 1;
    }

    Sphere(const Vector3f &center, float radius, Material *material) : Object3D(material) {
        this->center = center;
        this->radius = radius;
        this->material = material;
        // 
    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin, bool fl = 0) override {
        if (!fl && material->getDiffuseColor()[0] >= 15) {
            return false;
        }
        //
        Vector3f l = center - r.getOrigin();
        // if (l.length() <= radius) return false;
        Vector3f dir = r.getDirection();
        dir.normalize();
        float tp = Vector3f::dot(l, dir);
        if (l.length() > radius + eps && tp < eps) {
            return false;
        }
        float d = sqrt(max(l.length() * l.length() - tp * tp, (float)1e-5));
        if (d >= radius - eps) {
            return false;
        }
        float td = sqrt(max(radius * radius - d * d, (float)1e-5));
        float t;
        if (l.length() <= radius + eps) {
            t = tp + td;
        } else {
            t = tp - td;
        }
        if (t < tmin || t > h.getT()) {
            return false;
        }
        Vector3f inters = r.pointAtParameter(t);
        Vector3f normal = inters - center;
        // if (Vector3f::dot(normal, dir) > 0) normal = center - inters;
        normal.normalize();
        h.set(t, material, normal);
        /* if (t != t) {
            printf("%.5f %.5f\n", tp, td);
        } */
        assert(t == t); 
        return true;
    }
    
	Photon generatePhoton(unsigned short *Xi) override {
		float theta = M_PI * (erand48(Xi) - 1.0 / 2), alpha = 2 * M_PI * erand48(Xi);
        float z = cos(theta);
        Vector3f dir(z * cos(alpha), z * sin(alpha), sin(theta));
		Photon p(center + radius * dir, dir, material->getEmissionColor() * fabs(cos(theta))); // 
        // printf("%.5f %.5f\n", material->getEmissionColor()[0], cos(theta));
        return p;
	}

    Material* getMaterial() {
        return material;
    }
    
    Vector3f getCenter() { return center; }

protected:
    Vector3f center;
    float radius;
    Material* material;

};


#endif
