#ifndef PLANE_H
#define PLANE_H

#include "image.hpp"
#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() {
        this->image = nullptr;
        this->xdir = Vector3f::ZERO;
    }

    Plane(const Vector3f &normal, float d, Material *m = nullptr, Image *image = nullptr,
            Vector3f x_dir = Vector3f::ZERO, Vector3f y_dir = Vector3f::ZERO, float x_off = 0, float y_off = 0) : Object3D(m) {
        this->normal = normal;
        this->d = d;
        this->material = m;
        this->image = image;
        this->xdir = x_dir;
        this->ydir = y_dir;
        this->xoff = x_off;
        this->yoff = y_off;
    }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin, bool fl = 0) override {
        if (fabs(Vector3f::dot(normal, r.getDirection())) < 1e-4) return false;
        float t = -(-d + Vector3f::dot(normal, r.getOrigin())) / Vector3f::dot(normal, r.getDirection());
        /* if (t != t) {
            normal.print();
            printf("%.5f %.5f %.5f\n", d, Vector3f::dot(normal, r.getOrigin()), Vector3f::dot(normal, r.getDirection()));
        } */
        if (t < tmin || t > h.getT()) return false;
        Vector3f intr = r.pointAtParameter(t);
        int x = Vector3f::dot(intr, xdir) - xoff, y = Vector3f::dot(intr, ydir) - yoff;
        if (xdir != Vector3f::ZERO && 0 <= x && x < image->Width() && 0 <= y && y < image->Height()) {
            // printf("%.5f %.5f %.5f %d\n", ydir[0], ydir[1], ydir[2], image->Height());
            h.set(t, material, normal, image->GetPixel(x, y));
        } else {
            h.set(t, material, normal);
        }
        
        assert(t == t);
        return true;
    }
    
    bool intersect_box(const Ray &r, Hit &h, float tmin) {
        if (fabs(Vector3f::dot(normal, r.getDirection())) < 1e-4) return false;
        float t = -(-d + Vector3f::dot(normal, r.getOrigin())) / Vector3f::dot(normal, r.getDirection());
        h.set(t, material, normal);
        return true;
    }

    Vector3f getCenter() {}

protected:
    float d;
    Vector3f normal;
    Material* material;
    Image *image;
    Vector3f xdir, ydir;
    float xoff, yoff;

};

#endif //PLANE_H
		

