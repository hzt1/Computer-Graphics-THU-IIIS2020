#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include "plane.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>


class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up);
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point, unsigned short *Xi) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {

public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle, float foc) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        this->angle = angle;
        this->focal = foc;
    }

    Ray generateRay(const Vector2f &point, unsigned short *Xi) override {
        // 
        float z = height / 2 / tan(angle / 2);
        Vector3f dir(point - Vector2f(width, height) / 2, z);
        dir.normalize();
        Matrix3f R(horizontal, up, direction, true);
        if (focal > 1e6) return Ray(center, R * dir);
        Ray r = Ray(center, R * dir);
        Hit h;
        assert(Plane(direction, focal).intersect(r, h, 1e-2));
        Vector3f inters = r.pointAtParameter(h.getT());
        Vector3f orig = center + Vector3f(erand48(Xi), erand48(Xi), erand48(Xi)) * 4 - Vector3f(2, 2, 2);
    //    printf("%.5f %.5f %.5f\n", orig[0], orig[1], orig[2]);
    //    printf("%.5f %.5f %.5f\n", inters[0], inters[1], inters[2]);
    //    printf("-- %.5f %.5f %.5f\n", (inters - orig)[0], (inters - orig)[1], (inters - orig)[2]);
        return Ray(orig, inters - orig);
    }

protected:
    float angle;
    float focal;
};

#endif //CAMERA_H
