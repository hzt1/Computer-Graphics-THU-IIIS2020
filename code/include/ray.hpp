#ifndef RAY_H
#define RAY_H

#include <cassert>
#include <iostream>
#include <Vector3f.h>


// Ray class mostly copied from Peter Shirley and Keith Morley
class Ray {
public:

    Ray() = delete;
    Ray(const Vector3f &orig, const Vector3f &dir) {
        origin = orig;
        direction = dir.normalized();
    }

    Ray(const Ray &r) {
        origin = r.origin;
        direction = r.direction;
    }

    const Vector3f &getOrigin() const {
        return origin;
    }

    const Vector3f &getDirection() const {
        return direction;
    }

    Vector3f pointAtParameter(float t) const {
        return origin + direction * t;
    }

private:

    Vector3f origin;
    Vector3f direction;

};

class Photon : public Ray {
public:
    Photon() = delete;
    Photon(const Vector3f &orig, const Vector3f &dir, const Vector3f col) : Ray(orig, dir) {
        origin = orig;
        direction = dir;
        color = col;
    }

    Ray ray() { return Ray(origin, direction); }

    const Vector3f &getColor() const {
        return color;
    }

    
private:
    
    Vector3f origin;
    Vector3f direction;
    Vector3f color;
};

inline std::ostream &operator<<(std::ostream &os, const Ray &r) {
    os << "Ray <" << r.getOrigin() << ", " << r.getDirection() << ">";
    return os;
}

#endif // RAY_H
