#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>

// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
public:

    explicit Material(const Vector3f &d_color, const Vector3f &s_color = Vector3f::ZERO, int refl = 0, float s = 0) :
            diffuseColor(d_color), emissionColor(s_color), reflection(refl), shininess(s) {

    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    virtual Vector3f getEmissionColor() const {
        return emissionColor;
    }
    
    virtual int getReflection() const {
        return reflection;
    }

    float clamp(float x) { return x > 0 ? x : 0; }

    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor) {
        Vector3f shaded = Vector3f::ZERO;
        // 
        Vector3f n = hit.getNormal(), l = dirToLight, v = -ray.getDirection(), r;
        l.normalize(), v.normalize();
        r = 2 * Vector3f::dot(n, l) * n - l;
        shaded = diffuseColor * clamp(Vector3f::dot(l, n));
        shaded = shaded * lightColor;
        return shaded;
    }
    
    float shininess;

protected:
    Vector3f diffuseColor;
    Vector3f emissionColor;
    int reflection;
    // Vector3f* texture;
};


#endif // MATERIAL_H
