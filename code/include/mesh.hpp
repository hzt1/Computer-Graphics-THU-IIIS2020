#ifndef MESH_H
#define MESH_H

#include <vector>
#include "object3d.hpp"
#include "triangle.hpp"
#include "Vector2f.h"
#include "Vector3f.h"


class Mesh : public Object3D {

public:
    Mesh(const char *filename, Material *m, float scale, Vector3f offset, bool inv);

    struct TriangleIndex {
        TriangleIndex() {
            x[0] = 0; x[1] = 0; x[2] = 0;
        }
        int &operator[](const int i) { return x[i]; }
        // By Computer Graphics convention, counterclockwise winding is front face
        int x[3]{};
    };
    
    struct Kdtree {
        Kdtree() { lc = rc = 0; }
        int lc, rc; Vector3f mn, mx;
    };

    std::vector<Vector3f> v;
    std::vector<TriangleIndex> t;
    std::vector<Vector3f> n;
    std::vector<Kdtree> kdt;
    bool intersect(const Ray &r, Hit &h, float tmin, bool fl = 0) override;
    struct Cmp {
        bool operator()(TriangleIndex &a, TriangleIndex &b) const;

        Mesh *parent;

        Cmp() = delete;
        Cmp(Mesh *p) : parent(p) {}
    };

    Vector3f getCenter() {}

private:

    // Normal can be used for light estimation
    void computeNormal();
    int build(int l, int r, int d);
    bool qry(int x, const Ray &r, Hit &h, float tmin);
    void upd(int x, int y);
    bool checkInter(int x, const Ray &r, float tmin);
    int kdroot;
    bool inv;
};

#endif
