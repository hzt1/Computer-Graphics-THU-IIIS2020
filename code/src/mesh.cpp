#include "mesh.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <sstream>

using namespace std;

int D;

/*
bool Mesh::cmp(TriangleIndex &a, TriangleIndex &b) {
    return v[a[0]][D] < v[b[0]][D];
}
*/

bool Mesh::Cmp::operator()(TriangleIndex &a, TriangleIndex &b) const {
    return parent->v[a[0]][D] < parent->v[b[0]][D];
}

void Mesh::upd(int x, int y) {
    for (int i = 0; i < 3; i ++) {
        kdt[x].mx[i] = max(kdt[x].mx[i], kdt[y].mx[i]);
        kdt[x].mn[i] = min(kdt[x].mn[i], kdt[y].mn[i]);
    }
}

int Mesh::build(int l, int r, int d) {
    D = d;
    int x = (l + r) >> 1;
    sort(t.begin() + l, t.begin() + r + 1, Cmp(this));
    for (int i = 0; i < 3; i ++) {
        kdt[x].mx[i] = -1e30;
        kdt[x].mn[i] = 1e30;
        for (int j = 0; j < 3; j ++) {
            kdt[x].mx[i] = max(kdt[x].mx[i], v[t[x][j]][i]);
            kdt[x].mn[i] = min(kdt[x].mn[i], v[t[x][j]][i]);
        }
    }
    if (l < x) {
        kdt[x].lc = build(l, x - 1, (d + 1) % 3);
        //upd(x, kdt[x].lc);
    }
    if (x < r) {
        kdt[x].rc = build(x + 1, r, (d + 1) % 3);
        //upd(x, kdt[x].rc);
    }
    // printf("%d %d %d %d %d %d -------------------\n", x, l, r, d, kdt[x].lc, kdt[x].rc);
    /* for (int i = 0; i < 3; i ++) {
        printf("%.5f %.5f\n", kdt[x].mn[i], kdt[x].mx[i]);
    }
    puts(""); */
    if (l < x) {
        upd(x, kdt[x].lc);
    }
    if (x < r) {
        upd(x, kdt[x].rc);
    }
    return x;
}

bool Mesh::checkInter(int x, const Ray &r, float tmin) {
    Vector3f a = kdt[x].mn, b = kdt[x].mx;
    // printf("%d %.5f %.5f %.5f %.5f %.5f %.5f\n", x, a[0], a[1], a[2], b[0], b[1], b[2]);
    Hit h1, h2;
    float m1 = -1e30, m2 = 1e30;
    if (Plane(Vector3f(1, 0, 0), Vector3f::dot(Vector3f(1, 0, 0), a), material).intersect_box(r, h1, tmin)
        && Plane(Vector3f(1, 0, 0), Vector3f::dot(Vector3f(1, 0, 0), b), material).intersect_box(r, h2, tmin)) {
        m1 = max(m1, min(h1.getT(), h2.getT()));
        m2 = min(m2, max(h1.getT(), h2.getT()));
    }
    if (Plane(Vector3f(0, 1, 0), Vector3f::dot(Vector3f(0, 1, 0), a), material).intersect_box(r, h1, tmin)
        && Plane(Vector3f(0, 1, 0), Vector3f::dot(Vector3f(0, 1, 0), b), material).intersect_box(r, h2, tmin)) {
        m1 = max(m1, min(h1.getT(), h2.getT()));
        m2 = min(m2, max(h1.getT(), h2.getT()));
    }
    if (Plane(Vector3f(0, 0, 1), Vector3f::dot(Vector3f(0, 0, 1), a), material).intersect_box(r, h1, tmin)
        && Plane(Vector3f(0, 0, 1), Vector3f::dot(Vector3f(0, 0, 1), b), material).intersect_box(r, h2, tmin)) {
        m1 = max(m1, min(h1.getT(), h2.getT()));
        m2 = min(m2, max(h1.getT(), h2.getT()));
    }
    // printf("%.5f %.5f\n", m1, m2);
    return m1 < m2 + 1e-4;
}

bool Mesh::qry(int x, const Ray &r, Hit &h, float tmin) {
    if (!checkInter(x, r, tmin)) {
        // printf("%d\n", x);
        return false;
    }
    // printf("%d\n", x);
    bool result = false;
    TriangleIndex& triIndex = t[x];
    Triangle triangle(v[triIndex[0]],
                        v[triIndex[1]], v[triIndex[2]], material);
    triangle.normal = n[x];
    result |= triangle.intersect(r, h, tmin);
    if (kdt[x].lc) result |= qry(kdt[x].lc, r, h, tmin);
    if (kdt[x].rc) result |= qry(kdt[x].rc, r, h, tmin);
    return result;
}

bool Mesh::intersect(const Ray &r, Hit &h, float tmin, bool fl) {
    // Optional: Change this brute force method into a faster one.
    /* bool result = false;
    for (int triId = 0; triId < (int) t.size(); ++triId) {
        TriangleIndex& triIndex = t[triId];
        Triangle triangle(v[triIndex[0]],
                          v[triIndex[1]], v[triIndex[2]], material);
        triangle.normal = n[triId];
        result |= triangle.intersect(r, h, tmin);
    }
    return result; */
    return qry(kdroot, r, h, tmin);
}

Mesh::Mesh(const char *filename, Material *material, float scale, Vector3f offset, bool inv) : Object3D(material) {

    // Optional: Use tiny obj loader to replace this simple one.
    std::ifstream f;
    f.open(filename);
    if (!f.is_open()) {
        std::cout << "Cannot open " << filename << "\n";
        return;
    }
    std::string line;
    std::string vTok("v");
    std::string fTok("f");
    std::string texTok("vt");
    // std::string usemtl("usemtl");
    char bslash = '/', space = ' ';
    std::string tok;
    int texID;
    while (true) {
        std::getline(f, line);
        if (f.eof()) {
            break;
        }
        if (line.size() < 3) {
            continue;
        }
        if (line.at(0) == '#') {
            continue;
        }
        /* if (line.find(bslash) != std::string::npos) {
            cout << line.length() << endl;
            cout << (line[line.length() - 2] == ' ')<<endl;
        } */
        if (line.length() > 3 && line[line.length() - 2] == ' ') {
            line = line.substr(0, line.length() - 3) + line[line.length() - 1];
        }
        // while (line.length() > 3 && line[0] == ' ') line = line.substr(1, line.length());
        std::stringstream ss(line);
        ss >> tok;
        if (tok == vTok) {
            Vector3f vec;
            ss >> vec[0] >> vec[1] >> vec[2];
            vec = vec * scale + offset;
            if (inv) {
                vec[2] = 2 * offset[2] - vec[2];
            }
            v.push_back(vec);
            
            // if (vec[0] < -3.5) printf("v: %.5f %.5f %.5f %d\n", vec[0], vec[1], vec[2], v.size());
            // assert(vec[0] > -3);
        } else if (tok == fTok) {
            /* if (line.find(bslash) != std::string::npos) {
                cout << line + "abc" << endl;
            } */
            int nk = std::count(line.begin(), line.end(), space), nsl = 2;
            if (line.find(bslash) == std::string::npos) nsl = 0;
            for (int j = 0; j < line.size() - 1; j ++) {
                if (line[j] == '/' &&  line[j + 1] == '/') nsl = 1;
            }
            if (line.find(bslash) != std::string::npos && (std::count(line.begin(), line.end(), bslash) % nk != 0)) {
                printf("%d\n", line.length());
                printf("%d %d\n", nk, std::count(line.begin(), line.end(), bslash));
                while (1);
            }
            assert(std::count(line.begin(), line.end(), bslash) % nk == 0);
            std::replace(line.begin(), line.end(), bslash, space);
            std::stringstream facess(line);
            facess >> tok;
            vector<int> pv;
            for (int ii = 0; ii < nk; ii++) {
                // facess >> trig[ii] >> texID;
                int x;
                facess >> x;
                pv.push_back(x);
                /* if (x < 0 || x > 60532) {
                    printf("%d\n", x);
                }
                assert(x >= 0 && x < 60532); */
                for (int tt = 0; tt < nsl; tt ++) {
                    facess >> texID;
                }
            }
            
            for (int j = 1; j < nk - 1; j ++) {
                TriangleIndex trig;
                trig[0] = pv[0];
                trig[1] = pv[j];
                trig[2] = pv[j + 1];
                for (int ii = 0; ii < 3; ii++) {
                    trig[ii]--;
                }
                t.push_back(trig);
            }
            /* else {
                 puts("NNNNNNNNNNNN");
                TriangleIndex trig;
                for (int ii = 0; ii < 3; ii++) {
                    ss >> trig[ii];
                    trig[ii]--;
                }
                t.push_back(trig);
            } */
        } else if (tok == texTok) {
            Vector2f texcoord;
            ss >> texcoord[0];
            ss >> texcoord[1];
        }
    }
    puts("Build Starts!");
    printf("%d\n", t.size());
    printf("vsize: %d\n", v.size());

    kdt.resize(t.size());
    kdroot = build(0, t.size() - 1, 0);

    computeNormal();
    puts("Build Fin!");

    f.close();
}

void Mesh::computeNormal() {
    n.resize(t.size());
    for (int triId = 0; triId < (int) t.size(); ++triId) {
        TriangleIndex& triIndex = t[triId];
        Vector3f a = v[triIndex[1]] - v[triIndex[0]];
        Vector3f b = v[triIndex[2]] - v[triIndex[0]];
        b = Vector3f::cross(a, b);
        n[triId] = b / b.length();
    }
}
