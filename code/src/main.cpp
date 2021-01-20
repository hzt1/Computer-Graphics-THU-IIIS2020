#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <queue>

#include "scene_parser.hpp"
#include "image.hpp"
#include "camera.hpp"
#include "group.hpp"
#include "light.hpp"

#include <string>

#define sqr(x) (x * x)

#define eps (1e-3)

using namespace std;

SceneParser* sceneP;

Vector3f d[1050][1050];
bool fl = 0;

Vector3f radiance(Ray &r, int depth, unsigned short *Xi) {
    double t;
    Hit hit;
    Group* baseGroup = sceneP->getGroup();
    if (!baseGroup->intersect(r, hit, eps, (depth > 0 && fl)))
        return sceneP->getBackgroundColor();
    t = hit.getT();
    Vector3f x = r.pointAtParameter(hit.getT()), n = hit.getNormal();
    Vector3f rdir = r.getDirection();
    n.normalize();
    Vector3f nl = Vector3f::dot(n, rdir) < 0 ? n : (-1) * n;
    Material *m = hit.getMaterial();
    Vector3f f = m->getDiffuseColor();
    if (hit.color[0] >= -eps) {
        f = hit.color;
    }

    Vector3f ecolor = m->getEmissionColor();
    float p = max(f[0], max(f[1], f[2]));
    // cout << depth << endl;
    if (++ depth > 5) {
        if (erand48(Xi) < p) {
            f = f * (1 / p);
        } else {
            return ecolor;
        }
    }
    if (depth > 30) return ecolor;
    int reflection = m->getReflection();
    if (reflection == 0) {
        if (fl) return ecolor;
        float r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vector3f w = nl;
        Vector3f u = Vector3f::cross((fabs(w[0]) > 0.1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), w).normalized();
        Vector3f v = Vector3f::cross(w, u);
        Vector3f d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();
        Ray emer(x, d);
        return ecolor + f * radiance(emer, depth, Xi);
    } else if (reflection == 1) {
        Ray emer(x, rdir - n * 2 * Vector3f::dot(n, rdir));
        return ecolor + f * radiance(emer, depth, Xi);
    } else {
        Ray reflRay(x, rdir - n * 2 * Vector3f::dot(n, rdir));
        bool into = Vector3f::dot(n, nl) > 0;
        float nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = Vector3f::dot(rdir, nl), cos2t;
        if ((cos2t = 1 - sqr(nnt) * (1 - sqr(ddn))) < 0)
            return ecolor + f * radiance(reflRay, depth, Xi);
        Vector3f tdir = (rdir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
        float a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : Vector3f::dot(tdir, n));
        float Re = R0 + (1 - R0) * pow(c, 5), Tr = 1 - Re, P = 0.25 + 0.5 * Re, RP = Re / P, TP = Tr / (1 - P);
        Ray emer(x, tdir);
        return ecolor + f * (depth > 2 ? (erand48(Xi) < P ?
                radiance(reflRay, depth, Xi) * RP : radiance(emer, depth, Xi) * TP):
            radiance(reflRay, depth, Xi) * Re + radiance(emer, depth, Xi) * Tr);
    }
}

float clamp(float x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
float gamma(float x) { return pow(clamp(x), 1 / 2.2); }

Vector3f c[1050][1050];

void PathTracing(Image *image, Camera *camera) {
    // 循环屏幕空间的像素
    int samps = 300;
    int width = camera->getWidth(), height = camera->getHeight();

#pragma omp parallel for schedule(dynamic, 1) 
    for (int x = 0; x < camera->getWidth(); ++x) {
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps,100.*x/(width-1));
        for (unsigned short y = 0, Xi[3] = {0, 0, x*x*x}; y < camera->getHeight(); ++y) {
            Vector3f finalColor = Vector3f::ZERO;
            for (int s = 0; s < samps; s ++) {
                double r1 = 2 * erand48(Xi), dx = r1<1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
                double r2 = 2 * erand48(Xi), dy = r1<1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
                Ray r(camera->generateRay(Vector2f(x + dx, y + dy), Xi));
                finalColor = finalColor + radiance(r, 0, Xi) * (1.0 / samps);
            }
            if (!(finalColor[0] == finalColor[0])) {
                 printf("Bad Point: %d %d\n", x, y);
            } else {
                d[x][y] = finalColor;
                finalColor = Vector3f(clamp(finalColor[0]), clamp(finalColor[1]), clamp(finalColor[2]));
                finalColor = Vector3f(gamma(finalColor[0]), gamma(finalColor[1]), gamma(finalColor[2]));
                image->SetPixel(x, y, finalColor);
            }
        }
    }
}

Vector3f centr;

struct Kdtree {
    int lc, rc; Vector3f mn, mx;
};

vector<Photon> pm;
vector<Kdtree> kdt;

int kdroot;
int Dc;

/* int K = 20, K0 = 20, K_min = 5;
int Nph = 50000;
int Cph = 6000;
float w, w0 = 1, wt = Nph / Cph / 5;
float R, Rt0 = 5; */


int K = 200, K0 = 200, K_min = 1;
int Nph = 3000000;
int Cph = 100000;
float w, w0 = 1, wt = Nph / Cph / 5;
float R, Rt0 = 0.4;

bool cmpKdt(Photon a, Photon b) { return a.getOrigin()[Dc] < b.getOrigin()[Dc]; }

void upd(int x, int y) {
    for (int i = 0; i < 3; i ++) {
        kdt[x].mx[i] = max(kdt[x].mx[i], kdt[y].mx[i]);
        kdt[x].mn[i] = min(kdt[x].mn[i], kdt[y].mn[i]);
    }
}

int kdtBuild(int l, int r, int d) {
    Dc = d;
    int x = (l + r) >> 1;
    nth_element(pm.begin() + l, pm.begin() + x, pm.begin() + r + 1, cmpKdt);
    for (int i = 0; i < 3; i ++) {
        kdt[x].mx[i] = -1e30;
        kdt[x].mn[i] = 1e30;
        for (int j = 0; j < 3; j ++) {
            kdt[x].mx[i] = max(kdt[x].mx[i], pm[x].getOrigin()[i]);
            kdt[x].mn[i] = min(kdt[x].mn[i], pm[x].getOrigin()[i]);
        }
    }
    if (l < x) {
        kdt[x].lc = kdtBuild(l, x - 1, (d + 1) % 3);
        upd(x, kdt[x].lc);
    }
    if (x < r) {
        kdt[x].rc = kdtBuild(x + 1, r, (d + 1) % 3);
        upd(x, kdt[x].rc);
    }
    return x;
}

float DisToCentr(Photon a) { return (a.getOrigin() - centr).length(); }

struct cmpKmin {
    bool operator()(Photon a, Photon b) { return (a.getOrigin() - centr).length() < (b.getOrigin() - centr).length(); }
};

Photon tp(Vector3f::ZERO, Vector3f::ZERO, Vector3f::ZERO);

float kdtDis(int x) {
    Vector3f dis;
    for (int i = 0; i < 3; i ++) {
        if (kdt[x].mn[i] <= centr[i] && centr[i] <= kdt[x].mx[i]) dis[i] = 0;
        else dis[i] = min(fabs(kdt[x].mn[i] - centr[i]), fabs(kdt[x].mx[i] - centr[i]));
    }
    return dis.length();
}

void kdtQry(int x, std::priority_queue< Photon, vector<Photon>, cmpKmin >* que) {
    if (kdtDis(x) >= DisToCentr(tp)) return ;
    if (que->size() < K || DisToCentr(pm[x]) <= DisToCentr(tp)) {
        if (que->size() >= K) {
            que->pop();
        }
        que->push(pm[x]);
        tp = que->top();
    }
    int lc = kdt[x].lc, rc = kdt[x].rc;
    if (lc && rc) {
        if (kdtDis(lc) < kdtDis(rc)) {
            kdtQry(kdt[x].lc, que);
            kdtQry(kdt[x].rc, que);
        } else {
            kdtQry(kdt[x].rc, que);
            kdtQry(kdt[x].lc, que);
        }
    } else if (lc) kdtQry(kdt[x].lc, que);
    else if (rc) kdtQry(kdt[x].rc, que);
}

void irradiance(Photon &p, int depth, unsigned short *Xi, int fc) {
    if (fc == 1) return ;
    double t;
    Hit hit;
    Group* baseGroup = sceneP->getGroup();
    if (!baseGroup->intersect(p, hit, eps))
        return ;
    t = hit.getT();
    Vector3f x = p.pointAtParameter(hit.getT()), n = hit.getNormal();
    Vector3f rdir = p.getDirection();
    n.normalize();
    Vector3f nl = Vector3f::dot(n, rdir) < 0 ? n : (-1) * n;
    Material *m = hit.getMaterial();
    Vector3f f = m->getDiffuseColor();
    if (hit.color[0] >= -eps) {
        f = hit.color;
    }
    float pf = max(f[0], max(f[1], f[2]));
    // cout << depth << endl;
    if (erand48(Xi) < pf) {
        f = f * (1 / pf);
    } else {
        return ;
    }
    if (++ depth > 30) return ;
    int reflection = m->getReflection();
    if (reflection == 0) {
        if (fc != 2) {
            Photon pn(x, rdir, p.getColor());
            pm.push_back(pn);
        }
        float r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vector3f w = nl;
        Vector3f u = Vector3f::cross((fabs(w[0]) > 0.1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), w).normalized();
        Vector3f v = Vector3f::cross(w, u);
        Vector3f d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();
        Photon emer(x, d, p.getColor() * f);
        irradiance(emer, depth, Xi, fc);
        return ;
    } else if (reflection == 1) {
        Photon emer(x, rdir - n * 2 * Vector3f::dot(n, rdir), p.getColor() * f);
        irradiance(emer, depth, Xi, fc | 1);
        return ;
    } else {
        Photon reflp(x, rdir - n * 2 * Vector3f::dot(n, rdir), p.getColor() * f);
        bool into = Vector3f::dot(n, nl) > 0;
        float nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = Vector3f::dot(rdir, nl), cos2t;
        if ((cos2t = 1 - sqr(nnt) * (1 - sqr(ddn))) < 0) {
            irradiance(reflp, depth, Xi, fc | 1);
            return ;
        }
        Vector3f tdir = (rdir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
        float a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : Vector3f::dot(tdir, n));
        float Re = R0 + (1 - R0) * pow(c, 5), Tr = 1 - Re, P = 0.25 + 0.5 * Re, RP = Re / P, TP = Tr / (1 - P);
        Photon emer(x, tdir, p.getColor() * f);
        erand48(Xi) < P ? irradiance(reflp, depth, Xi, fc | 1) : irradiance(emer, depth, Xi, fc | 1);
    }
}

// bool cmp(Photon a, Photon b) { return (a.getOrigin() - centr).length() < (b.getOrigin() - centr).length(); }

Vector3f collectPhoton(Ray &r, int depth, unsigned short *Xi, float pw) {
    double t;
    Hit hit;
    Group* baseGroup = sceneP->getGroup();
    if (!baseGroup->intersect(r, hit, eps, (depth > 0 && fl)))
        return sceneP->getBackgroundColor();
    t = hit.getT();
    assert (t == t);
    Vector3f x = r.pointAtParameter(hit.getT()), n = hit.getNormal();
    Vector3f rdir = r.getDirection();
    n.normalize();
    Vector3f nl = Vector3f::dot(n, rdir) < 0 ? n : (-1) * n;
    Material *m = hit.getMaterial();
    Vector3f f = m->getDiffuseColor();
    Vector3f ecolor = m->getEmissionColor();

    Group* group = sceneP->getGroup();
    if (m->shininess > 1e-3) {
        for (auto it = group->list.begin(); it != group->list.end(); it ++) {
            Object3D* obj = (*it);
            Vector3f emission = obj->getMaterial()->getEmissionColor();
                // printf("%.5f\n", emission[0]);
            if (fabs(Vector3f::dot(emission, Vector3f(1, 1, 1))) > 1e-3) {
                Vector3f l = (obj->getCenter() - x).normalized(), rn;
                rn = 2 * Vector3f::dot(nl, l) * nl - l;
                ecolor += f * pow(clamp(Vector3f::dot((-rdir).normalized(), rn)), m->shininess) * 0.3;
            }
            //    printf("pm size: %d\n", pm.size());
        }
    }
    if (hit.color[0] >= -eps) {
        f = hit.color;
    }
    float p = max(f[0], max(f[1], f[2]));
    // cout << depth << endl;
    // ecolor = Vector3f::ZERO;
    if (++ depth > 5) {
        if (depth > 30) return ecolor;
        if (pw < 1e-5) return ecolor;
        /* if (erand48(Xi) < p) {
            f = f * (1 / p);
        } else {
            return ecolor;
        } */
    }
    int reflection = m->getReflection();
    if (reflection == 0) {
        Vector3f col = Vector3f::ZERO;
        if (fl) return col;
        centr = x;
        priority_queue< Photon, vector<Photon>, cmpKmin >* que;
        que = new(priority_queue< Photon, vector<Photon>, cmpKmin >);
        assert( que->empty() );
        que->push(Photon(centr + Vector3f(R, 0, 0), Vector3f::ZERO, Vector3f(-0.1, 0, 0)));
        tp = que->top();
        kdtQry(kdroot, que);
        if (que->top().getColor()[0] < 0) que->pop();
        if (que->size() < K_min && R == Rt0) {
            K = K_min * 2;
            que->push(Photon(centr + Vector3f(R * 1e6, 0, 0), Vector3f::ZERO, Vector3f(-0.1, 0, 0)));
            tp = que->top();
            kdtQry(kdroot, que);
            assert(que->size() == K);
            K = K0;
        }
        if (!que->empty() && que->top().getColor()[0] < 0) que->pop();
        if (que->empty()) {
            assert(0);
            return ecolor;
        }
        float rd = (que->top().getOrigin() - x).length();
        if (rd != rd) {
            printf("TopOrigin: %.5f %.5f %.5f\n", que->top().getOrigin()[0], que->top().getOrigin()[1], que->top().getOrigin()[2]);
            printf("TopOrigin: %.5f %.5f %.5f\n", x[0], x[1], x[2]);
        }
        if (rd <= 1e-3) {
            cout << que->size()<< endl;
            while (!que->empty()) {
                que->pop();
            }
            K = K_min * 2;
            que->push(Photon(centr + Vector3f(R * 1e6, 0, 0), Vector3f::ZERO, Vector3f(-0.1, 0, 0)));
            tp = que->top();
            fl = 1;
            kdtQry(kdroot, que);
            fl = 0;
            K = K0;
            while (1);
        }
        // printf("%d %.5f\n", que->size(), rd);
        int z = que->size();
        while (!que->empty()) {
            if (que->top().getColor()[0] >= 0)
                col += que->top().getColor();
            // else z ++;
            que->pop();
        }
        assert(z);
        assert(rd == rd);
        if (rd <= 1e-3) printf("%.5f\n", rd);
        assert(rd > 1e-4);
        col = col * w / (rd * rd * (Nph + Cph * wt)) * 1000 * K / z;
        delete que;
        return ecolor + f * col;
    } else if (reflection == 1) {
        Ray emer(x, rdir - n * 2 * Vector3f::dot(n, rdir));
        return ecolor + f * collectPhoton(emer, depth, Xi, pw);
    } else {
        Ray reflRay(x, rdir - n * 2 * Vector3f::dot(n, rdir));
        bool into = Vector3f::dot(n, nl) > 0;
        float nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = Vector3f::dot(rdir, nl), cos2t;
        if ((cos2t = 1 - sqr(nnt) * (1 - sqr(ddn))) < 0)
            return ecolor + f * collectPhoton(reflRay, depth, Xi, pw);
        Vector3f tdir = (rdir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
        float a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : Vector3f::dot(tdir, n));
        float Re = R0 + (1 - R0) * pow(c, 5), Tr = 1 - Re, P = 0.25 + 0.5 * Re, RP = Re / P, TP = Tr / (1 - P);
        Ray emer(x, tdir);
        return ecolor + f * (pw < 0.01 ? (erand48(Xi) < P ?
                collectPhoton(reflRay, depth, Xi, pw) * RP : collectPhoton(emer, depth, Xi, pw) * TP) :
                (collectPhoton(reflRay, depth, Xi, pw * Re) * Re + collectPhoton(emer, depth, Xi, pw * Tr) * Tr) );
        //
            
    }
}

void RayTracing(Image *image, Camera *camera, int srd) {
    // 循环屏幕空间的像素
    int samps = fl ? 20 : 1;
    int width = camera->getWidth(), height = camera->getHeight();

 // #pragma omp parallel for collapse(1) schedule(dynamic, 1) private(centr)
    for (int x = 0; x < camera->getWidth(); ++x) {
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps,100.*x/(width-1));
        for (unsigned short y = 0, Xi[3] = {0, 0, x*x*x*srd*srd}; y < camera->getHeight(); ++y) {
            Vector3f finalColor = Vector3f::ZERO;
            for (int s = 0; s < samps; s ++) {
                double r1 = 2 * erand48(Xi), dx = r1<1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
                double r2 = 2 * erand48(Xi), dy = r1<1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
                Ray r(camera->generateRay(Vector2f(x + dx, y + dy), Xi));
                finalColor = finalColor + collectPhoton(r, 0, Xi, 1.0) * (1.0 / samps);
            }
            if (!(finalColor[0] == finalColor[0])) {
                printf("Bad Point: %d %d\n", x, y);
            } else {
                c[x][y] += finalColor;
            }
        }
    } 
    
}

void PhotonMapping(int fc, int srd) {
    pm.clear();
    kdt.clear();
    int nn = fc ? Cph : Nph;
    unsigned short Xi[3] = {0, 0, 233 * srd * srd};
    Group* group = sceneP->getGroup();
    printf("%d\n", group->getGroupSize());
    int lst = 0;
    while(pm.size() < nn) {
        for (auto it = group->list.begin(); it != group->list.end(); it ++) {
            Object3D* obj = (*it);
            Vector3f emission = obj->getMaterial()->getEmissionColor();
            if (fabs(Vector3f::dot(emission, Vector3f(1, 1, 1))) > 1e-3 && fabs(emission[0] - 12) < 1e-3) {
            // for (int i = 1; i <= nn; i ++) {
                 if (pm.size() % 5000 == 0 && pm.size() != lst) printf("%d\n", pm.size()), lst = pm.size();
                Photon p(obj->generatePhoton(Xi));
                irradiance(p, 0, Xi, fc);
            }
        }
        //    printf("pm size: %d\n", pm.size());
    }
    if (pm.size() == 0) return ;
    kdt.resize(pm.size());
    kdroot = kdtBuild(0, pm.size() - 1, 0);
    printf("pm size: %d\n", pm.size());
}

void PPM(Image *image, Camera *camera, char *outputFile) {
    for (int i = 1;; i ++) {
        if (Rt0 > 0.2) Rt0 *= 0.8;
        else if (Rt0 > 0.1) Rt0 *= 0.9;
        else if (Rt0 > 0.05) Rt0 *= 0.95;
        else Rt0 *= 0.98;
        printf("Round: %d %.5f\n", i, Rt0);
        for (int x = 0; x < camera->getWidth(); ++x) {
            for (int y = 0; y < camera->getHeight(); ++y) {
                c[x][y] = c[x][y] * (i - 1);
            }
        }
        R = Rt0;
        K = K0;
        w = w0;
        puts("meow");
        PhotonMapping(0, i); 
        RayTracing(image, camera, i);
        
        for (int x = 0; x < camera->getWidth(); ++x) {
            for (int y = 0; y < camera->getHeight(); ++y) {
                Vector3f finalColor = c[x][y] / i + d[x][y];
                finalColor = Vector3f(clamp(finalColor[0]), clamp(finalColor[1]), clamp(finalColor[2]));
                finalColor = Vector3f(gamma(finalColor[0]), gamma(finalColor[1]), gamma(finalColor[2]));
                finalColor *= 1.5;
                image->SetPixel(x, y, finalColor);
            }
        }
        image->SaveBMP(outputFile);

        puts("meow");
        R = Rt0 / 2;
        K = K0 / 2;
        w = wt;
        PhotonMapping(2, i);
        if (pm.size() != 0) {
            RayTracing(image, camera, i);
        }
        
        for (int x = 0; x < camera->getWidth(); ++x) {
            for (int y = 0; y < camera->getHeight(); ++y) {
                c[x][y] = c[x][y] / i;
                Vector3f finalColor = c[x][y];
                finalColor = Vector3f(clamp(finalColor[0]), clamp(finalColor[1]), clamp(finalColor[2]));
                finalColor = Vector3f(gamma(finalColor[0]), gamma(finalColor[1]), gamma(finalColor[2]));
                finalColor *= 1.5;
                image->SetPixel(x, y, finalColor);
            }
        }
        image->SaveBMP(outputFile);
    }
}

int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }

    if (argc != 3) {
        cout << "Usage: ./bin/PA1 <input scene file> <output bmp file>" << endl;
        return 1;
    }
    char *inputFile = argv[1];
    char *outputFile = argv[2];  // only bmp is allowed.

    // TODO: Main RayCasting Logic
    // First, parse the scene using SceneParser.
    // Then loop over each pixel in the image, shooting a ray
    // through that pixel and finding its intersection with
    // the scene.  Write the color at the intersection to that
    // pixel in your output image.
    cout << "Hello! Computer Graphics!" << endl;

    SceneParser scene(inputFile);
    sceneP = &scene;
    Camera *camera = scene.getCamera();
    int width = camera->getWidth(), height = camera->getHeight();
    Image image(width, height);
    Group *group = scene.getGroup();

    printf("mmm GroupSize: %d\n", group->getGroupSize());

    // PathTracing(&image, camera);
     PPM(&image, camera, outputFile);

    /* for (int x = 0; x < camera->getWidth(); ++x) {
        for (int y = 0; y < camera->getHeight(); ++y) {
            Vector3f finalColor = c[x][y];
            image.SetPixel(x, y, c[x][y]);
        }
    } */
    image.SaveBMP(outputFile);

    return 0;
}

