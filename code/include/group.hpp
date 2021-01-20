#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include <vector>


// TODO: Implement Group - add data structure to store a list of Object*
class Group : public Object3D {

public:

    Group() {

    }

    explicit Group (int num_objects) {

    }

    ~Group() override {

    }

    bool intersect(const Ray &r, Hit &h, float tmin, bool fl = 0) override {
        bool f = 0;
        for (auto it = list.begin(); it != list.end(); it ++) {
            Object3D* obj = (*it);
            f |= obj->intersect(r, h, tmin, fl);
        }
        return f;
    }

    void addObject(int index, Object3D *obj) {
        list.insert(list.begin() + index, obj);
    }

    int getGroupSize() {
        return list.size();
    }
    
    Vector3f getCenter() {}

    
    std::vector<Object3D*> list;

};

#endif
	
