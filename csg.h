//
// Created on 04/11/2019.
//
#ifndef CODE_CSG_H
#define CODE_CSG_H

#include "object.h"

using namespace std;

class CSG : public Object {
public:
    Object *u;
    Object *v;
    string operation;

    CSG(Object *u, Object *v, string o) : CSG(u, v, o, Material()) {}

    CSG(Object *u, Object *v, string o, Material m) : Object(m) {
        this->u = u;
        this->v = v;
        operation = o;
    }

    void intersection(Ray ray, Hit &hit) override;
};

#endif //CODE_CSG_H