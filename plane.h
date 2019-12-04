//
// Created on 04/11/2019.
//
#ifndef CODE_PLANE_H
#define CODE_PLANE_H

#include "vertex.h"
#include "vector.h"
#include "object.h"
#include "phong.h"

class Plane : public Object {
    Vector normal;
    Vertex position;
public:
    Plane(Vertex p, Vector n) : Plane(p, n, new Phong()) {}

    Plane(Vertex p, Vector n, Material *material) : Object(material) {
        position = p;
        normal = n;
    }
    void intersection(Ray ray, Hit &hit);
};
#endif //CODE_PLANE_H