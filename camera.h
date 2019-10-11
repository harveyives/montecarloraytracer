//
// Created by Harvey Ives on 10/10/2019.
//
#pragma once

#include "vertex.h"
#include "vector.h"

class Camera {
public:
    Vertex eye;
    Vertex look;
    Vector up;
    double d;

    Vector w;
    Vector u;
    Vector v;

    Camera(Vertex eye, Vertex look, Vector up, double d);
};
