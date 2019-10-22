/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include <vector>
#include "vertex.h"
#include "object.h"

class Sphere : public Object {
    Vertex center;

    float  radius;

public:
    Sphere(Vertex c, float r) : Sphere(c, r, Colour(255,255,255)) {}
    Sphere(Vertex c, float r, Colour colour) : Object(colour) {
        center = c;
        radius = r;
    }
    void intersection(Ray ray, Hit &hit);
};
