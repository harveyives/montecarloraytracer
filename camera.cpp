//
// Created by Harvey Ives on 10/10/2019.
//
#include "camera.h"
#include "vector.h"
#include "vertex.h"

Camera::Camera(Vertex eye, Vertex look, Vector up, double d){
    Camera::eye = eye;
    Camera::look = look;
    Camera::up = up;
    Camera::d = d;

    Camera::w = Vector(look.x, look.y, look.z, eye.x, eye.y, eye.z);
    w.normalise();

    Camera::up.cross(Camera::w, Camera::u);
    u.normalise();

    Camera::w.cross(Camera::u, Camera::v);
}

Vector Camera::get_ray_direction(float xv, float yv) {
    Vector D = u * xv + v * yv - w * d;
    D.normalise();

    return D;
}
