//
// Created by Harvey Ives on 10/10/2019.
//
#include "camera.h"
#include "vector.h"
#include "vertex.h"

Camera::Camera(Vertex eye, Vertex look, Vector up, double d, float fov, int height, int width) {
    // TODO clear this stuff up
    this->eye = eye;
    this->look = look;
    this->up = up;
    this->d = d;

    w = Vector(look.x, look.y, look.z, eye.x, eye.y, eye.z);
    w.normalise();

    up.cross(w, u);
    u.normalise();

    w.cross(u, v);

    fov_radians = M_PI * (fov / 2) / 180;
    aspect_ratio = height / width;

    half_width = tan(fov_radians);
    half_height = aspect_ratio * half_width;

    pixel_width = half_width * 2 / (width - 1.0);
    pixel_height = half_height * 2 / (height - 1.0);
}

Vector Camera::get_ray_direction(int x, int y) {
    float xv = (x * pixel_width) - half_width;
    float yv = (y * pixel_height) - half_height;

    Vector D = u * xv + v * yv - w * d;
    D.normalise();

    return D;
}
