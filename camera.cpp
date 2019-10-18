//
// Created by Harvey Ives on 10/10/2019.
//
#include "camera.h"
#include "vector.h"
#include "vertex.h"

Camera::Camera(Vertex eye, Vertex look, Vector up, double d, float fov, int height, int width) {
    Camera::eye = eye;
    Camera::look = look;
    Camera::up = up;
    Camera::d = d;

    Camera::w = Vector(look.x, look.y, look.z, eye.x, eye.y, eye.z);
    w.normalise();

    Camera::up.cross(Camera::w, Camera::u);
    u.normalise();

    Camera::w.cross(Camera::u, Camera::v);

    Camera::fov_radians = M_PI * (fov / 2) / 180;
    Camera::aspect_ratio = height / width;

    Camera::half_width = tan(fov_radians);
    Camera::half_height = aspect_ratio * half_width;

    Camera::pixel_width = half_width * 2 / (width - 1.0);
    Camera::pixel_height = half_height * 2 / (height - 1.0);
}

Vector Camera::get_ray_direction(int x, int y) {
    float xv = (x * Camera::pixel_width) - Camera::half_width;
    float yv = (y * Camera::pixel_height) - Camera::half_height;


    Vector D = u * xv + v * yv - w * d;
    D.normalise();

    return D;
}
