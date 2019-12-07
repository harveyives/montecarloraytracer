#include "camera.h"
#include "vector.h"
#include "vertex.h"

// class to model the camera, calculates rays based on location, size of image, etc...

Camera::Camera(Vertex eye, Vertex look, Vector up, double d, float fov, int height, int width) {
    this->eye = eye;
    this->look = look;
    this->up = up;
    this->d = d;

    w = eye - look;
    w.normalise();

    w.cross(up, u);
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
