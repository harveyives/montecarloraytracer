#pragma once

#include "vertex.h"
#include "vector.h"


// class to model the camera, calculates rays based on location, size of image, etc...

class Camera {
public:
    Vertex eye;
    Vertex look;
    Vector up;
    double d;

    Vector w;
    Vector u;
    Vector v;


    float fov_radians;
    float aspect_ratio;

    float half_width;
    float half_height;

    float pixel_width;
    float pixel_height;

    Camera(Vertex eye, Vertex look, Vector up, double d, float fov, int height, int width);

    Vector get_ray_direction(int xv, int yv);
};
