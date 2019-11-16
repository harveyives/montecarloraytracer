//
// Created by harve on 04/11/2019.
//

#ifndef CODE_MATERIAL_H
#define CODE_MATERIAL_H

#include "vector.h"
#include "colour.h"

class Material {
public:
    Vector colour;
    float ks;
    float kd;
    float ior;
    float r;
    float t;

    Material(Vector c = {255, 255, 255}, float specular = 0.2, float diffuse = 0.5, float refractive_index = 1,
             float reflective = 0, float transparent = 0) {
        colour = c;
        ks = specular;
        kd = diffuse;
        ior = refractive_index;
        r = reflective;
        t = transparent;
    }


    static Material set_colour(Vector c) {
        return {c};
    }
};
#endif //CODE_MATERIAL_H
