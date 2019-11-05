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
    float kr;
    float kt;

    Material(Vector c = {255, 255, 255}, float specular = 0.2, float diffuse = 0.5, float reflection = 0, float transparency = 0) {
        colour = c;
        ks = specular;
        kd = diffuse;
        kr = reflection;
        kt = transparency;
    }


    static Material set_colour(Vector c) {
        return {c};
    }
};
#endif //CODE_MATERIAL_H
