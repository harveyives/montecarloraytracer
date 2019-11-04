//
// Created by harve on 04/11/2019.
//

#ifndef CODE_MATERIAL_H
#define CODE_MATERIAL_H

#include "vector.h"

class Material {
public:
    Vector colour;
    float specular;
    float diffuse;

    Material(Vector c = {255, 255, 255}, float s = 0.3, float d = 0.5) {
        colour = c;
        specular = s;
        diffuse = d;
    }
};
#endif //CODE_MATERIAL_H
