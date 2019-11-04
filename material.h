//
// Created by harve on 04/11/2019.
//

#ifndef CODE_MATERIAL_H
#define CODE_MATERIAL_H

#include "vector.h"

class Material {
public:
    Vector colour;
    float ks;
    float kd;

    Material(Vector c = {255, 255, 255}, float specular = 0.3, float diffuse = 0.5) {
        colour = c;
        ks = specular;
        kd = diffuse;
    }
};
#endif //CODE_MATERIAL_H
