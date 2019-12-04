//
// Created by harve on 04/11/2019.
//

#ifndef CODE_PHONG_H
#define CODE_PHONG_H

#include "vector.h"
#include "colour.h"
#include "material.h"

class Phong : public Material {
public:
    Phong(Vector c = {255, 255, 255}, float specular = 0.2, float diffuse = 0.5, float refractive_index = 1,
          float reflective = 0, float transparent = 0) {
        colour = c;
        ks = specular;
        kd = diffuse;
        ior = refractive_index;
        r = reflective;
        t = transparent;
    }

    Vector compute_base(Hit &hit) override {
        return colour;
    }
};

#endif //CODE_PHONG_H
