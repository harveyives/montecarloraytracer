#ifndef CODE_PHONG_H
#define CODE_PHONG_H

#include "vector.h"
#include "material.h"

// Class to model simple phong textures, implementations are in material superclass

class Phong : public Material {
public:
    Phong(Vector c = {255, 255, 255}, float specular = 0.2, float diffuse = 0.5, float refractive_index = 1,
          float reflective = 0, float transparent = 0, Vector emissive = Vector(0, 0, 0)) : Material() {
        colour = c;
        ks = specular;
        kd = diffuse;
        ior = refractive_index;
        r = reflective;
        t = transparent;
        e = emissive;
    }

    Vector compute_base_colour(Hit &hit) override {
        return colour;
    }
};

#endif //CODE_PHONG_H
