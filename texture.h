//
// Created by harve on 04/11/2019.
//

#ifndef CODE_TEXTURE_H
#define CODE_TEXTURE_H

#include "vector.h"
#include "colour.h"
#include "material.h"
#include "perlin.h"

#include "hit.h"

class Texture : public Material {
public:
    Perlin *perlin = new Perlin();
    Vertex centre;

    Texture(Vector c = {255, 255, 255}, float specular = 0.2, float diffuse = 0.5, float refractive_index = 1,
            float reflective = 0, float transparent = 0) {
        colour = c;
        ks = specular;
        kd = diffuse;
        ior = refractive_index;
        r = reflective;
        t = transparent;
    }

    Vector compute_base_colour(Hit &hit) override {
        Object *obj = hit.what;
        long double min_x = abs(obj->min_limit_x);
        long double min_y = abs(obj->min_limit_y);
        long double max_x = abs(obj->max_limit_x);
        long double max_y = abs(obj->max_limit_y);

        // find ratio of hit position with regards to pixel space
        long double ratio_x = fabs((hit.position.x + min_x) / (min_x + max_x));
        long double ratio_y = fabs((hit.position.y + min_y) / (min_y + max_y));

        // find coordinates within perlin noise corresponding to the ratio
        int i = (int) (ratio_x * perlin->n);
        int j = (int) (ratio_y * perlin->n);
        Vector perlin_colour = colour * perlin->matrix[i][j];
        // amount of perlin to contribute
        float kp = 0.6;
        return colour * (1 - kp) + perlin_colour * kp;
    }
};

#endif //CODE_TEXTURE_H
