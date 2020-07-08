#ifndef CODE_TEXTURE_H
#define CODE_TEXTURE_H

#include "vector.h"
#include "material.h"
#include "perlin.h"

#include "hit.h"

class Marble : public Material {
public:
    Perlin *perlin;

    Marble(Vector c = {255, 255, 255}, float specular = 0.2, float diffuse = 0.5, float refractive_index = 1,
           float reflective = 0, float transparent = 0, Vector emissive = Vector(0, 0, 0)) :
            Material(Vector(), 0, 0, 0, 0, 0,
                     Vector(0, 0, 0)) {
        colour = c;
        ks = specular;
        kd = diffuse;
        ior = refractive_index;
        r = reflective;
        t = transparent;
        e = emissive;

        // parameters that best represent marble texture
        perlin = new Perlin(2, 256, 5, 10);
    }

    Vector base_colour(Hit &hit) override {
        Object *obj = hit.what;
        long double min_x = abs(obj->min_limit_x);
        long double min_y = abs(obj->min_limit_y);
        long double max_x = abs(obj->max_limit_x);
        long double max_y = abs(obj->max_limit_y);

        // find ratio of hit position with respect to entire object
        long double scaled_x = fabs((hit.position.x + min_x) / (min_x + max_x));
        long double scaled_y = fabs((hit.position.y + min_y) / (min_y + max_y));

        // use ratio to find coordinates within perlin noise
        int dimension = perlin->n;
        long double u = scaled_x * dimension;
        long double v = scaled_y * dimension;

        // use coordinates to find the noise value
        Vector perlin_colour = colour * perlin->matrix[(int) u][(int) v];
        
        // amount of perlin to contribute
        float kp = 0.6;
        return colour * (1 - kp) + perlin_colour * kp;
    }
};

#endif //CODE_TEXTURE_H
