//
// Created by harve on 04/11/2019.
//

#ifndef CODE_MATERIAL_H
#define CODE_MATERIAL_H

#include "vector.h"
#include "colour.h"

class Material {
protected:
    Vector colour;
public:
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

    virtual Vector compute_light_colour(Vector &view, Vector &light_direction, Vector &normal, Vector base_colour) {
        Vector result = Vector();

        // diffuse
        float diffuse = normal.dot(light_direction);
        //thus self occlusion
        if (diffuse < 0) {
            return result;
        }
        result = result + base_colour * diffuse * kd;

        result = result + get_specular_component(view, light_direction, normal, base_colour);
        return result;
    }

    virtual Vector get_specular_component(Vector &view, Vector &light_direction, Vector &normal, Vector &base_colour) {
        Vector result;
        Vector reflection = Vector();
        normal.reflection(light_direction, reflection);
        float specular = reflection.dot(view);
        // thus no contribution
        if (specular < 0) {
            specular = 0.0;
        }

        result = result + base_colour * pow(specular, 128) * ks;
        return result;
    };

    virtual Vector compute_base_colour(Hit &hit) = 0;
};
#endif //CODE_MATERIAL_H
