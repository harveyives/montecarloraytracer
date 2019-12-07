#ifndef CODE_MATERIAL_H
#define CODE_MATERIAL_H

#include "vector.h"

// abstract class to model materials for objects (phong or textured).

class Material {
protected:
    Vector colour;
public:
    float ks;
    float kd;
    float ior;
    float r;
    float t;
    Vector e;

    Material(Vector c = {255, 255, 255}, float specular = 0.2, float diffuse = 0.5, float refractive_index = 1,
             float reflective = 0, float transparent = 0, Vector emissive = Vector(0, 0, 0)) {
        colour = c;
        ks = specular;
        kd = diffuse;
        ior = refractive_index;
        r = reflective;
        t = transparent;
        e = emissive;
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

        result = result + specular(view, light_direction, normal, base_colour);
        return result;
    }

    virtual Vector specular(Vector &view, Vector &light_direction, Vector &normal, Vector &base_colour) {
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

    virtual bool is_emissive() {
        return (e.x != 0 && e.y != 0 && e.z != 0);
    }

    virtual bool is_transmissive() {
        return (t != 0);
    }

    virtual bool is_refractive() {
        return (r != 0);
    }
};

#endif //CODE_MATERIAL_H
