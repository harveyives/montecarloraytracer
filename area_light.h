#ifndef CODE_AREA_LIGHT_H
#define CODE_AREA_LIGHT_H

#include "light.h"
#include "utils.h"

// light where each ray is spawned at a different location within given limits

class AreaLight : public Light {
public:
    float spread;
    float min_x;
    float min_z;
    float max_x;
    float max_z;

    AreaLight(Vertex p, Vector d, Vector i, float s, float x1, float z1, float x2, float z2) {
        position = p;
        direction = d;
        intensity = i;
        spread = s;
        min_x = x1;
        min_z = z1;
        max_x = x2;
        max_z = z2;
    }

    Vertex get_position() override {
        // generate random point on surface
        Vertex point = {Utils::get_random_number(min_x, max_x), position.y, Utils::get_random_number(min_z, max_z)};

        // shift so that it doesn't self occulude
        Vector shift_bias = 0.001 * direction;
        return point + shift_bias;
    }

    Vector get_direction(Vertex p) override {
        return Utils::random_direction(direction, spread);
    }
};

#endif //CODE_AREA_LIGHT_H
