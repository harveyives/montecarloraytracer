#ifndef CODE_POINT_LIGHT_H
#define CODE_POINT_LIGHT_H

#include "light.h"

class PointLight : public Light {
public:
    PointLight(Vertex p) : PointLight(p, Vector(0, -1, 0), Vector(255, 255, 255)) {}

    PointLight(Vertex p, Vector d, Vector i) {
        position = p;
        direction = d;
        intensity = i;
    }

    Vector get_direction(Vertex p) override {
        direction = position - p;
        direction.normalise();
        return direction;
    }
};

#endif //CODE_POINT_LIGHT_H
