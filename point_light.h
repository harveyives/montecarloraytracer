//
// Created by harve on 23/10/2019.
//

#ifndef CODE_POINT_LIGHT_H
#define CODE_POINT_LIGHT_H

#include "light.h"

class PointLight : public Light {
public:
    PointLight(Vertex p) {
        position = p;
    }

    Vector get_light_direction(Vertex p) override {
        direction = position - p;
        direction.normalise();
        return direction;
    }
};

#endif //CODE_POINT_LIGHT_H
