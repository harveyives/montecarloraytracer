//
// Created by harve on 23/10/2019.
//

#ifndef CODE_DIRECTIONAL_LIGHT_H
#define CODE_DIRECTIONAL_LIGHT_H

#include "light.h"

class DirectionalLight : public Light {
public:
    DirectionalLight(Vector d) {
        position = {0,0,0};
        direction = d;
    }

    Vector get_light_direction(Vertex p) override {
        return direction;
    }
};

#endif //CODE_DIRECTIONAL_LIGHT_H
