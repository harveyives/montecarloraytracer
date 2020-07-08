#ifndef CODE_DIRECTIONAL_LIGHT_H
#define CODE_DIRECTIONAL_LIGHT_H

#include "light.h"

// light with direction only - position is not accounted for

class DirectionalLight : public Light {
public:
    DirectionalLight(Vector d) {
        position = {0,0,0};
        direction = d;
    }

    Vector get_direction(Vertex p) override {
        return direction;
    }
};

#endif //CODE_DIRECTIONAL_LIGHT_H
