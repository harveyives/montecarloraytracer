//
// Created by harve on 22/10/2019.
//

#ifndef CODE_LIGHT_H
#define CODE_LIGHT_H

#include "vertex.h"

class Light {
public:
    Vertex position;
    Vector direction;
    Vector intensity;

    virtual Vector get_light_direction(Vertex point) = 0;
};

#endif //CODE_LIGHT_H
