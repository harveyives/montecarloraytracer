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

    Light()
    {
    }

    Light(Vertex p, Vector d)
    {
        position = p;
        direction = d;
        direction.normalise();
    }
};

#endif //CODE_LIGHT_H
