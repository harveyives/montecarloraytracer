#ifndef CODE_LIGHT_H
#define CODE_LIGHT_H

#include "vertex.h"

// abstract class for lights

class Light {
public:
    Vertex position;
    Vector direction;
    Vector intensity;

    virtual Vertex get_position() {
        return position;
    }

    virtual Vector get_direction(Vertex point) = 0;
};

#endif //CODE_LIGHT_H
