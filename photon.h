#ifndef CODE_PHOTON_H
#define CODE_PHOTON_H

#include "vector.h"
#include "vertex.h"
#include "ray.h"
#include <string>

using namespace std;

enum class photon_type {
    direct,
    indirect,
    shadow,
    caustic
};

class Photon {
public:
    Ray ray;
    Vector colour;
    photon_type type;

    Photon(Ray r, Vector c, photon_type t) {
        ray = r;
        colour = c;
        type = t;
    }
};

#endif //CODE_PHOTON_H