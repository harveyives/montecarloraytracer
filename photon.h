#ifndef CODE_PHOTON_H
#define CODE_PHOTON_H

#include "vector.h"
#include "vertex.h"
#include "ray.h"
#include <string>


using namespace std;

// Class to model photons that are used for photon mapping

class Photon {
public:
    Ray ray;
    Vector colour;
    // strings are used because they are easier serialise for saving
    string type;

    Photon(Ray r, Vector c, string t) {
        ray = r;
        colour = c;
        type = t;
    }
};

#endif //CODE_PHOTON_H