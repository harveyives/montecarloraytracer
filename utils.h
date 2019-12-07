#ifndef CODE_UTILS_H
#define CODE_UTILS_H

#include <math.h>
#include <string>
#include <iostream>
#include <cstring>
#include <random>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <float.h>
#include <vector>
#include "vector.h"
#include "vertex.h"
#include "camera.h"
#include "sphere.h"
#include "scene.h"
#include "alglib/stdafx.h"
#include "alglib/alglibmisc.h"

using namespace std;

class Utils {
public:
    // split string based on tabs or spaces
    static vector<string> split_string(string line) {
        istringstream ss(line);
        istream_iterator<string> begin(ss), end;
        vector<string> words(begin, end);
        return words;
    }

    // generate random float between two values
    static float get_random_number(float min, float max) {
        random_device device;
        mt19937 random(device());
        uniform_real_distribution<> distribution(min, max);
        auto n = distribution(random);
        return n;
    }

    // generate random direction in direction within cone (defined by size of theta)
    static Vector random_direction(Vector direction, float theta) {
        // make basis vectors based on direction
        Vector basis_z = direction;
        basis_z.normalise();

        Vector basis_x;
        Vector different = (0.5 > fabs(basis_z.x)) ? Vector(1, 0, 0) : Vector(0, 1, 0);
        basis_z.cross(different, basis_x);
        basis_x.normalise();

        Vector basis_y;
        basis_x.cross(basis_z, basis_y);

        // height
        float z = Utils::get_random_number(cos(theta), 1);

        // rotation
        float phi = Utils::get_random_number(-M_PI, +M_PI);
        float x = sqrt(1 - z * z) * cos(phi);
        float y = sqrt(1 - z * z) * sin(phi);

        // combine
        Vector result = x * basis_x + y * basis_y + z * basis_z;
        result.normalise();
        return result;
    }

    // linear interpolation
    static float lerp(double t, double a, double b) {
        return a + t * (b - a);
    }

    static float clamp(float n, float lower, float upper) {
        return max(lower, min(n, upper));
    }
};


#endif //CODE_UTILS_H
