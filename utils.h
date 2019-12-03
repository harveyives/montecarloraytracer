//
// Created by harve on 02/12/2019.
//

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
    static vector<string> split_string(string line) {
        istringstream ss(line);
        istream_iterator<string> begin(ss), end;
        vector<string> words(begin, end);
        return words;
    }

    static float get_random_number(int min, int max) {
        random_device device;
        mt19937 random(device());
        uniform_real_distribution<> distribution(min, max);
        auto n = distribution(random);
        return n;
    }

    static Vector random_direction(Vector direction, float theta) {
        // Make an orthogonal basis whose third vector is along `direction'
        Vector b3 = direction;
        b3.normalise();
        Vector different = (fabs(b3.x) < 0.5f) ? Vector(1.0f, 0.0f, 0.0f) : Vector(0.0f, 1.0f, 0.0f);
        Vector b1;
        b3.cross(different, b1);
        b1.normalise();
        Vector b2;
        b1.cross(b3, b2);
        //TODO rework this
        // Pick (x,y,z) randomly around (0,0,1)
        float z = Utils::get_random_number(cos(theta), 1);
        float r = sqrt(1.0f - z * z);
        float phi = Utils::get_random_number(-M_PI, +M_PI);
        float x = r * cos(phi);
        float y = r * sin(phi);

        // Construct the vector that has coordinates (x,y,z) in the basis formed by b1, b2, b3
        return x * b1 + y * b2 + z * b3;
    }
};


#endif //CODE_UTILS_H
