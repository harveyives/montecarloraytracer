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
};


#endif //CODE_UTILS_H
