//
// Created on 05/11/2019.
//

#ifndef CODE_COLOUR_H
#define CODE_COLOUR_H
#pragma once

#include "vertex.h"
#include "vector.h"
class Colour {
public:
    int r;
    int g;
    int b;
    float a;

    Colour() : Colour(255,255,255){}
    Colour(int red, int green, int blue) : Colour(red, green, blue, 1){}
    Colour(int red, int green, int blue, float alpha) {
        r = red;
        g = green;
        b = blue;
        a = alpha;
    }
};
#endif //CODE_COLOUR_H
