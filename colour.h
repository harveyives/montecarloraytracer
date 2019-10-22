//
// Created by harve on 22/10/2019.
//

#ifndef CODE_COLOUR_H
#define CODE_COLOUR_H
class Colour {
public:
    int red;
    int green;
    int blue;

    Colour() : Colour(255,255,255) {}
    Colour(int r, int g, int b) {
        red = r;
        green = g;
        blue = b;
    }
};

#endif //CODE_COLOUR_H
