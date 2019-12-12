#include "perlin.h"
#include "utils.h"

using namespace std;

// Class used to model marble, wood, clouds textures that can be sampled for texturing objects
// Perlin Noise by Ken Perlin
// https://mrl.nyu.edu/~perlin/noise/
// with advice from https://lodev.org/cgtutor/randomnoise.html
Perlin::Perlin(float turbulence_power, float turbulence_size, float x_period, float y_period) {
    // convert permutations to be within colour range 0-255
    for (int i = 0; i < permutations_original.size(); i++)
        permutations[i] = permutations_original[i % 256];

    // building noise array
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = sample(i, j, turbulence_power, turbulence_size, x_period / n, y_period / n);
        }
    }
}

Vector Perlin::sample(int i, int j, float turbulence_power, float turbulence_size, float x_period, float y_period) {
    // generate output value
    float val = fabs(sin(x_period * i / n + y_period * j / n + turbulence_power * turbulence(j, i, 0.8, turbulence_size) * M_PI));
    return Vector(val, val, val);
}

float Perlin::noise(float x, float y, float z) {
    float f_x, f_y, f_z;
    int i_x, i_y, i_z;

    // find unit cube that contains point xyz
    i_x = (int) floor(x) & 255;
    i_y = (int) floor(y) & 255;
    i_z = (int) floor(z) & 255;
    // find relative xyz in cube
    x -= floor(x);
    y -= floor(y);
    z -= floor(z);
    // compute fade curves for xyz
    f_x = fade(x);
    f_y = fade(y);
    f_z = fade(z);
    //hash coordinates of 8 cube corners
    vector<int> corners = {
            permutations[permutations[i_x] + i_y] + i_z,
            permutations[permutations[i_x + 1] + i_y] + i_z,
            permutations[permutations[i_x] + i_y + 1] + i_z,
            permutations[permutations[i_x + 1] + i_y + 1] + i_z
    };

    return Utils::lerp(f_z, Utils::lerp(f_y, Utils::lerp(f_x, gradient(permutations[corners[0]], x, y, z),
                                                         gradient(permutations[corners[1]], x - 1, y, z)),
                                        Utils::lerp(f_x, gradient(permutations[corners[2]], x, y - 1, z),
                                                    gradient(permutations[corners[3]], x - 1, y - 1, z))),
                       Utils::lerp(f_y, Utils::lerp(f_x, gradient(permutations[corners[0] + 1], x, y, z - 1),
                                                    gradient(permutations[corners[1] + 1], x - 1, y, z - 1)),
                                   Utils::lerp(f_x, gradient(permutations[corners[2] + 1], x, y - 1, z - 1),
                                               gradient(permutations[corners[3] + 1], x - 1, y - 1, z - 1))));
}

float Perlin::gradient(int hash, double x, double y, double z) {
    // convert lo 4 bits of hash into 12 gradient directions
    int h = hash & 15;
    double u = h < 8 ? x : y, v = h < 4 ? y : h == 12 || h == 14 ? x : z;
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

float Perlin::turbulence(float x, float y, float z, float scale) {
    float value = 0;
    float size = scale;

    // add samples for different scales
    while (scale >= 1) {
        value += fabs(noise(x / scale, y / scale, z) * scale);
        scale = scale / 2;
    }
    // normalise value
    return value / (2 * size);
}

float Perlin::fade(double t) {
    //smooth noise
    return pow(t, 3) * 10
           - pow(t, 4) * 15
           + pow(t, 5) * 6;
}
