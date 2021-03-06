/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */
#include <iostream>
#include <string>
#include <math.h>
#include "framebuffer.h"
#include "sphere.h"
#include "camera.h"
#include "vertex.h"
#include "vector.h"
#include "scene.h"

using namespace std;

// main method to raytrace an image

int main(int argc, char *argv[]) {
    int width = stoi(argv[1]);
    int height = stoi(argv[2]);
    bool generate_photon_map = stoi(argv[3]);

    // Create a framebuffer
    FrameBuffer *fb = new FrameBuffer(width, height);

    Vertex eye = Vertex(0, 0, 0);
    Vertex look = Vertex(0, 0, 1);
    Vector up = Vector(0, 1, 0);
    Camera *camera = new Camera(eye, look, up, 1, 50, height, width);
    Scene *scene = new Scene(0.6, generate_photon_map);

    // starting timer...
    clock_t start;
    double time;

    start = clock();
    cout << "Tracing... " << endl;
    for (int c = 0; c < width; c++) {
        for (int r = 0; r < height; r++) {
            Ray ray = Ray(eye, camera->get_ray_direction(c, r));

            Vector colour = scene->trace(ray, 4);
            fb->plotPixel(c, r, colour.x, colour.y, colour.z);
        }
    }
    cout << "DONE" << endl;
    time = (clock() - start) / ((double) CLOCKS_PER_SEC);

    std::cout << "printf: " << time << '\n';

    //   Output the framebuffer.
    cout << "Outputting..." << endl;
    fb->writeRGBFile((char *) "image.ppm");
    return 0;
}