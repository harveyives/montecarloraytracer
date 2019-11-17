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
#include "photon_map.h"

using namespace std;

Hit raytrace(Ray &ray, Hit &hit);

PhotonMap *pm = new PhotonMap();
int main(int argc, char *argv[])
{
  int width = 150;
  int height = 150;

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(0,0,1);
  Vector up = Vector(0,1,0);
  // TODO move this to the scene?
  Camera *camera = new Camera(eye, look, up, 1, 50, height, width);
  Scene *scene = new Scene(0.9);



  for(int c = 0; c < width; c++) {
      for(int r = 0; r < height; r++) {
          Ray ray = Ray(eye, camera->get_ray_direction(c, r));

          Vector colour = scene->compute_colour(ray, 3);
          fb->plotPixel(c, r, colour.x, colour.y, colour.z);
      }
  }

//   Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  return 0;
}

