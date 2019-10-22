/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */
#include <iostream>
#include "framebuffer.h"
#include "sphere.h"
#include "camera.h"
#include "vertex.h"
#include "vector.h"
#include "scene.h"
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
  int width = 100;
  int height = 100;

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(0,0,1);
  Vector up = Vector(0,-1,0);

  Camera *camera = new Camera(eye, look, up, 1, 55, height, width);
  Scene *scene = new Scene();

  for(int c = 0; c < width; c++) {
      for(int r = 0; r < height; r++) {
          Ray ray = Ray(eye, camera->get_ray_direction(c, r));
          Hit hit = Hit();
          hit.t = MAXFLOAT;

          for (Object *obj : scene->objects) {
              Hit obj_hit = Hit();
              obj->intersection(ray, obj_hit);
              if(obj_hit.flag) {
                  if(obj_hit.t < hit.t) {
                      hit = obj_hit;
                  }
              }
          }
          if (hit.flag) {
              fb->plotDepth(c, r, hit.position.magnitude());
          }
      }
  }

  // Output the framebuffer.
//  fb->writeRGBFile((char *)"test.ppm");
  fb->writeDepthFile((char *)"test.ppm");

  return 0;
}
