/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */
#include <iostream>
#include "framebuffer.h"
#include "linedrawer.h"
#include "polymesh.h"
#include "sphere.h"
#include "camera.h"
#include "vertex.h"
#include "vector.h"
#include <float.h>

using namespace std;

int main(int argc, char *argv[])
{
    int width = 100;
    int height = 100;
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f,0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 7.0f,0.0f,0.0f,0.0f,1.0f);

  // Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(0,0,1);
  Vector up = Vector(0,-1,0);

  float d = 1;

  float fov = 90;
  //TODO move this to the camera
  float fovRadians = M_PI * (fov / 2) / 180;
  float aspect_ratio = height / width;

  float half_width = tan(fovRadians);
  float half_height = aspect_ratio * half_width;

  float pixel_width = half_width * 2 / (width - 1.0);
  float pixel_height = half_height * 2 / (height - 1.0);

  Camera *camera = new Camera(eye, look, up, d);

  Sphere *sphere = new Sphere(Vertex(200,0,0), 160);
  for(int c = 0; c < width; c++)
  {
      for(int r = 0; r < height; r++)
      {
          float xv = (c * pixel_width) - half_width;
          float yv = (r * pixel_height) - half_height;

          Ray ray = Ray(eye, camera->get_ray_direction(xv,yv));
          Hit hit = Hit();

          pm->intersection(ray, hit);
          if(hit.flag) {
              fb->plotDepth(c, r, hit.position.magnitude());
          }
      }
  }


  // Output the framebuffer.
//  fb->writeRGBFile((char *)"test.ppm");
  fb->writeDepthFile((char *)"test.ppm");

  return 0;
  
}
