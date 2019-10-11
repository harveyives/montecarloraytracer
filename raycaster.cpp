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
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(2048,2048);

  // Read in the teapot model.
//  PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(1,0,0);
  Vector up = Vector(0,1,0);
  float d = 6;

  //change this to fov ?

  float s = 0.1;//tan(deg2rad(fov * 0.5));

  Camera *camera = new Camera(eye, look, up, d);

  std::cout << "w: ("<< camera->w.x << ", " << camera->w.y << ", " << camera->w.z << ")" << std::endl;
  std::cout << "u: ("<< camera->u.x << ", " << camera->u.y << ", " << camera->u.z << ")" << std::endl;
  std::cout << "v: ("<< camera->v.x << ", " << camera->v.y << ", " << camera->v.z << ")" << std::endl;
  int width = 2048;
  int height = 2048;

  Vertex sphereC = Vertex(10,0,0);
  Sphere *sphere = new Sphere(sphereC, 9);
  for(int c = 0; c < width; c++)
  {
      for(int r = 0; r < height; r++)
      {

          float xv = s * (c - width * 0.5);
          float yv = s * (r - height * 0.5);
          Vector D = camera->u * xv + camera->v * yv - camera->w * d;
          D.normalise();
          Ray ray = Ray(eye, D);
          Hit hit = Hit();
          sphere->intersection(ray, hit);
          //hit.position.magnitude();
          if(hit.flag == true) {
              fb->plotPixel(c, r, 1.0f, 1.0f, 1.0f);
          }
      }
  }


  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");

  return 0;
  
}
