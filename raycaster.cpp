/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */
#include "framebuffer.h"
#include "polymesh.h"
#include "sphere.h"
#include "camera.h"
#include "vertex.h"
#include "vector.h"

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

  float fov = 50;

  Camera *camera = new Camera(eye, look, up, 1, fov, height, width);

  Sphere *sphere = new Sphere(Vertex(200,0,0), 160);
  for(int c = 0; c < width; c++)
  {
      for(int r = 0; r < height; r++)
      {
          Ray ray = Ray(eye, camera->get_ray_direction(c,r));
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
