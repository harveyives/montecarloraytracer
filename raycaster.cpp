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

int main(int argc, char *argv[])
{
  int width = 100;
  int height = 100;

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(0,0,1);
  Vector up = Vector(0,-1,0);

  Camera *camera = new Camera(eye, look, up, 1, 65, height, width);
  Scene *scene = new Scene(0.2);

  for(int c = 0; c < width; c++) {
      for(int r = 0; r < height; r++) {
          Ray ray = Ray(eye, camera->get_ray_direction(c, r));
          Hit hit = Hit();
          hit.t = MAXFLOAT;

          for (Object *obj : scene->objects) {
              Hit obj_hit = Hit();

              obj_hit.t = MAXFLOAT;
              obj->intersection(ray, obj_hit);
              if(obj_hit.flag) {
                  if(obj_hit.t < hit.t) {
                      hit = obj_hit;
                  }
              }
          }
          if (hit.flag) {
              Vector colour = {0, 0, 0};
              int n = 128;
              Vector hit_colour = hit.what->material.colour;
              for(Light *light : scene->lights) {

                  Hit obj_shadow_hit = Hit();
                  obj_shadow_hit.t = MAXFLOAT;
                  for (Object *obj : scene->objects) {
                      obj_shadow_hit.flag = false;

                      Ray shadow_ray = Ray();
                      shadow_ray.position = hit.position;
                      shadow_ray.direction = light->position - hit.position;
                      shadow_ray.direction.normalise();
                      shadow_ray.position = shadow_ray.get_point(0.001);

                      obj->intersection(shadow_ray, obj_shadow_hit);
                      if(obj_shadow_hit.flag) {
                          break;
                      }
                  }
                  if(obj_shadow_hit.flag) {
                      continue;
                  }

                  //get light direction based on point if point light, otherwise get directional light
                  Vector light_direction = light->get_light_direction(hit.position);
                  light_direction.normalise();

                  float diffuse =  light_direction.dot(hit.normal);

                  //thus self occulusion
                  if (diffuse < 0)
                      continue;

                  // calculate specular component
                  Vector reflection = Vector();
                  hit.normal.reflection(light_direction, reflection);

                  float specular = reflection.dot(ray.direction);
                  if (specular < 0.0)
                      specular = 0.0;

                  colour = colour + hit_colour * diffuse * hit.what->material.diffuse + hit_colour * pow(specular, n) * hit.what->material.specular;
              }
              //TODO divide by number of lights?
              colour = colour + hit_colour * scene->ka;
              colour = colour / scene->lights.size();
              fb->plotPixel(c, r, colour.x, colour.y, colour.z);
          }
      }
  }

  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  return 0;
}
