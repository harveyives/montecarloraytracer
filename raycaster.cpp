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

bool
object_occluded(std::vector<Object *> &objects, Vertex &hit_position, Vertex &light_position);

int main(int argc, char *argv[])
{
  int width = 100;
  int height = 100;

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(0,0,1);
  Vector up = Vector(0,1,0);

  Camera *camera = new Camera(eye, look, up, 1, 90, height, width);
  Scene *scene = new Scene(0.2);

  for(int c = 0; c < width; c++) {
      for(int r = 0; r < height; r++) {
          Ray ray = Ray(eye, camera->get_ray_direction(c, r));
          Hit hit = Hit();

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
              Vector colour = {0, 0, 0};
              Vector hit_colour = hit.what->material.colour;

              for(Light *light : scene->lights) {
                  if(object_occluded(scene->objects, hit.position, light->position)) {
                      continue;
                  }

                  //get light direction based on point if point light, otherwise get directional light
                  Vector light_direction = light->get_light_direction(hit.position);
                  light_direction.normalise();

                  // diffuse
                  float diffuse =  light_direction.dot(hit.normal);

                  //thus self occulusion
                  if (diffuse < 0)
                      continue;

                  // specular
                  Vector reflection = Vector();
                  hit.normal.reflection(light_direction, reflection);

                  float specular = reflection.dot(ray.direction);
                  if (specular < 0)
                      specular = 0.0;

                  colour = colour + hit_colour * diffuse * hit.what->material.kd + hit_colour * pow(specular, 128) * hit.what->material.ks;
              }
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

bool object_occluded(std::vector<Object *> &objects, Vertex &hit_position, Vertex &light_position) {
    Hit obj_shadow_hit = Hit();
    for (Object *obj : objects) {
        obj_shadow_hit.flag = false;

        Ray shadow_ray = Ray();
        shadow_ray.position = hit_position;
        shadow_ray.direction = light_position - hit_position;
        shadow_ray.direction.normalise();
        shadow_ray.position = shadow_ray.get_point(0.001);

        obj->intersection(shadow_ray, obj_shadow_hit);
        if(obj_shadow_hit.flag) {
            return true;
        }
    }
    return false;
}
