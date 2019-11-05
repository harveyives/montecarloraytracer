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

bool object_occluded(std::vector<Object *> &objects, Vertex &hit_position, Vertex &light_position);

Vector raytrace(Scene *scene, Ray ray, int depth);

int main(int argc, char *argv[])
{
  int width = 256;
  int height = 256;

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(0,0,1);
  Vector up = Vector(0,1,0);
  // TODO move this to the scene?
  Camera *camera = new Camera(eye, look, up, 1, 90, height, width);
  Scene *scene = new Scene(0.2);

  for(int c = 0; c < width; c++) {
      for(int r = 0; r < height; r++) {
          Ray ray = Ray(eye, camera->get_ray_direction(c, r));

          Vector colour = raytrace(scene, ray, 2);
          fb->plotPixel(c, r, colour.x, colour.y, colour.z);
      }
  }

  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  return 0;
}

Vector raytrace(Scene *scene, Ray ray, int depth = 1) {
    Hit hit = Hit();
    Vector colour = {0,0,0};
    if(depth == 0) return colour;
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

            //thus self occlusion
            if (diffuse < 0) {
                continue;
            }

            // specular
            Vector reflection = Vector();
            hit.normal.reflection(light_direction, reflection);
            float specular = reflection.dot(ray.direction);
            // thus no contribution
            if (specular < 0) {
                specular = 0.0;
            }
            //TODO colour individual channels? 
            colour = colour + hit_colour * diffuse * hit.what->material.kd + hit_colour * pow(specular, 128) * hit.what->material.ks;
        }
        colour = colour + hit_colour * scene->ka;
        colour = colour / scene->lights.size();


        // tODO add conditional here to ignore reflection && refraction if coeffs really small
        Ray reflection_ray;
        hit.normal.reflection(ray.direction, reflection_ray.direction);
        reflection_ray.position = hit.position;
        // TODO find a cleaner way of doing this:
        reflection_ray.position = reflection_ray.get_point(0.001);
        colour = colour + hit.what->material.kr * raytrace(scene, reflection_ray, --depth);
//
//        Ray transparency_ray;
//        hit.normal.reflection(ray.direction, transparency_ray.direction);
//        transparency_ray.position = hit.position;
//        // TODO find a cleaner way of doing this:
//        transparency_ray.position = reflection_ray.get_point(0.001);
//        colour = colour + hit.what->material.kt * raytrace(scene, transparency_ray, --depth);

    }
    return colour;
}

bool object_occluded(std::vector<Object *> &objects, Vertex &hit_position, Vertex &light_position) {
    Ray shadow = Ray();
    shadow.position = hit_position;
    shadow.direction = light_position - hit_position;
    // don't cast shadows from objects that are further away than the light
    float distance_to_light = shadow.direction.magnitude();
    shadow.direction.normalise();

    // shifting so that the face does not self-shadow
    shadow.position = shadow.get_point(0.001);

    Hit obj_shadow_hit = Hit();
    for (Object *obj : objects) {
        obj_shadow_hit.flag = false;

        obj->intersection(shadow, obj_shadow_hit);
        if(obj_shadow_hit.flag && obj_shadow_hit.t < distance_to_light) {
            return true;
        }
    }
    return false;
}
