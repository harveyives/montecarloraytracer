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
Vector raytrace(Scene *scene, Ray &ray, int depth);
float fresnel(float refractive_index, float cos_i);
Vector refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i);

int main(int argc, char *argv[])
{
  int width = 200;
  int height = 200;

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  Vertex eye = Vertex(0,0,0);
  Vertex look = Vertex(0,0,1);
  Vector up = Vector(0,1,0);
  // TODO move this to the scene?
  Camera *camera = new Camera(eye, look, up, 1, 65, height, width);
  Scene *scene = new Scene(0.9);

  for(int c = 0; c < width; c++) {
      for(int r = 0; r < height; r++) {
          Ray ray = Ray(eye, camera->get_ray_direction(c, r));

          Vector colour = raytrace(scene, ray, 3);
          fb->plotPixel(c, r, colour.x, colour.y, colour.z);
      }
  }

  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  return 0;
}

Vector raytrace(Scene *scene, Ray &ray, int depth) {
    Hit hit = Hit();
    Vector colour = {0,0,0};
    if(depth <= 0) return colour;
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
            colour = colour + hit_colour * diffuse * hit.what->material.kd + hit_colour * pow(specular, 128) * hit.what->material.ks;
        }
        // multiply by how transparent it is
        colour = colour + hit_colour * scene->ka * (1 - hit.what->material.t);
        colour = colour / scene->lights.size();


        // cos(theta1)
        float cos_i = max(-1.f, min(ray.direction.dot(hit.normal), 1.f));
        // compute kr by Fresnel only if transparent
        float kr = (hit.what->material.t != 0) ? fresnel(hit.what->material.ior, cos_i) : hit.what->material.r;

        // Remove speckles TODO reuse this in shadow
        Vector shift_bias = 0.001 * hit.normal;

        // only compute if reflective material
        if(hit.what->material.r != 0){
            Ray reflection_ray;
            hit.normal.reflection(ray.direction, reflection_ray.direction);
            reflection_ray.direction.normalise();

            reflection_ray.position = cos_i < 0 ? hit.position + shift_bias : hit.position + -shift_bias;

            colour = colour +  kr * raytrace(scene, reflection_ray, depth - 1);
        }
        // if kr = 1 then total internal reflection, so only reflect
        if (hit.what->material.t != 0 && kr < 1) {
            Ray refraction_ray = Ray();
            refraction_ray.direction = refract(ray.direction, hit.normal, hit.what->material.ior, cos_i);
            refraction_ray.direction.normalise();

            refraction_ray.position = cos_i < 0 ? hit.position + -shift_bias : hit.position + shift_bias;

            colour = colour + (1 - kr) * raytrace(scene, refraction_ray, depth - 1);
        }
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

// Fresnel Equations as per wikipedia.
// https://en.wikipedia.org/wiki/Fresnel_equations
float fresnel(float refractive_index, float cos_i) {
    float n1, n2;

    if (cos_i < 0) {
        // outside
        n1 = 1; // ior of air
        n2 = refractive_index; //ior of medium
        cos_i = -cos_i;
    }
    else {
        n1 = refractive_index;
        n2 = 1;
    }

    float n = n1 / n2;
    float sin_t_squared = n * n * max((1.0 - cos_i * cos_i), 0.0);
    // as per snell's law, if >1, total internal reflection
    if (sqrt(sin_t_squared) >= 1) {
        return 1;
    }
    else {
        float cos_t = sqrt(1.0 - sin_t_squared);

        // Fresnel Equations
        // p-polarised
        float per = ((n2 * cos_i) - (n1 * cos_t)) /
                    ((n2 * cos_i) + (n1 * cos_t));
        // s-polarised
        float par = ((n2 * cos_t) - (n1 * cos_i)) /
                    ((n2 * cos_t) + (n1 * cos_i));

        return (per * per + par * par) / 2;
    }
}

Vector refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i) {
    float n1, n2;
    if (cos_i < 0) {
        // outside
        n1 = 1; // ior of air
        n2 = refractive_index; //ior of medium
        cos_i = -cos_i;
    }
    else {
        n1 = refractive_index;
        n2 = 1;
        normal.negate();
    }

    float n = n1 / n2;
    float sin_t_squared = n * n * (1.0 - cos_i * cos_i);
    float cos_t = sqrt(1.0 - sin_t_squared);

    Vector refracted_ray = n * incident_ray + (n * cos_i - cos_t) * normal;
    refracted_ray.normalise();
    return refracted_ray;
}