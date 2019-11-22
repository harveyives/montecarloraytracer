#include "vector.h"
#include "vertex.h"
#include "camera.h"
#include "sphere.h"
#include "framebuffer.h"
#include <math.h>
#include <string>
#include <iostream>
#include <cstring>
#include <random>
#include "scene.h"
#include "alglib/stdafx.h"
#include "alglib/alglibmisc.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace alglib;

Scene::Scene(float ambient) {
    ka = ambient;
    // Adding objects:

    // Polys
//    Transform *transform = new Transform(
//            1.5, 0, 0, 2,
//            0, 0, 1.5, -2,
//            0, 1.5, 0, 25.0,
//            0, 0, 1.5, 1);
//    PolyMesh *pm = new PolyMesh((char *) "teapot.ply", transform, Material({0, 255, 0}, 0.8, 0.1, 1, 0.4, 0));
//
//    Transform *transform_pyramid = new Transform(
//            5, 0, 0, -2.5,
//            0, 0, 5, 0,
//            0, 5, 0, 20.0,
//            0, 0, 5, 1);
//    PolyMesh *pyramid = new PolyMesh((char *) "cube.ply", transform_pyramid,
//                                     Material({255, 255, 255}, 0.8, 0.1, 1.01, 1.3, 0.99));

    // Spheres
    Sphere *sphere = new Sphere(Vertex(-5, 0, 50), 12, Material({255, 255, 255}, 0.1, 0.1, 1.7, 1, 1));
    Sphere *sphere_yellow = new Sphere(Vertex(5, 10, 110), 15, Material({255, 255, 0}, 0, 0.4));
    Sphere *sphere3 = new Sphere(Vertex(4, 4, 15), 1, Material({0, 0, 255}, 0.9, 0.6, 2, 1, 0));
    Sphere *sphere4 = new Sphere(Vertex(12, 8, 35), 12, Material({255, 0, 0}, 0, 0.4, 1));

    // Cornell Box
    Material wall = Material({255, 255, 255}, 0.4, 0.3, 1);
    Plane *top = new Plane(Vertex(0, 35, 0), Vector(0, -1, 0), wall);
    Plane *bottom = new Plane(Vertex(0, -35, 0), Vector(0, 1, 0), wall);
    Plane *left = new Plane(Vertex(-35, 0, 0), Vector(1, 0, 0), wall.set_colour({255, 0, 0}));
    Plane *right = new Plane(Vertex(35, 0, 0), Vector(-1, 0, 0), wall.set_colour({0, 255, 0}));
    Plane *back = new Plane(Vertex(0, 0, 150), Vector(0, 0, -1), wall);
    Plane *behind = new Plane(Vertex(0, 0, -50), Vector(0, 0, 1), wall);

    // Adding to list
//        objects.push_back(pm);
//    objects.push_back(pyramid);

//    objects.push_back(sphere);
    objects.push_back(sphere_yellow);
//        objects.push_back(sphere3);
//        objects.push_back(sphere4);

    objects.push_back(top);
    objects.push_back(bottom);
    objects.push_back(left);
    objects.push_back(right);
    objects.push_back(back);
//        objects.push_back(behind);

    // Adding lights:
    Light *l1 = new PointLight(Vertex({0, 30, 25}));
//        Light *l2 = new PointLight(Vertex(-9,10,-5));
    Light *l3 = new DirectionalLight(Vector(0, 0, -1));

    // Adding to list
    lights.push_back(l1);
//        lights.push_back(l3);

    emit_photons(500000);
    cout << "I'm running the scene\n";
    real_2d_array a;
    a.attach_to_ptr(points.size() / 3, 3, points.data());
    ae_int_t nx = 3;
    ae_int_t ny = 0;
    ae_int_t normtype = 2;
    real_1d_array x;
    integer_1d_array input;
    input.setcontent(points.size() / 3, tags.data());
    kdtreebuildtagged(a, input, nx, ny, normtype, kdt);
}

Hit Scene::check_intersections(Ray &ray, Hit &hit) {
    for (Object *obj : objects) {
        Hit obj_hit = Hit();

        obj->intersection(ray, obj_hit);
        if (obj_hit.flag) {
            if (obj_hit.t < hit.t) {
                hit = obj_hit;
            }
        }
    }
    return hit;
}

Vector Scene::compute_colour(Ray &ray, int depth) {
    Vector colour = {0, 0, 0};
    if (depth <= 0) return colour;
    Hit hit = Hit();
    check_intersections(ray, hit);

    if (hit.flag) {
//        return sample(hit.position, 2);
//        cout<<colour.x<<colour.y<<colour.z<<endl;
        // TODO change these ugly pointers
//        colour = pm->sample(hit.position, 6);
        Vector hit_colour = hit.what->material.colour;

        for (Light *light : lights) {
            // if area shaded
            if (object_occluded(objects, hit.position, light->position)) {
                continue;
            }
            hit_colour = sample(hit.position);

            //get light direction based on point if point light, otherwise get directional light
            Vector light_direction = light->get_light_direction(hit.position);
            light_direction.normalise();

            // diffuse
            float diffuse = light_direction.dot(hit.normal);
            //thus self occlusion
            if (diffuse < 0) {
                continue;
            }

            float specular = compute_specular_component(ray, hit, light_direction);

            colour = colour + hit_colour * diffuse * hit.what->material.kd +
                     hit_colour * pow(specular, 128) * hit.what->material.ks;
        }
        // multiply by how transparent it is
        colour = colour + hit_colour * ka * (1 - hit.what->material.t);
        colour = colour / lights.size();


        // cos(theta1)
        float cos_i = max(-1.f, min(ray.direction.dot(hit.normal), 1.f));
        // compute kr by Fresnel only if transparent
        float kr = (hit.what->material.t != 0) ? fresnel(hit.what->material.ior, cos_i) : hit.what->material.r;

        // Remove speckles TODO reuse this in shadow
        Vector shift_bias = 0.001 * hit.normal;

        // only compute if reflective material
        if (hit.what->material.r != 0) {
            Ray reflection_ray;
            hit.normal.reflection(ray.direction, reflection_ray.direction);
            reflection_ray.direction.normalise();
            reflection_ray.position = cos_i < 0 ? hit.position + shift_bias : hit.position + -shift_bias;

            colour = colour + kr * compute_colour(reflection_ray, depth - 1);
        }
        // if kr = 1 then total internal reflection, so only reflect
        if (hit.what->material.t != 0 && kr < 1) {
            Ray refraction_ray = Ray();
            refraction_ray.direction = refract(ray.direction, hit.normal, hit.what->material.ior, cos_i);
            refraction_ray.direction.normalise();


            refraction_ray.position = cos_i < 0 ? hit.position + -shift_bias : hit.position + shift_bias;

            colour = colour + (1 - kr) * compute_colour(refraction_ray, depth - 1);
        }
    }
    return colour;
}

float Scene::compute_specular_component(Ray &ray, Hit &hit, Vector &light_direction) const {
    Vector reflection = Vector();
    hit.normal.reflection(light_direction, reflection);
    float specular = reflection.dot(ray.direction);
    // thus no contribution
    if (specular < 0) {
        specular = 0.0;
    }
    return specular;
}

bool Scene::object_occluded(vector<Object *> &objects, Vertex &hit_position, Vertex &light_position) {
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
        if (obj_shadow_hit.flag && obj_shadow_hit.t < distance_to_light) {
            return true;
        }
    }
    return false;
}

// Fresnel Equations as per wikipedia.
// https://en.wikipedia.org/wiki/Fresnel_equations
float Scene::fresnel(float refractive_index, float cos_i) {
    float n1, n2;

    if (cos_i < 0) {
        // outside
        n1 = 1; // ior of air
        n2 = refractive_index; //ior of medium
        cos_i = -cos_i;
    } else {
        n1 = refractive_index;
        n2 = 1;
    }

    float n = n1 / n2;
    float sin_t_squared = n * n * max((1.0 - cos_i * cos_i), 0.0);
    // as per snell's law, if >1, total internal reflection
    if (sqrt(sin_t_squared) >= 1) {
        return 1;
    } else {
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

Vector Scene::refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i) {
    float n1, n2;
    if (cos_i < 0) {
        // outside
        n1 = 1; // ior of air
        n2 = refractive_index; //ior of medium
        cos_i = -cos_i;
    } else {
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

void Scene::emit_photons(int n) {
    cout << "Mapping photons..." << endl;
    for (Light *light : lights) {
        for (int i = 0; i < n; i++) {
            // TODO improve this for other lights
            Vector direction = get_random_direction();
            Ray ray = Ray(light->position, direction);
            Photon photon = Photon(ray, Vector(), photon_type::regular);
            trace_photon(photon, 5, true);
        }
    }

    //tree = KDTree(points);
    cout << "Mapping complete." << endl;
}

void Scene::trace_photon(Photon photon, int depth, bool first_intersection = false) {
    if (depth <= 0) return;

    Hit hit = Hit();
    check_intersections(photon.ray, hit);
    if (hit.flag) {
        photon.ray.position = hit.position;
        photon.colour = hit.what->material.colour;
        if (!first_intersection) photon.type = photon_type::shadow;
        photons.push_back(photon);
        tags.push_back(points.size() / 3);
        points.push_back(photon.ray.position.x);
        points.push_back(photon.ray.position.y);
        points.push_back(photon.ray.position.z);

        // Russian Roulette
        // TODO refactor this
        int p = get_random_number(0, 100);
        if (p <= 10) {
            // absorb
            return;
        } else if (p <= 70) {
            if (hit.what->material.r != 0) {
                // reflect
                Vector reflection = Vector();
                hit.normal.reflection(photon.ray.direction, photon.ray.direction);
                photon.type = photon_type::shadow;
                trace_photon(photon, depth - 1, false);
            }
        } else {
            if (hit.what->material.r != 0) {
                // transmit
                float cos_i = max(-1.f, min(photon.ray.direction.dot(hit.normal), 1.f));
                photon.ray.direction = refract(photon.ray.direction, hit.normal, hit.what->material.ior,
                                               cos_i);
                trace_photon(photon, depth - 1, false);
            }
        }
    }
}

Vector Scene::sample(Vertex query) {
    vector<double> point = {query.x, query.y, query.z};
    real_1d_array x;
    x.setcontent(3, point.data());

    ae_int_t k = kdtreequeryknn(kdt, x, 10);
    real_2d_array r = "[[]]";
    integer_1d_array output_tags = "[]";
    kdtreequeryresultstags(kdt, output_tags);

    // TODO fix library?
    Vector colour = Vector();
    for (int i = 0; i < 10; i++) {
        colour = colour + photons[output_tags[i]].colour;
    }
    colour = colour / 10;

    return colour;
}

Vector Scene::get_random_direction() {
    Vector direction = Vector(get_random_number(INT32_MIN, INT32_MAX), get_random_number(INT32_MIN, INT32_MAX),
                              get_random_number(INT32_MIN, INT32_MAX));
    direction.normalise();
    return direction;
}

int Scene::get_random_number(int min, int max) {
    random_device device;
    mt19937 random(device());
    uniform_int_distribution<mt19937::result_type> distribution(min, max);
    return distribution(random);
}