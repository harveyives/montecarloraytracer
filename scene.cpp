#include "vector.h"
#include "vertex.h"
#include "sphere.h"
#include <cmath>
#include <string>
#include <iostream>
#include <cstring>
#include "scene.h"
#include "alglib/stdafx.h"
#include "alglib/alglibmisc.h"
#include <vector>
#include "utils.h"
#include "phong.h"
#include "texture.h"
#include "area_light.h"
#include "photon_map.h"

using namespace alglib;

// Class containing all the logic for processing the scene and computing pixel colours
// including objects, lights & calculations, photon mapping ...

Scene::Scene(float ambient, bool generate_photon_map) {
    ka = ambient;

    // Adding objects:
    // Polys
    Marble *marble = new Marble({255, 255, 255}, 0, 0.9, 1, 0, 0);
    Transform *transform = new Transform(
            3.5, 0, 0, -3,
            0, 0, 3.5, -15,
            0, 3.5, 0, 50,
            0, 0, 3, 1
    );
    PolyMesh *teapot = new PolyMesh((char *) "teapot.ply", transform, marble);
    Transform *transform_cube = new Transform(
            5, 0, 0, -12,
            0, 0, 5, -12,
            0, 5, 0, 50,
            0, 0, 10, 1);
    PolyMesh *cube = new PolyMesh((char *) "cube.ply", transform_cube,
                                  new Phong({255, 255, 255}, 0.8, 0.1, 1.01, 1.3, 0.99));

    // area lights modelled as objects to leverage pre-existing code
    Transform *transform_light = new Transform(
            6, 0, 0, 0,
            0, 0, 1, 14.999,
            0, 6, 0, 40,
            0, 0, 6, 1
    );
    PolyMesh *area_light_obj = new PolyMesh((char *) "area_light.ply", transform_light,
                                            new Phong({255, 255, 255}, 0, 0, 1, 0, 0, {255, 255, 255}));

    // Spheres
    Sphere *mirror_ball = new Sphere(Vertex(8, 0, 52), 6, new Phong({0, 0, 0}, 0.2, 0.8, 1, 1));
    Sphere *glass_ball = new Sphere(Vertex(11, -8, 38), 4, new Phong({255, 255, 255}, 0.5, 0.0, 1.7, 1, 0.99));

    // Cornell Box
    Material *wall = new Phong({255, 255, 255}, 0.0, 0.6, 1);
    Material *red_wall = new Phong({255, 0, 0}, 0.0, 0.6, 1);
    Material *green_wall = new Phong({0, 255, 0}, 0.0, 0.6, 1);
    Plane *top = new Plane(Vertex(0, 15, 0), Vector(0, -1, 0), wall);
    Plane *bottom = new Plane(Vertex(0, -15, 0), Vector(0, 1, 0), wall);
    Plane *left = new Plane(Vertex(-15, 0, 0), Vector(1, 0, 0), red_wall);
    Plane *right = new Plane(Vertex(15, 0, 0), Vector(-1, 0, 0), green_wall);
    Plane *back = new Plane(Vertex(0, 0, 60), Vector(0, 0, -1), wall);

    // Adding to list
    objects.push_back(teapot);
    objects.push_back(area_light_obj);
//    objects.push_back(cube);

    objects.push_back(glass_ball);
    objects.push_back(mirror_ball);

    objects.push_back(top);
    objects.push_back(bottom);
    objects.push_back(left);
    objects.push_back(right);
    objects.push_back(back);
    // Adding light:
    // (just below ceiling)
    AreaLight *area_light = new AreaLight({0, 14, 40}, {0, -1, 0}, {255, 255, 255}, M_PI,
                                          area_light_obj->min_limit_x,
                                          area_light_obj->min_limit_z,
                                          area_light_obj->max_limit_x,
                                          area_light_obj->max_limit_z);

    // Adding to list
    lights.push_back(area_light);
    pm = new PhotonMap(objects, lights);
    if (generate_photon_map) {
        // emit photons into scene, build tree from intersections, then save to a file
        cout << "Generating new photon map...  " << endl;
        emit(50000, 50, points_global, photons_global, tags_global, false);
        emit(10000, 50, points_caustic, photons_caustic, tags_caustic, true);
        pm->build_kd_tree(points_global, tree_global, tags_global);
        pm->build_kd_tree(points_caustic, tree_caustic, tags_caustic);
        cout << "DONE" << endl;
        cout << "Writing to file... " << endl;
        pm->save_map_to_file(tree_global, "tree_global", photons_global, "photons_global");
        pm->save_map_to_file(tree_caustic, "tree_caustic", photons_caustic, "photons_caustic");
        cout << "DONE" << endl;
    } else {
        // load a prebuild map
        cout << "Loading pre-built maps... " << endl;
        pm->load_map_from_file(tree_global, "tree_global", photons_global, "photons_global");
        pm->load_map_from_file(tree_caustic, "tree_caustic", photons_caustic, "photons_caustic");
        cout << "DONE" << endl;
    }

    return;
}

// function that determines the colour of the pixel
Vector Scene::trace(Ray &ray, int depth) {
    Vector colour = {0, 0, 0};
    if (depth <= 0) return colour;
    Hit hit = Hit();
    check_intersections(ray, hit);

    Vector emissive, global, caustic;

    int k = 700;
    vector<Photon *> local_photons = pm->gather_photons(hit.position, k, tree_global, photons_global);
    int d = 200;
    vector<Photon *> local_photons_caustic = pm->gather_photons(hit.position, d, tree_caustic, photons_caustic);

    Material *m = hit.what->material;
    if (m->is_emissive()) {
        emissive += m->e;
    }
    global = pm->estimate_radiance(ray, hit, tree_global, photons_global, local_photons);
    caustic = pm->estimate_radiance(ray, hit, tree_caustic, photons_caustic, local_photons_caustic);

    // add proportions of light from difference sources, scaled to even brightness
    Vector base_colour = emissive * pow(10, 3.5) +
                         global +
                         caustic * pow(10, -1.15);

    if (hit.flag) {
        for (Light *light : lights) {
            // if not many shadow photons, don't need to test for shadows
            if (pm->in_shadow(hit, k, local_photons)) {
                if (object_occluded(objects, hit.position, light->position)) {
                    continue;
                }
            }

            //get light direction based on point if point light, otherwise get directional light
            Vector light_direction = light->get_direction(hit.position);
            light_direction.normalise();

            // only compute specular since diffuse lighting pre-calculated during radiance estimate
            colour += m->specular(ray.direction, light_direction, hit.normal, base_colour);
        }
        // multiply by how transparent it is to get a small amount of colour of the object
        colour = colour + base_colour * ka * (1 - m->t);
        colour = colour / lights.size();

        // cos(incident angle)
        float cos_i = Utils::clamp(ray.direction.dot(hit.normal), -1, 1);
        // compute kr by Fresnel only if transparent
        float kr = (m->t != 0) ? fresnel(m->ior, cos_i) : m->r;

        // Remove speckles by shifting ray
        Vector shift_bias = 0.001 * hit.normal;

        // only compute if reflective material
        if (m->is_refractive()) {
            Ray reflection_ray;
            hit.normal.reflection(ray.direction, reflection_ray.direction);
            reflection_ray.direction.normalise();
            reflection_ray.position = cos_i < 0 ? hit.position + shift_bias : hit.position + -shift_bias;

            colour += kr * trace(reflection_ray, depth - 1);
        }
        // if kr = 1 then total internal reflection, so only reflect
        if (m->is_transmissive() && kr < 1) {
            Ray refraction_ray = Ray();
            refraction_ray.direction = refract(ray.direction, hit.normal, m->ior, cos_i);
            refraction_ray.direction.normalise();

            refraction_ray.position = cos_i < 0 ? hit.position + -shift_bias : hit.position + shift_bias;

            colour += (1 - kr) * trace(refraction_ray, depth - 1);
        }
    }
    return colour;
}


void Scene::trace_photon(Photon p, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags) {
    if (depth <= 0) return;

    Hit hit = Hit();
    check_intersections(p.ray, hit);
    if (hit.flag) {
        p.ray.position = hit.position;

        Material *m = hit.what->material;
        // only store if it's a diffuse object
        if (m->kd != 0) {
            pm->store_photon(p, points, photons, tags);
        }

        // only record shadow photons in the global map
        if (p.type != "caustic") {
            // test intersections with each object and store shadow photons there
            pm->distribute_shadow_photons(p, points, photons, tags);
        }

        // Russian Roulette
        float probability = Utils::get_random_number(0, 1);
        if (probability <= m->kd) {
            // According to grammar, caustic photons terminate at diffuse bounce
            if (p.type == "caustic") return;

            // diffuse reflection
            p.ray.direction = Utils::random_direction(hit.normal, M_PI);
            p.colour = m->base_colour(hit) / m->kd;
            trace_photon(p, depth - 1, points, photons, tags);
        } else if (probability <= m->kd + m->ks) {
            // specular reflection
            hit.normal.reflection(p.ray.direction, p.ray.direction);
            p.colour = m->base_colour(hit) / m->ks;
            trace_photon(p, depth - 1, points, photons, tags);
        } else if (probability <= m->kd + m->ks + m->t) {
            // transmit
            if (p.type != "caustic") return;

            float cos_i = max(-1.f, min(p.ray.direction.dot(hit.normal), 1.f));
            Ray refraction_ray = Ray();
            refraction_ray.direction = refract(p.ray.direction, hit.normal, hit.what->material->ior, cos_i);
            refraction_ray.direction.normalise();

            Vector shift_bias = 0.001 * hit.normal;
            refraction_ray.position = cos_i < 0 ? hit.position + -shift_bias : hit.position + shift_bias;
            p.ray = refraction_ray;
            p.colour = m->base_colour(hit) / m->t;
            trace_photon(p, depth - 1, points, photons, tags);
        } else {
            // absorbed
            return;
        }
    }
}

// emit n photons from light source
void Scene::emit(int n, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags,
                 bool is_caustic) {
    for (Light *light : lights) {
        for (int i = 0; i < n; i++) {
            if (is_caustic) {
                trace_caustic_photon(depth, points, photons, tags, light);
            } else {
                trace_default_photon(depth, points, photons, tags, light);
            }
        }
        for (Photon p: photons) {
            // scale by number of photons from light
            p.colour = p.colour / n;
        }
    }
}

// generates photon in random direction based on the surface of the light
void Scene::trace_default_photon(int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags,
                                 Light *light) {
    Vector direction = Utils::random_direction(light->direction, M_PI / 2);
    Ray ray = Ray(light->get_position(), direction);
    Photon photon = Photon(ray, light->intensity, "default");
    trace_photon(photon, depth, points, photons, tags);
}

// traces caustic photons, which behave slighly differently because they're projection mapped
void Scene::trace_caustic_photon(int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags,
                                 Light *light) {
    for (Object *obj : objects) {
        // if object not transmissive, skip
        if (obj->material->t == 0) continue;

        // calculate vector to object
        Vector to_obj = obj->centre - light->position;
        float distance = to_obj.magnitude();
        to_obj.normalise();

        //calculate angle of right angle triangle from light to object
        float cone_angle = atan(obj->radius / distance);

        // use angle to generate vector within the cone from light to object, containing the object
        Vector direction = Utils::random_direction(to_obj, cone_angle);
        Ray ray = Ray(light->get_position(), direction);

        // specify caustic photon to change behaviour on types of bounce possible
        Photon photon = Photon(ray, light->intensity, "caustic");
        trace_photon(photon, depth, points, photons, tags);
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

// Fresnel Equations as per wikipedia.
// https://en.wikipedia.org/wiki/Fresnel_equations
// with advice from https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
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

Hit Scene::check_intersections(Ray &ray, Hit &hit) {
    // for each object, check if a given ray intersects, return closes
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

bool Scene::object_occluded(vector<Object *> &objects, Vertex &hit_position, Vertex &light_position) {
    Ray shadow = Ray();
    shadow.position = hit_position;
    shadow.direction = light_position - hit_position;
    // don't cast shadows from objects that are further away than the light
    float distance_to_light = shadow.direction.magnitude();
    shadow.direction.normalise();

    // shifting so that the face does not self-shadow
    shadow.direction.normalise();
    shadow.position = shadow.get_point(0.001);

    Hit obj_shadow_hit = Hit();
    for (Object *obj : objects) {
        // emmissive objects can't cast simple shadows
        if (obj->material->is_emissive()) continue;
        obj_shadow_hit.flag = false;

        obj->intersection(shadow, obj_shadow_hit);
        if (obj_shadow_hit.flag && obj_shadow_hit.t < distance_to_light)
            return true;
    }
    return false;
}