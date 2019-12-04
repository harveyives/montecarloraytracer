#include "vector.h"
#include "vertex.h"
#include "camera.h"
#include "sphere.h"
#include <cmath>
#include <string>
#include <iostream>
#include <cstring>
#include <random>
#include "scene.h"
#include "alglib/stdafx.h"
#include "alglib/alglibmisc.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <cfloat>
#include <vector>
#include <algorithm>
#include "vector.h"
#include "utils.h"
#include "csg.h"

using namespace alglib;

Scene::Scene(float ambient, bool mapping, bool generate_photon_map) {
    ka = ambient;
    photon_mapping = mapping;
    // Adding objects:

    // Polys
    Transform *transform = new Transform(
            1.5, 0, 0, 2,
            0, 0, 1.5, -2,
            0, 1.5, 0, 45.0,
            0, 0, 1.5, 1);
    PolyMesh *pm = new PolyMesh((char *) "teapot.ply", transform, Material({0, 255, 0}, 0.8, 0.1, 1, 0.4, 0));
//
    Transform *transform_pyramid = new Transform(
            5, 0, 0, -2.5,
            0, 0, 5, 0,
            0, 5, 0, 20.0,
            0, 0, 5, 1);
    PolyMesh *pyramid = new PolyMesh((char *) "cube.ply", transform_pyramid,
                                     Material({255, 255, 255}, 0.8, 0.1, 1.01, 1.3, 0.99));

    // Spheres
    Sphere *sphere = new Sphere(Vertex(6, -6, 36), 5, Material({255, 255, 255}, 0.1, 0.0, 1.7, 1, 0.99));
    Sphere *sphere_yellow = new Sphere(Vertex(-8, -10, 50), 5, Material({255, 255, 255}, 0.0, 0.9));
    Sphere *sphere3 = new Sphere(Vertex(-3, -8, 41), 5, Material({255, 255, 255}, 0.0, 0.8));
    Sphere *sphere4 = new Sphere(Vertex(3, -8, 39), 5, Material({255, 255, 255}, 0.0, 0.8));


    // Cornell Box
    Material wall = Material({255, 255, 255}, 0.0, 0.6, 1);
    Material red_wall = Material({255, 0, 0}, 0.0, 0.6, 1);
    Material green_wall = Material({0, 255, 0}, 0.0, 0.6, 1);
    Plane *top = new Plane(Vertex(0, 15, 0), Vector(0, -1, 0), wall);
    Plane *bottom = new Plane(Vertex(0, -15, 0), Vector(0, 1, 0), wall);
    Plane *left = new Plane(Vertex(-15, 0, 0), Vector(1, 0, 0), red_wall);
    Plane *right = new Plane(Vertex(15, 0, 0), Vector(-1, 0, 0), green_wall);
    Plane *back = new Plane(Vertex(0, 0, 60), Vector(0, 0, -1), wall);
    Plane *behind = new Plane(Vertex(0, 0, -50), Vector(0, 0, 1), wall);

    // Adding to list
//        objects.push_back(pm);
//    objects.push_back(pyramid);

    objects.push_back(sphere);
    objects.push_back(sphere_yellow);
//    Object *difference = new CSG(sphere4, sphere3, "difference");
//    objects.push_back(difference);
//    objects.push_back(sphere3);
//    objects.push_back(sphere4);

    objects.push_back(top);
    objects.push_back(bottom);
    objects.push_back(left);
    objects.push_back(right);
    objects.push_back(back);
    // Adding lights:
    Light *l1 = new PointLight(Vertex({0, 13, 40}));
    Light *l2 = new PointLight(Vertex({-3, 13, 40}));
    Light *l3 = new PointLight(Vertex({3, 13, 45}));
    Light *l4 = new PointLight(Vertex({-3, 13, 45}));
    Light *l5 = new PointLight(Vertex({10, 10, 0}));
    // Adding to list
    lights.push_back(l1);
//    lights.push_back(l2);
//    lights.push_back(l3);
//    lights.push_back(l4);
//    lights.push_back(l5);
    // TODO cleanup
    if (!photon_mapping) return;

    if (generate_photon_map) {
        cout << "Generating new photon map...  " << endl;
        emit_photons(200000, 50, points_global);
        build_kd_tree(points_global, tree_global, tags_global);
//        build_kd_tree(points_caustic, tree_caustic, tags_caustic);
        cout << "DONE" << endl;
        cout << "Writing to file... " << endl;
        save_map_to_file(tree_global, "tree_global", photons_global, "photons_global");
        cout << "DONE" << endl;
    } else {
        cout << "Loading pre-built maps... " << endl;
        load_map_from_file(tree_global, "tree_global", photons_global, "photons_global");
        cout << "DONE" << endl;
    }

    return;
}

void Scene::load_map_from_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons,
                               const char *photons_filename) {
    ifstream tree_file(tree_filename);
    stringstream buffer;
    buffer << tree_file.rdbuf();
    tree_file.close();

    //deserialise file into tree
    kdtreeunserialize(buffer, tree);

    ifstream photons_file(photons_filename, ios::in);
    string line;
    while (getline(photons_file, line)) {
        vector<string> photon_line = Utils::split_string(line);
        Vector colour = Vector(stof(photon_line[0]), stof(photon_line[1]), stof(photon_line[2]));
        Vector direction = Vector(stof(photon_line[3]), stof(photon_line[4]), stof(photon_line[5]));
        Vertex position = Vertex(stof(photon_line[6]), stof(photon_line[7]), stof(photon_line[8]));
        string type = photon_line[9];
        Photon photon = Photon(Ray(position, direction), colour, type);
        photons.push_back(photon);
    }
}

void Scene::save_map_to_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons,
                             const char *photons_filename) {
    stringstream ss;
    kdtreeserialize(tree, ss);
    ofstream tree_file;
    tree_file.open(tree_filename, ios::out);
    tree_file << ss.str();
    tree_file.close();

    //serialise photon
    ofstream photons_file(photons_filename, ios::out);
    for (size_t i = 0; i < photons.size(); i++) {
        photons_file << photons[i].colour.x << "\t" << photons[i].colour.y << "\t" << photons[i].colour.z << "\t"
                     << photons[i].ray.direction.x << "\t" << photons[i].ray.direction.y << "\t"
                     << photons[i].ray.direction.z << "\t"
                     << photons[i].ray.position.x << "\t" << photons[i].ray.position.y << "\t"
                     << photons[i].ray.position.z
                     << "\t"
                     << photons[i].type << "\r\n";
    }
    photons_file.close();
}

void Scene::build_kd_tree(vector<double> &points, kdtree &tree, vector<long> &tags) {
    // converting points so that they can be stored in tree through ALGLIB
    real_2d_array matrix;
    matrix.attach_to_ptr(points.size() / 3, 3, points.data());
    ae_int_t nx = 3;
    ae_int_t ny = 0;
    ae_int_t normtype = 2;
    real_1d_array x;
    integer_1d_array input;
    input.setcontent(points.size() / 3, tags.data());
    kdtreebuildtagged(matrix, input, nx, ny, normtype, tree);
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

    Vector base_colour;
    if (photon_mapping) {
//        base_colour = estimate_radiance(ray, hit, 5, "direct", tree_caustic, 900, photons_caustic);
        base_colour = estimate_radiance(ray, hit, 5, "indirect", tree_global, 900, photons_global);
    } else {
        base_colour = hit.what->material.colour;
    }
    if (hit.flag) {
//        // TODO change these ugly pointers
        for (Light *light : lights) {
            // if area shaded
            if (object_occluded(objects, hit.position, light->position)) {
                continue;
            }
            //get light direction based on point if point light, otherwise get directional light
            Vector light_direction = light->get_light_direction(hit.position);
            light_direction.normalise();

            if (photon_mapping) {
                Vector v = light_direction;
                v.negate();
                Vector u = ray.direction;
                u.negate();
                colour = colour + hit.what->material.get_specular_component(u, v, hit.normal,
                                                                            base_colour);
            } else {

                colour = colour + hit.what->material.compute_colour(ray.direction, light_direction, hit.normal,
                                                                    base_colour);
            }

        }
        // multiply by how transparent it is
        colour = colour + base_colour * ka * (1 - hit.what->material.t);
        colour = colour / lights.size();

//        return colour;
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

Vector
Scene::estimate_radiance(Ray &ray, Hit &hit, float cone_k, const char *type, kdtree &tree, int neighbours,
                         vector<Photon> &photons) {
    Vector colour = Vector();
    vector<Photon *> local_photons = gather_photons(hit.position, neighbours, tree, photons);
    float max_dist = -1;
    //calculate distance first
    for (Photon *p: local_photons) {
//        if (p.type != type) continue;
        float dist = (p->ray.position - hit.position).magnitude();
        if (dist > max_dist) max_dist = dist;
    }

    for (Photon *p: local_photons) {
//        if (p.type != type) continue;
        float dist = (p->ray.position - hit.position).magnitude();
        float filter = 1 - (dist / (cone_k * max_dist));
        Vector v = p->ray.direction;
        v.normalise();
        v.negate();

        Vector u = ray.direction;
        u.negate();
        colour = colour + p->colour * hit.what->material.compute_colour(u, v, hit.normal,
                                                                        hit.what->material.colour) * filter;
    }
    colour = colour / ((M_PI * max_dist * max_dist) * (1 - (2 / (3 * cone_k))));
    return colour;
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

void Scene::emit_photons(int n, int depth, vector<double> &points) {
    for (Light *light : lights) {
        for (int i = 0; i < n; i++) {
            Vector direction = Utils::random_direction({0, -1, 0}, M_PI / 2);
            Ray ray = Ray(light->position, direction);
            Photon photon = Photon(ray, light->intensity, "direct");
            trace_photon(photon, depth, points, photons_global, tags_global);
        }
        // scale by number of photons_global from light
        for (Photon p: photons_global) {
            p.colour = p.colour / n;
        }

//        for (int i = 0; i < n; i++) {
//            // TODO improve this for other lights
//            for (Object *obj : objects) {
//                if (obj->material.t == 0) continue;
//                Vector to_obj = obj->centre - light->position;
//                float distance = to_obj.magnitude();
//                to_obj.normalise();
//                float cone_angle = atan(obj->radius / distance);
//
//                Vector direction = Utils::random_direction(to_obj, cone_angle);
//                Ray ray = Ray(light->position, direction);
//                Photon photon = Photon(ray, direction, "direct");
//                trace_photon(photon, depth, points_caustic, photons_caustic, tags_caustic);
//            }
//        }
//        // scale by number of photons_global from light
//        for (Photon p: photons_caustic) {
//            p.colour = p.colour / n;
//        }
    }
}

void
Scene::trace_photon(Photon photon, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags) {
    if (depth <= 0) return;

    Hit hit = Hit();
    check_intersections(photon.ray, hit);
    if (hit.flag) {
        photon.ray.position = hit.position;
//        if (photon.type == "direct") {
//            photon.colour = hit.what->material.colour;
//        }

        Material m = hit.what->material;
        if (m.ks == 0) {
            photons.push_back(photon);
            tags.push_back(points.size() / 3);
            points.push_back(photon.ray.position.x);
            points.push_back(photon.ray.position.y);
            points.push_back(photon.ray.position.z);
        }
//        photon.type = "indirect";

        // Russian Roulette
        float p = Utils::get_random_number(0, 1);
        if (p <= m.kd) {
            //diffuse
            photon.ray.direction = Utils::random_direction(hit.normal, M_PI);
            photon.type = "indirect";
            Vector col = m.colour / m.kd;
            Photon test = Photon(photon.ray, col, "indirect");
            trace_photon(test, depth - 1, points, photons, tags);
        } else if (p <= m.kd + m.ks) {
            hit.normal.reflection(photon.ray.direction, photon.ray.direction);
            photon.type = "indirect";


            Vector col = m.colour / m.ks;
            Photon test = Photon(photon.ray, col, "indirect");
            trace_photon(test, depth - 1, points, photons, tags);
        } else if (p <= m.kd + m.ks + m.t) {
//            Vector col;
//            col.x = photon.colour.x * hit.what->material.colour.x / m.t;
//            col.y = photon.colour.y * hit.what->material.colour.y / m.t;
//            col.z = photon.colour.z * hit.what->material.colour.z / m.t;
//            Photon test = Photon(photon.ray, col, "indirect");
//            // cos(theta1)
//            float cos_i = max(-1.f, min(photon.ray.direction.dot(hit.normal), 1.f));
//            Ray refraction_ray = Ray();
//            refraction_ray.direction = refract(photon.ray.direction, hit.normal, hit.what->material.ior, cos_i);
//            refraction_ray.direction.normalise();
//
//            // Remove speckles TODO reuse this in shadow
//            Vector shift_bias = 0.001 * hit.normal;
//            refraction_ray.position = cos_i < 0 ? hit.position + -shift_bias : hit.position + shift_bias;
//            photon.ray = refraction_ray;
//            photon.type = "caustic";
//            trace_photon(photon, depth - 1, points, photons, tags);
        } else {
            // absorbed
            return;
        }
    }
}

vector<Photon *> Scene::gather_photons(Vertex p, int k, kdtree &tree, vector<Photon> &photons) {
    vector<double> point = {p.x, p.y, p.z};
    real_1d_array x;
    x.setcontent(3, point.data());

    kdtreequeryknn(tree, x, k);
    real_2d_array r = "[[]]";
    integer_1d_array output_tags = "[]";
    kdtreequeryresultstags(tree, output_tags);

    vector<Photon *> local_photons;
    for (int i = 0; i < k; i++) {
        local_photons.push_back(&photons[output_tags[i]]);
    }
    return local_photons;
}