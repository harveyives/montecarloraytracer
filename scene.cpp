#include "vector.h"
#include "vertex.h"
#include "camera.h"
#include "sphere.h"
#include <math.h>
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
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <float.h>
#include <vector>
#include <algorithm>
#include "vector.h"

using namespace alglib;

Scene::Scene(float ambient, bool generate_photon_map) {
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
    Sphere *sphere = new Sphere(Vertex(6, -10, 40), 5, Material({255, 255, 255}, 0.1, 0.1, 1.7, 1, 0));
    Sphere *sphere_yellow = new Sphere(Vertex(-8, -10, 50), 5, Material({255, 255, 255}, 0.0, 0.8));
    Sphere *sphere3 = new Sphere(Vertex(4, 4, 15), 1, Material({0, 0, 255}, 0.4, 0.6, 2, 1, 0));
    Sphere *sphere4 = new Sphere(Vertex(12, 8, 35), 12, Material({255, 0, 0}, 0, 0.4, 1));

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
//        objects.push_back(sphere3);
//        objects.push_back(sphere4);

    objects.push_back(top);
    objects.push_back(bottom);
    objects.push_back(left);
    objects.push_back(right);
    objects.push_back(back);
//        objects.push_back(behind);

    // Adding lights:
    Light *l1 = new PointLight(Vertex({3, 10, 40}));
    Light *l2 = new PointLight(Vertex({-3, 10, 40}));
    Light *l3 = new PointLight(Vertex({3, 10, 45}));
    Light *l4 = new PointLight(Vertex({-3, 10, 45}));
//        Light *l2 = new PointLight(Vertex(-9,10,-5));
//    Light *l3 = new DirectionalLight(Vector(0, 0, -1));

    // Adding to list
    lights.push_back(l1);
    lights.push_back(l2);
    lights.push_back(l3);
    lights.push_back(l4);
    // TODO cleanup
    if (generate_photon_map) {
        cout << "Generating new photon map...\n";
        emit_photons(50000, 50);
        real_2d_array matrix;
        matrix.attach_to_ptr(points.size() / 3, 3, points.data());
        ae_int_t nx = 3;
        ae_int_t ny = 0;
        ae_int_t normtype = 2;
        real_1d_array x;
        integer_1d_array input;
        input.setcontent(points.size() / 3, tags.data());
        kdtreebuildtagged(matrix, input, nx, ny, normtype, kdt);
        cout << "Writing to file...\n";
        std::stringstream ss;
        kdtreeserialize(kdt, ss);
        std::ofstream outFile;
        outFile.open("kdtree", ios::out);
        outFile << ss.str();
        outFile.close();

        std::ifstream file("kdtree");
        std::stringstream buffer;

        buffer << file.rdbuf();

        file.close();
        kdtreeunserialize(buffer, kdt);

        ofstream ofs("photons", ios::out);

        for (size_t i = 0; i < photons.size(); i++) {
            ofs << photons[i].colour.x << "\t" << photons[i].colour.y << "\t" << photons[i].colour.z << "\t"
                << photons[i].ray.direction.x << "\t" << photons[i].ray.direction.y << "\t"
                << photons[i].ray.direction.z << "\t"
                << photons[i].ray.position.x << "\t" << photons[i].ray.position.y << "\t" << photons[i].ray.position.z
                << "\t"
                << photons[i].type << "\r\n";
        }
        cout << "Done writing to file." << endl;
    } else {
        cout << "Loading pre-built map...\n";
        std::ifstream file("kdtree");
        std::stringstream buffer;

        buffer << file.rdbuf();

        file.close();
        kdtreeunserialize(buffer, kdt);

        ifstream ifs("photons", ios::in);

        string line;
        while (getline(ifs, line)) {
            vector<string> photon_line = split_string(line);
            Vector colour = Vector(stoi(photon_line[0]), stoi(photon_line[1]), stoi(photon_line[2]));
            Vector direction = Vector(stof(photon_line[3]), stof(photon_line[4]), stof(photon_line[5]));
            Vertex position = Vertex(stof(photon_line[6]), stof(photon_line[7]), stof(photon_line[8]));
            string type = photon_line[9];
            Photon photon = Photon(Ray(position, direction), colour, type);
            photons.push_back(photon);
        }
        cout << "Done loading map." << endl;
    }
}

//TODO move this to utils
vector<string> Scene::split_string(string line) {
    istringstream ss(line);
    istream_iterator<string> begin(ss), end;
    vector<string> words(begin, end);
    return words;
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

//    colour = colour + hit.what->material.colour * ka;
    Vector hit_colour = approximate_indirect(ray, hit);
    Vector indirect = approximate_indirect(ray, hit);

    if (hit.flag) {
        float max_value = max({indirect.x, indirect.y, indirect.z});
        Vector scaled_indirect = {255 * indirect.x / max_value, 255 * indirect.y / max_value,
                                  255 * indirect.z / max_value};
//        colour = colour + scaled_indirect;
//        // TODO change these ugly pointers
//
//        Vector hit_colour = hit.what->material.colour;
//
        for (Light *light : lights) {
            // if area shaded
            if (object_occluded(objects, hit.position, light->position)) {
                continue;
            }
            //get light direction based on point if point light, otherwise get directional light
            Vector light_direction = light->get_light_direction(hit.position);
            light_direction.normalise();
            colour = colour + hit.what->material.compute_colour(ray.direction, light_direction, hit.normal,
                                                                hit.what->material.colour);

        }
        // multiply by how transparent it is
//        colour = colour + hit_colour * ka * (1 - hit.what->material.t);
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


Vector Scene::approximate_emmissive(Ray &ray, Hit &hit) {
    Vector indirect_diffuse = Vector();
    vector<Photon> local_photons = gather_photons(hit.position, 800);
    float max_dist = -1;
    //calculate distance first
    for (Photon p: local_photons) {
        if (p.type != "direct") continue;
        float dist = (p.ray.position - hit.position).magnitude();
        if (dist > max_dist) max_dist = dist;
    }

    for (Photon p: local_photons) {
        if (p.type != "direct") continue;

        float alpha = 0.918;
        float beta = 1.953;
        float dist = (p.ray.position - hit.position).magnitude();
        float gaussian = alpha * (1 - (1 - exp(-beta * (dist * dist) / (2 * max_dist * max_dist))) / (1 - exp(-beta)));

        indirect_diffuse = indirect_diffuse + p.colour * gaussian;
    }

    Vector colour = indirect_diffuse;
    return colour;
}

Vector Scene::approximate_indirect(Ray &ray, Hit &hit) {
    Vector indirect_diffuse = Vector();
    vector<Photon> local_photons = gather_photons(hit.position, 100);
    float max_dist = -1;
    //calculate distance first
    for (Photon p: local_photons) {
//        if (p.type != "indirect") continue;
        float dist = (p.ray.position - hit.position).magnitude();
        if (dist > max_dist) max_dist = dist;
    }

    for (Photon p: local_photons) {
//        if (p.type != "indirect") continue;


        float alpha = 0.918;
        float beta = 1.953;
        float dist = (p.ray.position - hit.position).magnitude();
        float gaussian =
                alpha * (1 - (1 - exp(-beta * ((dist * dist) / (2 * max_dist * max_dist)))) / (1 - exp(-beta)));


        Vector v = p.ray.direction;
        v.negate();
        float diffuse = v.dot(hit.normal);
        //thus self occlusion
        if (diffuse < 0) {
            continue;
        }

//        indirect_diffuse = indirect_diffuse + hit.what->material.colour * diffuse;
//        indirect_diffuse = indirect_diffuse + hit.what->material.compute_colour(ray.direction, v, hit.normal, p.colour);
        indirect_diffuse =
                indirect_diffuse + hit.what->material.compute_colour(ray.direction, v, hit.normal, p.colour) * gaussian;
    }

    Vector colour = indirect_diffuse;
//    Vector colour = indirect_diffuse / (M_PI * max_dist * max_dist);
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

void Scene::emit_photons(int n, int depth) {
    cout << "Mapping photons..." << endl;
    for (Light *light : lights) {
        for (int i = 0; i < n; i++) {
            // TODO improve this for other lights
            Vector direction = get_random_direction();
            Ray ray = Ray(light->position, direction);
            Photon photon = Photon(ray, direction, "direct");
            trace_photon(photon, depth);
        }
        // scale by number of photons from light
        for (Photon p: photons) {
            p.colour = p.colour / n;
        }
    }
    cout << "Mapping complete." << endl;
}

void Scene::trace_photon(Photon photon, int depth) {
    if (depth <= 0) return;

    Hit hit = Hit();
    check_intersections(photon.ray, hit);
    if (hit.flag) {
        photon.ray.position = hit.position;
        if (photon.type == "direct") {
            photon.colour = hit.what->material.colour;
        }

        Material m = hit.what->material;
        if (m.ks == 0) {
            photons.push_back(photon);
            tags.push_back(points.size() / 3);
            points.push_back(photon.ray.position.x);
            points.push_back(photon.ray.position.y);
            points.push_back(photon.ray.position.z);
        }

        // Russian Roulette
        float p = get_random_number(0, 1);
        if (p <= m.kd) {
            //diffuse
            photon.ray.direction = get_random_direction(hit.normal);
            photon.colour = photon.colour;
            photon.type = "indirect";
            trace_photon(photon, depth - 1);
        } else if (p <= m.kd + m.ks) {
            hit.normal.reflection(photon.ray.direction, photon.ray.direction);
            photon.colour = photon.colour;
            photon.type = "indirect";
            trace_photon(photon, depth - 1);
        } else {
            return;
        }
    }
}

Photon Scene::get_nearest_photon(Vertex query) {
    vector<double> point = {query.x, query.y, query.z};
    real_1d_array x;
    x.setcontent(3, point.data());

    ae_int_t k = kdtreequeryknn(kdt, x, 1);
    real_2d_array r = "[[]]";
    integer_1d_array output_tags = "[]";
    kdtreequeryresultstags(kdt, output_tags);
    return photons[output_tags[0]];
}

vector<Photon> Scene::gather_photons(Vertex query, int k) {
    vector<double> point = {query.x, query.y, query.z};
    real_1d_array x;
    x.setcontent(3, point.data());

    kdtreequeryknn(kdt, x, k);
    real_2d_array r = "[[]]";
    integer_1d_array output_tags = "[]";
    kdtreequeryresultstags(kdt, output_tags);

    vector<Photon> local_photons;
    for (int i = 0; i < k; i++) {
        local_photons.push_back(photons[output_tags[i]]);
    }
    return local_photons;
}

photon_type Scene::get_majority_type(Vertex query) {
//    vector<Photon> gathered_photons = gather_photons(query, 10);
//    map<photon_type, int> map;
//    map.insert(make_pair(photon_type::direct, 0));
//    map.insert(make_pair(photon_type::indirect, 0));
//    map.insert(make_pair(photon_type::shadow, 0));
//    map.insert(make_pair(photon_type::caustic, 0));
//    for (Photon p: gathered_photons) {
//        map[p.type]++;
//    }
//
//    auto max = max_element(map.begin(), map.end(),
//                           [](const pair<photon_type, int> &a, const pair<photon_type, int> &b) {
//                               return a.second < b.second;
//                           });
//    return max->first;
}

Vector Scene::get_random_direction(Vector normal) {
    Vector random = get_random_direction();
    Vector n = normal;
    n.normalise();
//    if(random.dot(n) <= 0) random.negate();
    n.normalise();
    Vector direction = n + random;
    direction.normalise();
    return direction;
}

Vector Scene::get_random_direction() {
    Vector random = Vector(get_random_number(-10000, 10000), get_random_number(-10000, 10000),
                           get_random_number(-10000, 10000));
    random.normalise();
    return random;
}

float Scene::get_random_number(int min, int max) {
    random_device device;
    mt19937 random(device());
    uniform_real_distribution<> distribution(min, max);
    auto n = distribution(random);
    return n;
}