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

using namespace alglib;

// Class containing all the logic for processing the scene and computing pixel colours
// including objects, lights & calculations, photon mapping ...

Scene::Scene(float ambient, bool mapping, bool generate_photon_map) {
    ka = ambient;
    photon_mapping = mapping;

    // Adding objects:
    // Polys
    Texture *perlin = new Texture({255, 255, 255}, 0.0, 0.8, 1, 0, 0);
    Transform *transform = new Transform(
            3, 0, 0, -5,
            0, 0, 3, -15,
            0, 3, 0, 48,
            0, 0, 3, 1
    );
    PolyMesh *teapot = new PolyMesh((char *) "teapot.ply", transform, perlin);
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
    Sphere *mirror_ball = new Sphere(Vertex(9, -8, 55), 5, new Phong({255, 255, 255}, 0.2, 0.8, 1, 1));
    Sphere *glass_ball = new Sphere(Vertex(11, -9, 38), 4, new Phong({255, 255, 255}, 0.5, 0.0, 1.7, 1, 0.99));

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

    // if not photon mapping, skip mapping step
    if (!photon_mapping) return;

    if (generate_photon_map) {
        // emit photons into scene, build tree from intersections, then save to a file
        cout << "Generating new photon matrix...  " << endl;
        emit(50000, 50, points_global, photons_global, tags_global);
        emit_caustic(10000, 50, points_caustic, photons_caustic, tags_caustic);
        build_kd_tree(points_global, tree_global, tags_global);
        build_kd_tree(points_caustic, tree_caustic, tags_caustic);
        cout << "DONE" << endl;
        cout << "Writing to file... " << endl;
        save_map_to_file(tree_global, "tree_global", photons_global, "photons_global");
        save_map_to_file(tree_caustic, "tree_caustic", photons_caustic, "photons_caustic");
        cout << "DONE" << endl;
    } else {
        // load a prebuild map
        cout << "Loading pre-built maps... " << endl;
        load_map_from_file(tree_global, "tree_global", photons_global, "photons_global");
        load_map_from_file(tree_caustic, "tree_caustic", photons_caustic, "photons_caustic");
        cout << "DONE" << endl;
    }

    return;
}

Vector Scene::trace(Ray &ray, int depth) {
    Vector colour = {0, 0, 0};
    if (depth <= 0) return colour;
    Hit hit = Hit();
    check_intersections(ray, hit);

    Vector base_colour;

    // if photon mapping then sample radiance, otherwise use regular base colour
    Material *m = hit.what->material;
    if (m->is_emissive()) {
        base_colour += m->e * pow(10, 3.5);
    }
    if (photon_mapping) {
        base_colour += estimate_radiance(ray, hit, tree_global, 500, photons_global);
        base_colour += estimate_radiance(ray, hit, tree_caustic, 200, photons_caustic) * pow(10, -1.15);
    } else {
        base_colour = m->compute_base_colour(hit);
    }

    if (hit.flag) {
        for (Light *light : lights) {
            // if area shaded
            if (object_occluded(objects, hit.position, light->position))
                continue;

            //get light direction based on point if point light, otherwise get directional light
            Vector light_direction = light->get_direction(hit.position);
            light_direction.normalise();

            if (photon_mapping) {
                colour += m->specular(ray.direction, light_direction, hit.normal, base_colour);
            } else {
                colour += m->compute_light_colour(ray.direction, light_direction, hit.normal, base_colour);
            }
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

Vector Scene::estimate_radiance(Ray &ray, Hit &hit, kdtree &tree, int neighbours, vector<Photon> &photons) {
    Vector colour = Vector();
    vector<Photon *> local_photons = gather_photons(hit.position, neighbours, tree, photons);

    // find max distance
    float max_dist = -1;
    for (Photon *p: local_photons) {
        float dist = (p->ray.position - hit.position).magnitude();
        if (dist > max_dist) max_dist = dist;
    }

    for (Photon *p: local_photons) {
        float dist = (p->ray.position - hit.position).magnitude();

        Material *m = hit.what->material;
        Vector base_colour = hit.what->material->compute_base_colour(hit);

        // gaussian filtering as per Henrik Jensen's recommendations.
        float alpha = 0.918;
        float beta = 1.953;
        float gaussian =
                alpha * (1 - (1 - exp(-beta * ((dist * dist) / (2 * max_dist * max_dist)))) / (1 - exp(-beta)));

        // calculate photon contribution based on brdf and photon intensity
        colour += p->colour * m->compute_light_colour(ray.direction, p->ray.direction, hit.normal, base_colour) *
                  gaussian;
    }
    // divide through by max volume of sphere sampled within
    colour = colour / (M_PI * max_dist * max_dist);
    return colour;
}

void Scene::emit(int n, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags) {
    for (Light *light : lights) {
        // emit photons from light in random direction downwards from light
        for (int i = 0; i < n_global; i++) {
            Vector direction = Utils::random_direction(light->direction, M_PI / 2);
            Ray ray = Ray(light->get_position(), direction);
            Photon photon = Photon(ray, light->intensity, "default");
            trace_photon(photon, depth, points, photons, tags);
        }

        for (Photon p: photons) {
            // scale by number of photons from light
            p.colour = p.colour / n;
        }
    }
}

void Scene::emit_caustic(int n, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags) {
    for (Light *light : lights) {
        for (int i = 0; i < n; i++) {
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
        // scale by number of photons from light
        for (Photon p: photons) {
            p.colour = p.colour / n;
        }
    }
}

void Scene::trace_photon(Photon p, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags) {
    if (depth <= 0) return;

    Hit hit = Hit();
    check_intersections(p.ray, hit);
    if (hit.flag) {
        p.ray.position = hit.position;

        Material *m = hit.what->material;
        // only store if it's a diffuse object
        if (m->ks == 0) {
            p.ray.direction.negate();
            photons.push_back(p);
            tags.push_back(points.size() / 3);

            // store points individually due to limitation with ALGLIB
            points.push_back(p.ray.position.x);
            points.push_back(p.ray.position.y);
            points.push_back(p.ray.position.z);
        }

        // Russian Roulette
        float probability = Utils::get_random_number(0, 1);
        if (probability <= m->kd) {
            // According to grammar, caustic photons terminate at diffuse bounce
            if (p.type == "caustic") return;

            // diffuse reflection
            p.ray.direction = Utils::random_direction(hit.normal, M_PI);
            //store  / prob
            p.colour = m->compute_base_colour(hit) / m->kd;
            trace_photon(p, depth - 1, points, photons, tags);
        } else if (probability <= m->kd + m->ks) {
            // specular reflection
            hit.normal.reflection(p.ray.direction, p.ray.direction);
            p.colour = m->compute_base_colour(hit) / m->ks;
            trace_photon(p, depth - 1, points, photons, tags);
        } else if (probability <= m->kd + m->ks + m->t) {
            if (p.type != "caustic") return;
            // TODO use refraction method?


            float cos_i = max(-1.f, min(p.ray.direction.dot(hit.normal), 1.f));
            Ray refraction_ray = Ray();
            refraction_ray.direction = refract(p.ray.direction, hit.normal, hit.what->material->ior, cos_i);
            refraction_ray.direction.normalise();

            Vector shift_bias = 0.001 * hit.normal;
            refraction_ray.position = cos_i < 0 ? hit.position + -shift_bias : hit.position + shift_bias;
            p.ray = refraction_ray;
            p.colour = m->compute_base_colour(hit) / m->t;
            trace_photon(p, depth - 1, points, photons, tags);
        } else {
            // absorbed
            return;
        }
    }
}

// gather k nearest neighbour photons
vector<Photon *> Scene::gather_photons(Vertex p, int k, kdtree &tree, vector<Photon> &photons) {
    // boilerplate code to query kd tree from ALGLIB
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
    Vector shift_bias = 0.001 * shadow.direction;
    shadow.position = shadow.position + shift_bias;

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

void Scene::build_kd_tree(vector<double> &points, kdtree &tree, vector<long> &tags) {
    // boilerplate code
    // converting points vector so that it can be stored in tree through ALGLIB
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
