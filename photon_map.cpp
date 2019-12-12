#include "area_light.h"
#include "texture.h"
#include "phong.h"
#include "utils.h"
#include <vector>
#include "alglib/alglibmisc.h"
#include "alglib/stdafx.h"
#include "scene.h"
#include <cstring>
#include <iostream>
#include <string>
#include <cmath>
#include "sphere.h"
#include "vertex.h"
#include "vector.h"
#include "photon_map.h"

// class containing functions to manipulate the data in the photon map

PhotonMap::PhotonMap(vector<Object *> o, vector<Light *> l) {
    objects = o;
    lights = l;
}

// estimate radiance for a given point
Vector PhotonMap::estimate_radiance(Ray &ray, Hit &hit, kdtree &tree, vector<Photon> &photons,
                             vector<Photon *> local_photons) {
    Vector colour = Vector();
    Material *m = hit.what->material;

    // find max distance
    float max_dist = -1;
    for (Photon *p: local_photons) {
        float dist = (p->ray.position - hit.position).magnitude();
        if (dist > max_dist) max_dist = dist;
    }

    for (Photon *p: local_photons) {
        float dist = (p->ray.position - hit.position).magnitude();

        Vector base_colour = m->base_colour(hit);

        // gaussian filtering as per Henrik Jensen's recommendations.
        float alpha = 0.918;
        float beta = 1.953;
        float gaussian =
                alpha * (1 - (1 - exp(-beta * ((dist * dist) / (2 * max_dist * max_dist)))) / (1 - exp(-beta)));

        // calculate photon contribution based on brdf and photon intensity
        colour += p->colour * m->light_colour(ray.direction, p->ray.direction, hit.normal, base_colour) * gaussian;
    }
    // divide through by max volume of disc sampled within
    colour = colour / (M_PI * max_dist * max_dist);
    return colour;
}


// gather k nearest neighbour photons
vector<Photon *> PhotonMap::gather_photons(Vertex p, int k, kdtree &tree, vector<Photon> &photons) {
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

// function to test if photons surrounding hit are shadow photons
bool PhotonMap::in_shadow(const Hit &hit, int k, vector<Photon *> &local_photons) {
    int counter = 0;
    for (Photon *p : local_photons) {
        if (p->type == "shadow") counter++;
    }
    if (counter >= (k / 2)) {
        return true;
    }
    return false;
}

void
PhotonMap::distribute_shadow_photons(const Photon &p, vector<double> &points, vector<Photon> &photons,
                                     vector<long> &tags) {
    for (Object *obj : objects) {
        Hit shadow_photon_hit = Hit();
        obj->intersection(p.ray, shadow_photon_hit);
        if (shadow_photon_hit.flag) {
            Photon shadow = Photon(Ray(shadow_photon_hit.position, p.ray.direction), Vector(), "shadow");
            store_photon(shadow, points, photons, tags);
        }
    }
}

void
PhotonMap::store_photon(Photon p, vector<double> &points, vector<Photon> &photons, vector<long> &tags) {
    p.ray.direction.negate();
    photons.push_back(p);
    // finding the index of the photon, and then pushing to array so it can be accessed later
    tags.push_back(points.size() / 3);

    // store points individually due to limitation with ALGLIB
    points.push_back(p.ray.position.x);
    points.push_back(p.ray.position.y);
    points.push_back(p.ray.position.z);
}

void PhotonMap::build_kd_tree(vector<double> &points, kdtree &tree, vector<long> &tags) {
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

void PhotonMap::load_map_from_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons,
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

void PhotonMap::save_map_to_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons,
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