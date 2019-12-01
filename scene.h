#ifndef CODE_SCENE_H
#define CODE_SCENE_H

#include "transform.h"
#include "polymesh.h"
#include "vector.h"
#include "light.h"
#include "point_light.h"
#include "directional_light.h"
#include "plane.h"
#include "sphere.h"
#include "photon.h"
#include "alglib/stdafx.h"
#include "alglib/alglibmisc.h"

using namespace std;
using namespace alglib;

class Scene {
public:
    vector<Object *> objects;
    vector<Light *> lights;
    float ka;
    bool photon_mapping;
    vector<Photon> photons;
    vector<double> points;
    vector<long> tags;
    kdtree kdt;

    Scene(float ambient, bool mapping, bool generate_photon_map);

    Hit check_intersections(Ray &ray, Hit &hit);

    Vector compute_colour(Ray &ray, int depth);

    static bool object_occluded(vector<Object *> &objects, Vertex &hit_position, Vertex &light_position);

    float fresnel(float refractive_index, float cos_i);

    Vector refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i);

    Vector get_random_vector_in_direction(Vector direction);

    float get_random_number(int min, int max);

    Vector sample(Vertex query, int k);

    void trace_photon(Photon photon, int depth);

    void emit_photons(int n, int depth);

    Photon get_nearest_photon(Vertex query);

    vector<Photon> gather_photons(Vertex query, int k);

    photon_type get_majority_type(Vertex query);

    Vector gather_photons_with_type(Vertex query, photon_type t);

    Vector approximate_indirect(Ray &ray, Hit &hit);

    float BRDF(Hit hit, Vector incoming, Vector view);

    Vector approximate_emmissive(Ray &ray, Hit &hit);

    vector<string> split_string(string line);

    Vector get_random_point_on_hemisphere(Vector normal);

    Vector get_random_vector();

    void build_kd_tree();

    void save_map_to_file();
};
#endif //CODE_SCENE_H
