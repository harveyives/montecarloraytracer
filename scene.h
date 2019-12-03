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

    void trace_photon(Photon photon, int depth);

    void emit_photons(int n, int depth);

    vector<Photon> gather_photons(Vertex query, int k);
    Vector approximate_indirect(Ray &ray, Hit &hit);

    void build_kd_tree();

    void save_map_to_file();

    void load_map_from_file();
};
#endif //CODE_SCENE_H
