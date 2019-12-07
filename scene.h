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
    vector<Photon> photons_global;
    vector<Photon> photons_caustic;
    vector<double> points_global;
    vector<double> points_caustic;
    vector<long> tags_global;
    vector<long> tags_caustic;
    kdtree tree_global;
    kdtree tree_caustic;

    Scene(float ambient, bool mapping, bool generate_photon_map);

    Hit check_intersections(Ray &ray, Hit &hit);

    Vector trace(Ray &ray, int depth);

    static bool object_occluded(vector<Object *> &objects, Vertex &hit_position, Vertex &light_position);

    float fresnel(float refractive_index, float cos_i);

    Vector refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i);

    void trace_photon(Photon p, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags);

    void emit(int n, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags);

    vector<Photon *> gather_photons(Vertex p, int k, kdtree &tree, vector<Photon> &photons);

    Vector estimate_radiance(Ray &ray, Hit &hit, kdtree &tree, int neighbours, vector<Photon> &photons);

    void build_kd_tree(vector<double> &points, kdtree &tree, vector<long> &tags);

    static void
    save_map_to_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons, const char *photons_filename);

    void load_map_from_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons,
                            const char *photons_filename);

    void emit_caustic(int n, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags);
};
#endif //CODE_SCENE_H
