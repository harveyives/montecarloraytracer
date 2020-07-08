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
#include "photon_map.h"
#include "alglib/stdafx.h"
#include "alglib/alglibmisc.h"

using namespace std;
using namespace alglib;

class Scene {
public:
    vector<Object *> objects;
    vector<Light *> lights;
    float ka;
    vector<Photon> photons_global;
    vector<Photon> photons_caustic;
    vector<double> points_global;
    vector<double> points_caustic;
    vector<long> tags_global;
    vector<long> tags_caustic;
    kdtree tree_global;
    kdtree tree_caustic;
    PhotonMap *pm;

    Scene(float ambient, bool generate_photon_map);

    Hit check_intersections(Ray &ray, Hit &hit);

    Vector trace(Ray &ray, int depth);

    static bool object_occluded(vector<Object *> &objects, Vertex &hit_position, Vertex &light_position);

    float fresnel(float refractive_index, float cos_i);

    Vector refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i);

    void emit(int n, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags, bool is_caustic);

    void trace_photon(Photon p, int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags);

    void
    trace_default_photon(int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags, Light *light);

    void
    trace_caustic_photon(int depth, vector<double> &points, vector<Photon> &photons, vector<long> &tags, Light *light);
};
#endif //CODE_SCENE_H
