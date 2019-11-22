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
    vector<Photon> photons;
    vector<double> points;
    vector<long> tags;
    kdtree kdt;

    Scene(float ambient = 1);

    Hit check_intersections(Ray &ray, Hit &hit);

    Vector compute_colour(Ray &ray, int depth);

    static bool object_occluded(vector<Object *> &objects, Vertex &hit_position, Vertex &light_position);

    float fresnel(float refractive_index, float cos_i);

    Vector refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i);

    float compute_specular_component(Ray &ray, Hit &hit, Vector &light_direction) const;

    Vector get_random_direction();

    int get_random_number(int min, int max);

    Vector sample(Vertex query);

    void trace_photon(Photon photon, int depth, bool first_intersection);

    void emit_photons(int n);
};
#endif //CODE_SCENE_H
