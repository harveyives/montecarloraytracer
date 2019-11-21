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
#include "KDTree.h"
#include "photon.h"

using namespace std;

class Scene {
public:
    vector<Object *> objects;
    vector<Light *> lights;
    float ka;
    vector<Photon> photons;
    pointVec points;
    KDTree tree;

    Scene(float ambient = 1);

    Hit check_intersections(Ray &ray, Hit &hit);

    Vector compute_colour(Ray &ray, int depth);

    static bool object_occluded(vector<Object *> &objects, Vertex &hit_position, Vertex &light_position);

    float fresnel(float refractive_index, float cos_i);

    Vector refract(Vector incident_ray, Vector normal, float refractive_index, float cos_i);

    float compute_specular_component(Ray &ray, Hit &hit, Vector &light_direction) const;

    Vector get_random_direction();

    int get_random_number(int min, int max);

    Vector sample(Vertex query_pt, float radius);

    void trace_photon(Photon photon, int depth, bool first_intersection);

    void emit_photons(int n);
};
#endif //CODE_SCENE_H
