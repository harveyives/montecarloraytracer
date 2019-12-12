#ifndef CODE_PHOTON_MAP_H
#define CODE_PHOTON_MAP_H

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

class PhotonMap {
public:
    vector<Object *> objects;
    vector<Light *> lights;

    PhotonMap(vector<Object *> o, vector<Light *> l);

    vector<Photon *> gather_photons(Vertex p, int k, kdtree &tree, vector<Photon> &photons);

    Vector estimate_radiance(Ray &ray, Hit &hit, kdtree &tree, vector<Photon> &photons, vector<Photon *> local_photons);

    void build_kd_tree(vector<double> &points, kdtree &tree, vector<long> &tags);

    static void
    save_map_to_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons, const char *photons_filename);

    void load_map_from_file(kdtree &tree, const char *tree_filename, vector<Photon> &photons,
                            const char *photons_filename);

    void store_photon(Photon p, vector<double> &points, vector<Photon> &photons, vector<long> &tags);

    bool in_shadow(const Hit &hit, int k, vector<Photon *> &local_photons);

    void
    distribute_shadow_photons(const Photon &p, vector<double> &points, vector<Photon> &photons, vector<long> &tags);
};

#endif //CODE_PHOTON_MAP_H
