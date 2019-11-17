#ifndef CODE_PHOTON_MAP_H
#define CODE_PHOTON_MAP_H

#include <random>
#include <iostream>
#include "vertex.h"
#include "photon.h"
#include "Object.h"
#include <cstdlib>

/**
 * Nanoflann KD tree library
 * https://github.com/jlblancoc/nanoflann
 */
#include "nanoflann.hpp"

using namespace std;
using namespace nanoflann;
/**
 * Nanoflann KD tree library
 * https://github.com/jlblancoc/nanoflann
 */
template<typename T>
struct PointCloud {
    typedef T coord_t;

    std::vector<Photon> pts;
};

/**
 * Nanoflann KD tree library
 * adapted from examples
 * https://github.com/jlblancoc/nanoflann
 */
template<typename Derived>
struct PointCloudAdaptor {
    typedef typename Derived::coord_t coord_t;

    const Derived &obj;

    PointCloudAdaptor(const Derived &obj_) : obj(obj_) {}

    inline const Derived &derived() const { return obj; }

    inline size_t kdtree_get_point_count() const { return derived().pts.size(); }

    inline coord_t kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return derived().pts[idx].ray.position.x;
        else if (dim == 1) return derived().pts[idx].ray.position.y;
        else return derived().pts[idx].ray.position.z;
    }

    template<class BBOX>
    bool kdtree_get_bbox(BBOX &) const { return false; }

};

class PhotonMap {
public:
    PointCloud<float> cloud;

    void trace(int n, Vertex start, vector<Object *> objects) {
        vector<Photon> photons;
        Vertex position = start;


        for (int i = 0; i < n; i++) {
            Vector direction = get_random_direction();
            Photon photon = Photon(Ray(position, direction), Vector(), photon_type::regular);

            Hit hit = Hit();
            for (Object *obj : objects) {
                Hit obj_hit = Hit();

                obj->intersection(photon.ray, obj_hit);
                if (obj_hit.flag) {
                    if (obj_hit.t < hit.t) {
                        hit = obj_hit;
                    }
                }
            }
            if (hit.flag) {
                // Russian Roulette
                float p = get_random_number(0, 100);
                if (p <= 33) {
                    // reflect
                    hit.normal.reflection(direction, photon.ray.direction);
                } else if (p <= 66) {
                    // transmit
                    // temp
                    hit.normal.reflection(direction, photon.ray.direction);
                } else {
                    // absorb
                    continue;
                }

                position = hit.position;
                //TODO maybe change this so it's not a ray
                photon.ray.position = hit.position;
                photon.colour = hit.what->material.colour;
                cloud.pts.push_back(photon);
            }
        }
    }

    Vector sample(Vertex query_pt, size_t points) {
        float test[3] = {query_pt.x, query_pt.y, query_pt.z};

        typedef PointCloudAdaptor<PointCloud<float> > KDAdaptor;
        const KDAdaptor adaptor(cloud);

        typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, KDAdaptor>, KDAdaptor, 3> tree;

        tree index(3, adaptor, KDTreeSingleIndexAdaptorParams(10));
        index.buildIndex();

        std::vector<size_t> ret_index(points);
        std::vector<float> out_dist_sqr(points);

        points = index.knnSearch(&test[0], points, &ret_index[0], &out_dist_sqr[0]);

        ret_index.resize(points);
        out_dist_sqr.resize(points);
        Vector colour = Vector();
        for (int i = 0; i < points; i++) {
            colour = colour + cloud.pts[ret_index[i]].colour;

//            cout << "\nidx["<< i << "]=" << ret_index[i] << " dist["<< i << "]=" << out_dist_sqr[i] << "\n"<< endl;
        }
        colour = colour / points;
        return colour;
    }

private:
    static Vector get_random_direction() {
        Vector direction = Vector(get_random_number(INT32_MIN, INT32_MAX), get_random_number(INT32_MIN, INT32_MAX),
                                  get_random_number(INT32_MIN, INT32_MAX));
        direction.normalise();
        return direction;
    }

    static int get_random_number(int min, int max) {
        random_device device;
        mt19937 random(device());
        uniform_int_distribution<mt19937::result_type> distribution(min, max);
        return distribution(random);
    };
};
#endif //CODE_PHOTON_MAP_H

