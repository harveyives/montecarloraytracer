#include "vector.h"
#include "ray.h"
#include "hit.h"
#include "csg.h"

void CSG::intersection(Ray ray, Hit &hit) {
    if (operation == "difference") {
        Hit u_hit = Hit();
        v->intersection(ray, u_hit);
        Hit v_hit = Hit();
        u->intersection(ray, v_hit);

        if (u_hit.flag && !v_hit.flag) {
            hit = u_hit;
        }
    }
}

