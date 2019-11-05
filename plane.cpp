#include "plane.h"
#include "vector.h"
#include "ray.h"
#include "hit.h"

void Plane::intersection(Ray ray, Hit &hit)
{
    ray.direction.normalise();
    normal.normalise();

    float denominator = normal.dot(ray.direction);
    if (fabs(denominator) <= 0.0001)
        return;

    Vector plane_position_vector = position.toVector();
    Vector ray_position_vector = ray.position.toVector();

    float t = -(normal.dot(ray_position_vector) + -(normal.dot(plane_position_vector))) / denominator;
    if(t < hit.t && t > 0) {
        hit.t = t;
        hit.position = ray.position + t * ray.direction;
        hit.normal = normal;
        hit.normal.normalise();
        hit.what = this;
        hit.flag = true;
    }
}