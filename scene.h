#ifndef CODE_SCENE_H
#define CODE_SCENE_H

#include "transform.h"
#include "polymesh.h"
#include "vector.h"
#include "light.h"
#include "point_light.h"
#include "directional_light.h"
#include "plane.h"

class Scene {
public:
    std::vector<Object*> objects;
    std::vector<Light*> lights;
    float ka;

    Scene(float ambient = 1) {
        ka = ambient;

        // Adding objects:

        // Polys
        Transform *transform = new Transform(
                1.5, 0, 0, 2,
                0, 0, 1.5, -2,
                0, 1.5, 0, 25.0,
                0, 0, 1.5, 1);
        PolyMesh *pm = new PolyMesh((char *) "teapot.ply", transform, Material({0, 255, 0}, 0.8, 0.1, 1, 0.4, 0));
//        Transform *transform_pyramid = new Transform(
//                5, 0, 0, -2.5,
//                0, 0, 5, -2,
//                0, 5, 0, 5.0,
//                0, 0, 5, 1);
        Transform *transform_pyramid = new Transform(
                3, 0, 0, 0,
                0, 0, 3, 0,
                0, 3, 0, 15.0,
                0, 0, 3, 1);
        PolyMesh *pyramid = new PolyMesh((char *) "pyramid.ply", transform_pyramid, Material({0, 255, 0}, 0.8, 0.1, 1.01, 1.3, 0.98));

        // Spheres
        Sphere *sphere = new Sphere(Vertex(-5, 0, 50), 12, Material({255, 255, 255}, 0.1, 0.1, 1.7, 1, 1));
        Sphere *sphere_yellow = new Sphere(Vertex(5, 10, 110), 15, Material({255, 255, 0}, 0, 0.4));
        Sphere *sphere3 = new Sphere(Vertex(4, 4, 15), 1, Material({0, 0, 255}, 0.9, 0.6, 2, 1, 0));
        Sphere *sphere4 = new Sphere(Vertex(12, 8, 35), 12, Material({255, 0, 0}, 0, 0.4, 1));

        // Cornell Box
        Material wall = Material({255, 255, 255}, 0.4, 0.3, 1);
        Plane *top = new Plane(Vertex(0,35,0),Vector(0,-1,0), wall.set_colour({0,0,255}));
        Plane *bottom = new Plane(Vertex(0,-35,0),Vector(0,1,0), wall);
        Plane *left = new Plane(Vertex(-35,0,0),Vector(1,0,0), wall.set_colour({255,0,0}));
        Plane *right = new Plane(Vertex(35,0,0),Vector(-1,0,0), wall.set_colour({0,255,0}));
        Plane *back = new Plane(Vertex(0,0,150),Vector(0,0,-1), wall);
        Plane *behind = new Plane(Vertex(0,0,-50),Vector(0,0,1), wall);

        // Adding to list
//        objects.push_back(pm);
        objects.push_back(pyramid);

//        objects.push_back(sphere);
        objects.push_back(sphere_yellow);
//        objects.push_back(sphere3);
//        objects.push_back(sphere4);

        objects.push_back(top);
        objects.push_back(bottom);
        objects.push_back(left);
        objects.push_back(right);
        objects.push_back(back);
//        objects.push_back(behind);

        // Adding lights:
        Light *l1 = new PointLight(Vertex(-10,10,-5));
//        Light *l2 = new PointLight(Vertex(-9,10,-5));
        Light *l3 = new DirectionalLight(Vector(0,0,-1));

        // Adding to list
        lights.push_back(l1);
//        lights.push_back(l3);
    };
};
#endif //CODE_SCENE_H
