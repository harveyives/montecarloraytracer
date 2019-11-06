//
// Created by Harvey Ives on 18/10/2019.
//

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
        Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f,0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 7.0f,0.0f,0.0f,0.0f,1.0f);
        PolyMesh *pm = new PolyMesh((char *) "teapot.ply", transform, Material({0, 255, 0}, 0.6, 0.4, 1.3, true, true));

        // Spheres
        // TODO add ambient override
        Sphere *sphere = new Sphere(Vertex(0, 0, 15), 8, Material({255, 255, 255}, 0.8, 0.0, 1.3, true, true));
        Sphere *sphere2 = new Sphere(Vertex(0, 0, 400), 350, Material({255, 255, 255}, 0, 0.4, 1, false, false));
        Sphere *sphere3 = new Sphere(Vertex(-4, 6, 8), 1, Material({0, 0, 255}, 0.9, 0.6, 2, true, false));
        Sphere *sphere4 = new Sphere(Vertex(12, 8, 35), 12, Material({255, 0, 0}, 0, 0.4, 1, false, false));

        // Cornell Box
        Material wall = Material({255, 255, 255}, 0.1, 0.3, 1);
        Plane *top = new Plane(Vertex(0,35,0),Vector(0,-1,0), wall);
        Plane *bottom = new Plane(Vertex(0,-35,0),Vector(0,1,0), wall.set_colour({0,0,255}));
        Plane *left = new Plane(Vertex(-35,0,0),Vector(1,0,0), wall.set_colour({255,0,0}));
        Plane *right = new Plane(Vertex(35,0,0),Vector(-1,0,0), wall.set_colour({0,255,0}));
        Plane *back = new Plane(Vertex(0,0,50),Vector(0,0,-1), wall);

        // Adding to list
//        objects.push_back(pm);

        objects.push_back(sphere);
        objects.push_back(sphere3);
        objects.push_back(sphere4);

        objects.push_back(top);
        objects.push_back(bottom);
        objects.push_back(left);
        objects.push_back(right);
        objects.push_back(back);

        // Adding lights:
        Light *l1 = new PointLight(Vertex(-10,10,-5));
        Light *l2 = new PointLight(Vertex(-10,10,-5));
//        Light *l3 = new DirectionalLight(Vector(0,0,-1));

        // Adding to list
        lights.push_back(l1);
        lights.push_back(l2);
//        lights.push_back(l3);
    };
};
#endif //CODE_SCENE_H
