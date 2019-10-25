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

class Scene {
public:
    std::vector<Object*> objects;
    std::vector<Light*> lights;
    float ka;

    Scene(float ambient = 1) {
        ka = ambient;

        // adding objects:
        Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f,0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 7.0f,0.0f,0.0f,0.0f,1.0f);
        PolyMesh *pm = new PolyMesh((char *) "teapot.ply", transform, Vector(0,255, 0));
        Sphere *sphere = new Sphere(Vertex(4, 12, 30), 2, Vector(255,0,0));
        Sphere *sphere2 = new Sphere(Vertex(0, 0, 400), 300);
        Sphere *sphere3 = new Sphere(Vertex(5, 12, 20), 2, Vector(0,0,255));
        Sphere *sphere4 = new Sphere(Vertex(8, 8, 30), 2, Vector(255,0,0));
        objects.push_back(sphere);
        objects.push_back(sphere2);
        objects.push_back(sphere3);
        objects.push_back(sphere4);
        objects.push_back(pm);

        //adding lights:
        Light *l1 = new PointLight(Vertex(-5,5,-5));

        lights.push_back(l1);
    };
};
#endif //CODE_SCENE_H
