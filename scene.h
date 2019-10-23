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
        Sphere *sphere = new Sphere(Vertex(3, 3, 10), 2, Vector(255,0,0));
        Sphere *sphere2 = new Sphere(Vertex(0, 3, 10), 2);

        objects.push_back(pm);
        objects.push_back(sphere2);
        objects.push_back(sphere);

        //adding lights:
        Light *l1 = new PointLight(Vertex(15,15,0));

        lights.push_back(l1);
    };
};
#endif //CODE_SCENE_H
