//
// Created by Harvey Ives on 18/10/2019.
//

#ifndef CODE_SCENE_H
#define CODE_SCENE_H

#include "transform.h"
#include "polymesh.h"

class Scene {
public:
    std::vector<Object*> objects;

    Scene() {
        // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
        Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f,0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 7.0f,0.0f,0.0f,0.0f,1.0f);
        PolyMesh *pm = new PolyMesh((char *) "teapot.ply", transform, Colour());
        Sphere *sphere = new Sphere(Vertex(3, 3, 10), 2);
        Sphere *sphere2 = new Sphere(Vertex(0, 3, 10), 2);

        objects.push_back(pm);
        objects.push_back(sphere2);
        objects.push_back(sphere);
    };
};
#endif //CODE_SCENE_H
