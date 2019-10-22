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

        // Read in the teapot model.
        PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);

        objects.push_back(pm);
    };
};
#endif //CODE_SCENE_H
