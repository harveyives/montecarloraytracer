#CM30075 Advanced Computer Graphics
The majority of the lighting and photon mapping code can be found in the scene.cpp file.
##Features 
- supports ply files
- bounding spheres for polymesh objects
- Photon mapping
    - transmitted caustics
    - supports map serialisation and saving
- Phong shading
- Perlin noise texture generation
## Compilation & Execution 
Compile using:  
>g++  -std=c++11 -o raytrace main.cpp camera.cpp framebuffer.cpp linedrawer.cpp polymesh.cpp sphere.cpp plane.cpp scene.cpp perlin.cpp alglib/*.cpp -lm -O3

Then run by: 
>./raytrace <width> <height> <generate_photon_map>

eg. `./raytracer 128 128 1` should take less than 30 seconds on modern computer hardware.


