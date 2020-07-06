# Monte Carlo Raytracer
<p align="middle">
    <img width="90%" src="img/finalimage.png">
</p>

## Features
-  Support for basic .ply files
-  Bounding spheres for mesh objects
-  Photon mapping
    - Transmitted caustics
    - Photon map serialisation & loading
-  Smooth shading
-  Perlin noise texture generation
-  Emissive objects

## Description

This raytracer was built as part of a 100% coursework advanced computer graphics module that I completed during my final year at university.

The program was built almost entirely from scratch, starting with only a few basic classes supplied by my lecturer. 
The program was then built up over a series of lab submissions, followed by a final task of implementing more complicated features such as photon mapping.

The code quality is by no means exemplar. This project was written under time constraints for a 100% coursework advanced computer graphics module in my final year at university, and so there are some questionable design decisions littered throughout. 

That said, I am proud to reveal it received a perfect score of 100/100.

The practical application of this code is limited but I showcase it simply because it's by far the most frustratingly enjoyable piece of work I've ever completed. 10/10, would code again.

I plan to rewrite this entirely at some point using the skills and mistakes I have made in this program, perhaps in Rust.

If the image above interests you, then please feel free to read the accompanying report.pdf that details how I went about creating it.

## Compilation & Execution 
Compile using:  
`g++ -std=c++11 -o raytrace *.cpp alglib/*.cpp -lm -O3`

Then run with: 
`./raytrace <width> <height> <generate_photon_map>`

>eg. `./raytrace 128 128 1` should take less than 30 seconds on modern computer hardware.


