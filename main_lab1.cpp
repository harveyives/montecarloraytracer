/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab one.
 * You should not need to modify this code.
 * It creates a framebuffer, calls the drawing function and then outputs the framebuffer as a ppm file.
 *
 * On linux.bath.ac.uk:
 *
 * Compile the code using g++ -o lab1executable main_lab1.cpp framebuffer.cpp linedrawer.cpp -lm
 *
 * Execute the code using ./lab1executable
 *
 * This will produce an image file called test.ppm. You can convert this a png file for viewing using
 *
 * pbmropng test.ppm > test.png
 *
 * You are expected to replace the line drawing function with one that only uses integers and no floats.
 */

#include <math.h>

#include "framebuffer.h"
#include "linedrawer.h"

int main(int argc, char *argv[])
{
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(512,512);

  // Generate 64 radial lines around the center of the image.
  for (float i = 0; i < M_PI * 2; i += M_PI / 32.0f)
  {
    // Generate a simple 2D vector
    float x = cos(i);
    float y = sin(i);

    // draw a line
    draw_line(fb, 256 + (int)(48.0f*x), 256 + (int)(48.0f*y), 256 + (int)(240.0f*x), 256 + (int)(240.0f*y));
  }

  // output the framebuffer
  fb->writeRGBFile((char *)"test.ppm");

  // we're done.
  return 0;
}
