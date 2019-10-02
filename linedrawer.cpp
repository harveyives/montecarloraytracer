/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the code you will need to replace for Lab 1.
 *
 * It contains two simple implementations that loop over the longest axis adding the gradient to the position on the other axis at each step.
 * The objective is for you to remove the use of floating point variables and work entirely with integers.
 * You should use Bresenhams Line Algorithm to achieve this.
 */

#include <iostream>
#include "linedrawer.h"

int draw_x_line(FrameBuffer *fb, int sx, int sy, int ex, int ey)
{ 
  int x = sx;
  int y = sy;

  int dx = ex - sx;
  int dy = ey - sy;

  int xdir = 1;
  int ydir = 1;
  if(dy < 0) {
    dy = -dy;
    ydir = -1;
  }
  if (dx < 0) {
    xdir = -1;
  }

  int param = 2 * dy - dx;

  while (x != ex) {
    fb->plotPixel(x, y, 1.0f, 1.0f, 1.0f);
    if(param > 0) {
        param = param + 2 * dy - 2 * dx;
        y += ydir;
    }
    else {
        param = param + 2 * dy;
    }
    x += xdir;
  }
}

int draw_y_line(FrameBuffer *fb, int sx, int sy, int ex, int ey)
{ 
  int x = sx;
  int y = sy;

  int dx = ex - sx;
  int dy = ey - sy;

  int xdir = 1;
  int ydir = 1;
  if(dx < 0) {
    dx = -dx;
    xdir = -1;
  }
  if (dy < 0) {
    ydir = -1;
  }

  int param = 2 * dy - dx;
  while (y != ey) {
    fb->plotPixel(x, y, 1.0f, 1.0f, 1.0f);
    if(param > 0) {
        param = param + 2 * dx - 2 * dy;
        x += xdir;
    }
    else {
        param = param + 2 * dx;
    }
    y += ydir;
  }
}

int draw_line(FrameBuffer *fb, int sx, int sy, int ex, int ey)
{
  if ((sx == ex) && (sy==ey))
  {
    return fb->plotPixel(sx, sy, 1.0f, 1.0f, 1.0f);
    
  } else if (((ex-sx)* (ex-sx)) >= ((ey-sy)* (ey-sy)))
  {
    if(sx > ex) {
      return draw_x_line(fb, ex, ey, sx, sy);
    }
    else {
      return draw_x_line(fb, sx, sy, ex, ey);
    }
    
  } else
  {
    if(sy > ey) {
      return draw_y_line(fb, ex, ey, sx, sy);
    }
    else {
      return draw_y_line(fb, sx, sy, ex, ey);
    }
  }

}
