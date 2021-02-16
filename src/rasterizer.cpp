#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

      //added for task2
      int sample_rate = this->sample_rate;
      int dots_per_side = sqrt(sample_rate);
      float space_bet_dots = 1.0 / dots_per_side;
      float starting_interval = space_bet_dots / 2.0;
      for (int inner_i = 0; inner_i < dots_per_side; inner_i++) {
          for (int inner_j = 0; inner_j < dots_per_side; inner_j++) {
              sample_buffer[sample_rate * (x * height + y) + inner_i * dots_per_side + inner_j] = c;
          }
      }
      //added for task2
    
    //task 1 stuff commented out
    //sample_buffer[y * width + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  bool RasterizerImp::inside(float pt_x, float pt_y,
      float x0, float y0,
      float x1, float y1,
      float x2, float y2) {

      float L1 = -(pt_x - x0) * (y1 - y0) + (pt_y - y0) * (x1 - x0);
      float L2 = -(pt_x - x1) * (y2 - y1) + (pt_y - y1) * (x2 - x1);
      float L3 = -(pt_x - x2) * (y0 - y2) + (pt_y - y2) * (x0 - x2);

      return (L1 >= 0 && L2 >= 0 && L3 >= 0) || (L1 <= 0 && L2 <= 0 && L3 <= 0);
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {

    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
      //for every point the field (pointx, pointy):
      //   if inside(tri, pointx, pointy)
      //        call rasterize_point()

      //opt. for extra credit
      // find min and max x values to optimize
      float min_x = min({ x0, x1, x2 });
      float max_x = max({ x0, x1, x2 });

      float min_y = min({ y0, y1, y2 });
      float max_y = max({ y0, y1, y2 });


      // find min and max y values to optimize for loops

      /*for (int i = 0; i < width; i++) {
          if (i < min_x-1 || max_x+1 < i) continue;
          for (int j = 0; j < height; j++) {
              if (j < min_y-1 || max_y+1 < j) continue;
              float pt_x = i + 0.5;
              float pt_y = j + 0.5;
              bool in_triangle = inside(pt_x, pt_y, x0, y0, x1, y1, x2, y2);
              if (in_triangle) {
                  rasterize_point(pt_x, pt_y, color);
              }

          }
      }*/
    // TODO: Task 2: Update to implement super-sampled rasterization
 
      int sample_rate = this->sample_rate;
      int dots_per_side = sqrt(sample_rate);
      float space_bet_dots = 1.0 / dots_per_side;
      float starting_interval = space_bet_dots/2.0;

      for (int i = 0; i < width; i++) {
          if (i < min_x - 1 || max_x + 1 < i) continue;
          for (int j = 0; j < height; j++) {
              if (j < min_y - 1 || max_y + 1 < j) continue;

              for (int inner_i = 0; inner_i < dots_per_side; inner_i++) {
                  for (int inner_j = 0; inner_j < dots_per_side; inner_j++) {
                      float pt_x = i + space_bet_dots * inner_i + starting_interval;
                      float pt_y = j + space_bet_dots * inner_j + starting_interval;
                      bool in_triangle = inside(pt_x, pt_y, x0, y0, x1, y1, x2, y2);
                      if (in_triangle) {
                          this->sample_buffer[sample_rate * (i * height + j) + inner_i * dots_per_side + inner_j] = color;
                      }
                  }
              }
          }
      }

  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

      float min_x = min({ x0, x1, x2 });
      float max_x = max({ x0, x1, x2 });
      float min_y = min({ y0, y1, y2 });
      float max_y = max({ y0, y1, y2 });

      int sample_rate = this->sample_rate;
      int dots_per_side = sqrt(sample_rate);
      float space_bet_dots = 1.0 / dots_per_side;
      float starting_interval = space_bet_dots / 2.0;

      for (int i = 0; i < width; i++) {
          if (i < min_x - 1 || max_x + 1 < i) continue;
          for (int j = 0; j < height; j++) {
              if (j < min_y - 1 || max_y + 1 < j) continue;
              for (int inner_i = 0; inner_i < dots_per_side; inner_i++) {
                  for (int inner_j = 0; inner_j < dots_per_side; inner_j++) {
                      float pt_x = i + space_bet_dots * inner_i + starting_interval;
                      float pt_y = j + space_bet_dots * inner_j + starting_interval;
                      bool in_triangle = inside(pt_x, pt_y, x0, y0, x1, y1, x2, y2);
                      if (in_triangle) {
                          float alpha = (-(pt_x - x1) * (y2 - y1) + (pt_y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                          float beta = (-(pt_x - x2) * (y0 - y2) + (pt_y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                          float gamma = 1 - alpha - beta;
                          this->sample_buffer[sample_rate * (i * height + j) + inner_i * dots_per_side + inner_j] = alpha*c0 + beta*c1 + gamma*c2;
                      }
                  }
              }
          }
      }


  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
      
      /*
      float min_x = min({ x0, x1, x2 });
      float max_x = max({ x0, x1, x2 });
      float min_y = min({ y0, y1, y2 });
      float max_y = max({ y0, y1, y2 });

      int sample_rate = this->sample_rate;
      int dots_per_side = sqrt(sample_rate);
      float space_bet_dots = 1.0 / dots_per_side;
      float starting_interval = space_bet_dots / 2.0;

      for (int i = 0; i < width; i++) {
          if (i < min_x - 1 || max_x + 1 < i) continue;
          for (int j = 0; j < height; j++) {
              if (j < min_y - 1 || max_y + 1 < j) continue;
              for (int inner_i = 0; inner_i < dots_per_side; inner_i++) {
                  for (int inner_j = 0; inner_j < dots_per_side; inner_j++) {
                      float pt_x = i + space_bet_dots * inner_i + starting_interval;
                      float pt_y = j + space_bet_dots * inner_j + starting_interval;
                      bool in_triangle = inside(pt_x, pt_y, x0, y0, x1, y1, x2, y2);
                      if (in_triangle) {
                          float alpha = (-(pt_x - x1) * (y2 - y1) + (pt_y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                          float beta = (-(pt_x - x2) * (y0 - y2) + (pt_y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                          float gamma = 1 - alpha - beta;
                          Vector2D uv = Vector2D(u0, v0) * alpha + Vector2D(u1, v1) * beta + Vector2D(u2, v2) * gamma;

                          //uv point nearest to pt_x
                          Color color = tex.sample_nearest(uv, 0);
                          //Color color = tex.sample_bilinear(uv, 0);
                          this->sample_buffer[sample_rate * (i * height + j) + inner_i * dots_per_side + inner_j] = color;
                      }
                  }
              }
          }
      }*/
      /*SampleParams sp = SampleParams();

      sp.p_uv = ;
      sp.psm = P_NEAREST;
      sp.lsm = L_ZERO;

      tex.sample(sp);*/

    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle


      float min_x = min({ x0, x1, x2 });
      float max_x = max({ x0, x1, x2 });
      float min_y = min({ y0, y1, y2 });
      float max_y = max({ y0, y1, y2 });

      int sample_rate = this->sample_rate;
      int dots_per_side = sqrt(sample_rate);
      float space_bet_dots = 1.0 / dots_per_side;
      float starting_interval = space_bet_dots / 2.0;

      for (int i = 0; i < width; i++) {
          if (i < min_x - 1 || max_x + 1 < i) continue;
          for (int j = 0; j < height; j++) {
              if (j < min_y - 1 || max_y + 1 < j) continue;
              for (int inner_i = 0; inner_i < dots_per_side; inner_i++) {
                  for (int inner_j = 0; inner_j < dots_per_side; inner_j++) {
                      float pt_x = i + space_bet_dots * inner_i + starting_interval;
                      float pt_y = j + space_bet_dots * inner_j + starting_interval;
                      bool in_triangle = inside(pt_x, pt_y, x0, y0, x1, y1, x2, y2);
                      if (in_triangle) {
                          // getting p_uv
                          float alpha = (-(pt_x - x1) * (y2 - y1) + (pt_y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                          float beta = (-(pt_x - x2) * (y0 - y2) + (pt_y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                          float gamma = 1 - alpha - beta;
                          Vector2D uv = Vector2D(u0, v0) * alpha + Vector2D(u1, v1) * beta + Vector2D(u2, v2) * gamma;
                          SampleParams sp;
                          sp.p_uv = uv;
                          sp.lsm = this->lsm;
                          sp.psm = this->psm;

                          // getting p_dx_uv
                          if (pt_x + 1 < width) {
                              float pt_dx = pt_x + 1;
                              float alpha = (-(pt_dx - x1) * (y2 - y1) + (pt_y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                              float beta = (-(pt_dx - x2) * (y0 - y2) + (pt_y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                              float gamma = 1 - alpha - beta;
                              Vector2D dx_uv = Vector2D(u0, v0) * alpha + Vector2D(u1, v1) * beta + Vector2D(u2, v2) * gamma;
                              sp.p_dx_uv = dx_uv;
                          }

                          // getting p_dy_uv
                          if (pt_y + 1 < height) {
                              float pt_dy = pt_y + 1;
                              float alpha = (-(pt_x - x1) * (y2 - y1) + (pt_dy - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                              float beta = (-(pt_x - x2) * (y0 - y2) + (pt_dy - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                              float gamma = 1 - alpha - beta;
                              Vector2D dy_uv = Vector2D(u0, v0) * alpha + Vector2D(u1, v1) * beta + Vector2D(u2, v2) * gamma;
                              sp.p_dy_uv = dy_uv;
                          }

                          //uv point nearest to pt_x
                          Color color = tex.sample(sp);
                          //Color color = tex.sample_bilinear(uv, 0);
                          this->sample_buffer[sample_rate * (i * height + j) + inner_i * dots_per_side + inner_j] = color;
                      }
                  }
              }
          }
      }



  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    //this->sample_buffer.resize(width * height, Color::White);
    
    //mie updated for task2
    this->sample_buffer.resize(width * height * rate, Color::White);
    //
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    //mie updated for task2
    unsigned int rate = this->sample_rate;
    this->sample_buffer.resize(width * height * rate, Color::White);
    //

    //this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
      
    //updated
    // assuming 4 sub pixels per pixel
    int sample_rate = this->sample_rate;
    int dots_per_side = sqrt(sample_rate);
    float space_bet_dots = 1.0 / dots_per_side;
    float starting_interval = space_bet_dots / 2.0;

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            Color col;
            for (int inner_i = 0; inner_i < dots_per_side; inner_i++) {
                for (int inner_j = 0; inner_j < dots_per_side; inner_j++) {
                    Color new_col = this->sample_buffer[sample_rate * (i * height + j) + inner_i * dots_per_side + inner_j];
                    if (inner_i == 0 && inner_j == 0) {
                        col = new_col;
                    }
                    else {
                        col += new_col;
                    }
                }
            }
            col = col * (1.0 / sample_rate);
            for (int k = 0; k < 3; ++k) {
                this->rgb_framebuffer_target[3 * (j * width + i) + k] = (&col.r)[k] * 255;
            }
        }
    }
    // end added for task2

   /* for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }*/
  }

  Rasterizer::~Rasterizer() { }


}// CGL
