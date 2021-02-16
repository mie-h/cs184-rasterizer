#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
      float level = get_level(sp); //6.3  6 7  (0.3* lvl6 color + 0.7* lvl 7)
      if (level < 0) {
          level = 0;
      }
      int nearest_level = (int) round(level);
      int lower_level = (int) floor(level);
      int upper_level = (int) ceil(level);

      if (nearest_level >= mipmap.size()) {
          nearest_level = mipmap.size()-1;
      }
      if (upper_level >= mipmap.size()) {
          upper_level = mipmap.size() - 1;
      }

    if (sp.lsm == L_ZERO) {
        if (sp.psm == P_NEAREST) {
            return sample_nearest(sp.p_uv, 0);
        }
        else {
            return sample_bilinear(sp.p_uv, 0);
        }
    }
    else if (sp.lsm == L_NEAREST) {
        //cout << level << " " <<nearest_level << endl;
        if (sp.psm == P_NEAREST) {
            return sample_nearest(sp.p_uv, nearest_level);
        }
        else {
            return sample_bilinear(sp.p_uv, nearest_level);
        }

    }
    else if (sp.lsm == L_LINEAR) {
        float upper_diff = upper_level - level;
        float lower_diff = level - lower_level;
        if (sp.psm == P_NEAREST) {
            if (upper_diff == lower_diff) {
                return sample_nearest(sp.p_uv, upper_level);
            }
            else {
                return upper_diff * sample_nearest(sp.p_uv, upper_level) + sample_nearest(sp.p_uv, lower_level) * lower_diff;
            }
         
        }
        else {
            if (upper_diff == lower_diff) {
                return sample_bilinear(sp.p_uv, upper_level);
            }
            else {
                return upper_diff * sample_bilinear(sp.p_uv, upper_level) + sample_bilinear(sp.p_uv, lower_level) * lower_diff;
            }
        }
    }

// return magenta for invalid level
    return Color(1, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    Vector2D diff1 = sp.p_dx_uv - sp.p_uv;
    Vector2D diff2 = sp.p_dy_uv - sp.p_uv;
    diff1.x = diff1.x * this->width;
    diff1.y = diff1.y * this->height;
    diff2.x = diff2.x * this->width;
    diff2.y = diff2.y * this->height;

    float first = sqrt(pow(diff1.x, 2) + pow(diff1.y, 2));
    float second = sqrt(pow(diff2.x, 2) + pow(diff2.y, 2));

    float level = log2(max(first, second));
    if (level < 0) {
        return 0;
    }
    else if (level >= mipmap.size()) {
        return mipmap.size() - 1;
    }
    return level;
    //return 0;
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];

    // find point that is closest to x and y passed in
    float u = uv.x * (mip.width-1);
    float v = uv.y * (mip.height-1);
    if (level < 0) {
        level = 0;
    }
    else if (level >= mipmap.size()) {
        level = mipmap.size() - 1;
    }
    //float lvl_adjuster = pow(2, level);
    float lvl_adjuster = 1;
    float level_adj_u = u / lvl_adjuster;
    float lvl_adj_v = v / lvl_adjuster;


    // find closest 1
    //cout << "width: " << mip.width << " hieght: " << mip.height << "\n";
    int tx = (int) round(level_adj_u);
    int ty = (int) round(lvl_adj_v);
    if (tx >= width || tx < 0 || ty >= height || ty < 0) {
        cout << "tx ty outside bounds\n";
        if (tx < 0 || ty < 0) {
            cout << "under 0\n";
            return Color(1, 0, 1);
        }
        else if (tx >= width || ty >= height) {
            cout << "tx ty too big\n";
            return Color(1, 0, 1);
        }
        return Color(1, 0, 1);
    }

    // call get_texel to get the color that matches with the texel?

    //cout << "got past mip\n";
    return mip.get_texel(tx, ty);

    // return magenta for invalid level
    //commented out for task5
    //return Color(1, 0, 1);
  }

  Color Texture::lerp(float x, Color v0, Color v1) {
      return v0 + (x * v1) + ((-x) * v0);
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];

 /*   float u = uv.x * mip.width;
    float v = uv.y * mip.height;*/

    float u = uv.x * (mip.width - 1);
    float v = uv.y * (mip.height - 1);
    if (level < 0) {
        level = 0;
    }
    else if (level >= mipmap.size()) {
        level = mipmap.size() - 1;
    }
    //float lvl_adjuster = pow(2, level);
    //float lvl_adjuster = 1;

    /*float level_adj_u = u / lvl_adjuster;
    float lvl_adj_v = v / lvl_adjuster;

    int tx = (int)round(level_adj_u);
    int ty = (int)round(lvl_adj_v);*/


    float s = u - floor(u);
    float t = v - floor(v);
    Color u00 = mip.get_texel(floor(u), floor(v));
    Color u01 = mip.get_texel(floor(u), ceil(v));
    Color u10 = mip.get_texel(ceil(u), floor(v));
    Color u11 = mip.get_texel(ceil(u), ceil(v));
    Color u0 = lerp(s, u00, u10);
    Color u1 = lerp(s, u01, u11);
    Color f = lerp(t, u0, u1);
    return f;

    /////////end for task5

    // call get_texel to get the color that matches with the texel?
    //commented out for task5
    //return mip.get_texel(tx, ty);

    // return magenta for invalid level
    //commented out for task5
    //return Color(1, 0, 1);
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
