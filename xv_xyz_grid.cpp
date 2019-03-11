#include <vector>

#include <fstream>
#include <iostream>
#include <locale.h>
#include <math.h>
#include <iomanip>
#include <cstdint>
#include <chrono>

#include <glm/glm.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtx/intersect.hpp>
#include <cmdline/cmdline.h>

///////////////////////////////////////////////////////////////////////////////
// STB
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include <stb_image_write.h>
///////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
  #include <io.h>
  #define access    _access_s
#else
  #include <unistd.h>
#endif
typedef unsigned char uint8;

#if defined(SINGLE)
  typedef float real;
  typedef glm::vec3 vec3;
  #define REAL float
  #define SCANF_FORMAT "%f %f %f"
  #define REAL_MAX FLT_MAX
#else
  typedef double real;
  typedef glm::dvec3 vec3;
  #define REAL double
  #define SCANF_FORMAT " %lf%*[ \t,;] %lf%*[ \t,;] %lf %*[\n]"
  #define REAL_MAX DBL_MAX
#endif

#include "AABB.h"


template <typename OUT, typename T>
bool write_binary(std::ofstream& out, T& value) {
  OUT v = (OUT)value;
  return out.write((char*)&v, sizeof(v)).good();
}


bool file_exists(const std::string &Filename) {
  return access(Filename.c_str(), 0) == 0;
}

struct process_t {
  process_t() : point_count(0), bitmap_height(0), bitmap_width(0) {

  }
  std::vector<std::string>  inputs;
  int                       simplifySplit;
  bool                      simplifySplitUseRatio;
  bool                      alwaysNegate;
  bool                      generateNormal;
  AABB                      aabb;
  int                       point_count;

  std::vector<real>         bitmap;
  int                       bitmap_width;
  int                       bitmap_height;

  std::vector<real>         previous_bitmap;
  int                       previous_bitmap_width;
  int                       previous_bitmap_height;

  void swap_bitmap() {
    this->previous_bitmap.swap(this->bitmap);
    std::swap(this->previous_bitmap_height, this->bitmap_height);
    std::swap(this->previous_bitmap_width,  this->bitmap_width);
  }
};

void parse_cmd_line(int ac, char** av, process_t& process) {
  cmdline::parser cmdparser;
  cmdparser.add<int>("simplify-split",
    'a',
    "Number of split used to generate the simplification grid",
    false,
    256);

  cmdparser.add<bool>("simplify-split-use-ratio",
    'b',
    "Compute the ratio of dataset, apply it to number of split",
    false,
    true);

  cmdparser.add<bool>("always-negate",
    'f',
    "Depth value will always be set to negative value",
    false,
    true);

  cmdparser.add<bool>("generate-normal",
    'n',
    "Generate a normal map",
    false,
    false);


  cmdparser.parse_check(ac, av);

  process.inputs = cmdparser.rest();
  process.simplifySplit = cmdparser.get<int>("simplify-split");
  process.simplifySplitUseRatio = cmdparser.get<bool>("simplify-split-use-ratio");
  process.alwaysNegate = cmdparser.get<bool>("always-negate");
  process.generateNormal = cmdparser.get<bool>("generate-normal");
}


void process_xyz_to_bin(process_t& process, const std::string& inputFileName, const std::string& outputFileName) {
  if (!file_exists(outputFileName)) {
    std::ifstream input(inputFileName);
    std::ofstream output(outputFileName, std::ios::binary | std::ios::out);

    process.point_count = 0;
    for (std::string line; std::getline(input, line); ) {
      size_t comment = line.find("#");
      if (comment > 0) {
        if (comment != std::string::npos) {
          line = line.substr(0, comment);
        }
        vec3  p;
        if (std::sscanf(line.c_str(), SCANF_FORMAT, &p.x, &p.y, &p.z) == 3) {
          output.write((char*)&p, sizeof(p));
          process.point_count++;
          process.aabb.add(p);
        }
      }
    }
  }
}

void process_bin_get_aabb(process_t& process, const std::string& input_as_bin) {
  std::ifstream input(input_as_bin, std::ios::binary | std::ios::in);
  vec3 p;
  while (input.read((char*)&p, sizeof(p))) {
    process.point_count++;
    process.aabb.add(p);
  }
}

real process_bin_to_bitmap(process_t& process, const std::string& input_as_bin) {

  vec3                len    = process.aabb.len();
  int                 width  = process.bitmap_width;
  int                 height = process.bitmap_height;
  vec3                ratio  = vec3(process.bitmap_width - 1, process.bitmap_height - 1, 1) / len;
  vec3                p;
  std:: ifstream      input(input_as_bin, std::ios::binary | std::ios::in);
  std:: vector<real>  count;
  real                z_len    = len.z > 0 ? len.z : 1;
  real                z_factor = 1.0 / z_len;

  if (glm::isinf(ratio.x)) {
    ratio.x = 0;
  }
  if (glm::isinf(-ratio.y)) {
    ratio.y = 0;
  }
  if (glm::isinf(-ratio.z)) {
    ratio.z = 0;
  }
  std::vector<real> bitmap;
  bitmap.resize(process.bitmap.size(), 0);
  count.resize(process.bitmap.size(), 0);

  while (input.read((char*)&p, sizeof(p))) {
    real x = (p.x - process.aabb.mMin.x) * ratio.x;
    real y = (p.y - process.aabb.mMin.y) * ratio.y;
    real z = (p.z - process.aabb.mMin.z) * z_factor ;

    int   px = (int)x;
    int   py = (int)y;

    int   nx = px + 1;
    int   ny = py + 1;

    real nfx = (x - (real)px);
    real nfy = (y - (real)py);

    real pfx = 1 - nfx;
    real pfy = 1 - nfy;

    real s = pfx * pfy;

    count[px + py * width] += s;
    bitmap[px + py * width] += (z * s);

    s = nfx * pfy;
    if (nx < width) {
      count[nx + py * width] += s;
      bitmap[nx + py * width] += (z * s);
    }
    
    s = pfx * nfy;
    if (ny < height) {
      count[px + ny * width] += s;
      bitmap[px + ny * width] += (z * s);
    }

    s = nfx * nfy;
    if ( (ny < height) && (nx < width) ) {
      count[nx + ny * width] += s;
      bitmap[nx + ny * width] += (z * s);
    }
  }

  size_t  size = process.bitmap.size();
  size_t  empty = 0;
  for (size_t i = 0; i < size; ++i) {
    real c = count[i];
    if (c > 0) {
      process.bitmap[i] = (bitmap[i] / count[i]) * z_len + process.aabb.mMin.z;
    } else {
      empty++;
    }
  }
  return 1.0 - ((real) empty / (real) size);
}
void process_bitmap_to_png(process_t& process, const std::string& outputFile) {
  std::vector<unsigned char> data;

  data.reserve(process.bitmap.size());
  for (auto v : process.bitmap) {
    if (v == REAL_MAX) {
      data.push_back(0);
    }
    else {
      real current = (v - process.aabb.mMin.z) / (process.aabb.mMax.z - process.aabb.mMin.z);
      data.push_back((unsigned char)(current * (real)(255.0f)));
    }
  }
  stbi_write_png(outputFile.c_str(), process.bitmap_width, process.bitmap_height, 1, &(data[0]), process.bitmap_width);
  // std::vector<float> data;
  // data.resize(process.bitmap.size());
  // std::copy(process.bitmap.begin(), process.bitmap.end(), data.begin());
  // stbi_write_hdr(outputFile.c_str(), process.bitmap_width, process.bitmap_height, 1, &(data[0]));
}

void process_bitmap_negate(process_t& process) {
  if (process.alwaysNegate) {
    for (auto& v : process.bitmap) {
      if (v != REAL_MAX) {
        v = (v > 0) ? -v : v;
      }
    }
    real minZ = (process.aabb.mMin.z > 0) ? -process.aabb.mMin.z : process.aabb.mMin.z;
    real maxZ = (process.aabb.mMax.z > 0) ? -process.aabb.mMax.z : process.aabb.mMax.z;
    process.aabb.mMin.z = glm::min(minZ, maxZ);
    process.aabb.mMax.z = glm::max(minZ, maxZ);
  }
}

void process_bitmap_to_xvb(process_t& process, const std::string& outputFile) {

  real min_z = REAL_MAX;
  real max_z = -REAL_MAX;
  for (auto& v : process.bitmap) {
    if (v != REAL_MAX) {
      min_z = glm::min(min_z, v);
      max_z = glm::max(max_z, v);
    }
  }

  std::ofstream output(outputFile.c_str(), std::ios_base::binary | std::ios_base::trunc);

  write_binary<int>(output, process.bitmap_width);
  write_binary<int>(output, process.bitmap_height);

  write_binary<float>(output, process.aabb.mMin.x);
  write_binary<float>(output, process.aabb.mMin.y);
  write_binary<float>(output, min_z);

  write_binary<float>(output, process.aabb.mMax.x);
  write_binary<float>(output, process.aabb.mMax.y);
  write_binary<float>(output, max_z);

  for (auto& v : process.bitmap) {
    real d = (v == REAL_MAX) ? (process.aabb.mMin.z - 1) : v;
    write_binary<float>(output, d);
  }
}

void process_bitmap_double(process_t& process) {

  int width  = process.bitmap_width * 2;
  int height = process.bitmap_height * 2;

  std::vector<real> bitmap;

  bitmap.resize(width * height);
  // this sucks monkey balls !
  for (int y=0,dy=0; y < process.bitmap_height; ++y, dy+=2) {
    for (int x=0,dx=0; x < process.bitmap_width; ++x, dx+=2) {
      real value = process.bitmap[y * process.bitmap_width + x];
      bitmap[dy     * width + dx]     = value;
      bitmap[dy     * width + dx + 1] = value;
      bitmap[(dy+1) * width + dx]     = value;
      bitmap[(dy+1) * width + dx + 1] = value;
    }
  }

  process.swap_bitmap();

  process.bitmap.swap(bitmap);
  process.bitmap_width = width;
  process.bitmap_height= height;
}


uint8 real_to_byte(real value) {
  return (uint8)  glm::max( 
                    0.0, 
                    glm::min(
                      (value + 1.0) * (255.0 / 2.0), 
                      255.0
                    )
                  );
}

vec3 compute_normal_sobel(const std::vector<double>& image, int w, int h, int u, int v)
{
    // Value from trial & error.
    // Seems to work fine for the scales we are dealing with.
    real u0 = glm::max(u - 1, 0);
    real v0 = glm::max(v - 1, 0) * w;

    real u2 = glm::min(u + 1, w - 1);
    real v2 = glm::min(v + 1, h - 1) * w;

    real u1 = u;
    real v1 = v * w;

    real tl = glm::abs(image[u0 + v0]);
    real l  = glm::abs(image[u0 + v1]);
    real bl = glm::abs(image[u0 + v2]);
    real b  = glm::abs(image[u1 + v2]);
    real br = glm::abs(image[u2 + v2]);
    real r  = glm::abs(image[u2 + v1]);
    real tr = glm::abs(image[u2 + v0]);
    real t  = glm::abs(image[u1 + v0]);

    real dX = tr + 2 * r + br - tl - 2 * l - bl;
    real dY = bl + 2 * b + br - tl - 2 * t - tr;

    return glm::normalize(vec3(dX, dY, 32.0f));
}

void process_generate_normalmap(const process_t& process, const std::string& file_name) {
  const int width  = process.bitmap_width;
  const int height = process.bitmap_height;
  const int maxx   = width - 1;
  const int maxy   = height - 1;  

  // Matrix filter = get_gaussian(5, 5, 10.0);

  // std::vector<double> bitmap = apply_filter(process.bitmap, width, height, filter);
  const auto& bitmap = process.bitmap;
  std::vector<uint8> data;
  data.resize(width * height * 3, 255);

  for (int y = 0; y < height; y++) {
    // int y1 = glm::min(y + 1, maxy);
    // glm::vec3 p, px, py, n;

//    glm::vec3* n = normals.data() + y * width;
    for (int x = 0; x < width; x++) {
      // int x1 = glm::min(x + 1, maxx);
      // glm::vec3 p0(x,  y,  bitmap[x  + y  * width]);
      // glm::vec3 p1(x1, y, bitmap[x1 + y  * width]);
      // glm::vec3 p2(x,  y1, bitmap[x  + y1 * width]);
      
      // glm::vec3 a = p1 - p0;
      // glm::vec3 b = p2 - p0;
      // n = glm::normalize(glm::cross(a, b));
      vec3 n = compute_normal_sobel(bitmap, width, height, x, y);
      data[y*width*3 + x*3 + 0] = real_to_byte(n.x);
      data[y*width*3 + x*3 + 1] = real_to_byte(n.y);
      data[y*width*3 + x*3 + 2] = real_to_byte(n.z);
    }
  }
  
  // for (int x = 0; x < width; ++x) {
  //     for (int y = 0; y < height; ++y) {
  //       const glm::vec3 n = compute_normal(process.bitmap, width, height, x, y);
  //       // convert to rgb
  //       data[y*width*3 + x*3 + 0] = float_to_char(n.x);
  //       data[y*width*3 + x*3 + 1] = float_to_char(n.y);
  //       data[y*width*3 + x*3 + 2] = float_to_char(n.z);
  //     }
  // }
  stbi_write_png(file_name.c_str(), width, height, 3, data.data(), 0);
}

void process_normalmap(process_t& process, const std::string& input, const std::string& input_as_bin) {
  int  split            = 128;
  vec3 len              = process.aabb.len();
  real bitmap_ratio     = len.x / len.y;

  process.bitmap_width  = (int)glm::max<real>((split * bitmap_ratio), 4);
  process.bitmap_height = (int)glm::max<real>((split * (1 / bitmap_ratio)), 4);
  process.bitmap.resize(process.bitmap_width * process.bitmap_height, REAL_MAX);

  while (true) {
    std::cout << "xv_xyz_grid: process_bin_to_bitmap " << split << std::endl;

    real coverage = process_bin_to_bitmap(process, input_as_bin);

    // we plotted less than 90% of the bitmap, bail out.
    if (coverage < 0.9) {
      if (!process.previous_bitmap.empty()) {
        process.swap_bitmap();
      }
      break;
    }

    {
      std::stringstream name;
      name << input << "." << split << ".normal.png";
      process_generate_normalmap(process, name.str());
    }
    {
      std::stringstream name;
      name << input << "." << split << ".png";
      process_bitmap_to_png(process, name.str());
    }
    if (split >= 2048) {
      break;
    }
    split *= 2;
    process_bitmap_double(process);
  }
  std::cout << "xv_xyz_grid: process_bitmap_negate" << std::endl;
  process_bitmap_negate(process);

  std::cout << "xv_xyz_grid: process_generate_normalmap" << std::endl;
  process_generate_normalmap(process, input + ".normal.png");
}

void process_xvb(process_t& process, const std::string& input, const std::string& input_as_bin) {
  vec3 len           = process.aabb.len();
  real bitmap_ratio  = len.x / len.y;
  int  split         = process.simplifySplit;
  
  process.bitmap_width  = (int)glm::max<real>((split * bitmap_ratio), 4);
  process.bitmap_height = (int)glm::max<real>((split * (1 / bitmap_ratio)), 4);

  process.bitmap.resize(process.bitmap_width * process.bitmap_height, REAL_MAX);
  
  std::cout << "xv_xyz_grid: process_bin_to_bitmap" << std::endl;
  process_bin_to_bitmap(process, input_as_bin);
  std::cout << "xv_xyz_grid: process_bitmap_negate" << std::endl;
  process_bitmap_negate(process);
  std::cout << "xv_xyz_grid: process_bitmap_to_png" << std::endl;
  process_bitmap_to_png(process, input + ".before.png");
  // process_fill_bitmap(process);
  std::cout << "xv_xyz_grid: process_bitmap_to_xvb" << std::endl;
  process_bitmap_to_xvb(process, input + ".xvb");
}

int main(int ac, char **av) {
  setlocale(LC_ALL, "C");

  process_t process;
  parse_cmd_line(ac, av, process);
  for (size_t arg = 0; arg < process.inputs.size(); arg++) {
    const std::string& input        = process.inputs[arg];
    const std::string& input_as_bin = input + ".bin";
    std::cout << "xv_xyz_grid: process_xyz_to_bin" << std::endl;
    process_xyz_to_bin(process, input, input_as_bin);
    std::cout << "xv_xyz_grid: process_bin_get_aabb" << std::endl;
    process_bin_get_aabb(process, input_as_bin);

    if (process.point_count == 0) {
      std::cout << "Was not able to read any data." << std::endl;
      return -1;
    }

    if (process.generateNormal) {
      process_normalmap(process, input, input_as_bin);
    } else {
      process_xvb(process, input, input_as_bin);
    }

  }
  return 0;
}
