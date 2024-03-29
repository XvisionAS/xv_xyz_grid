#include <vector>

#include <fstream>
#include <iostream>
#include <locale.h>
#include <math.h>
#include <iomanip>
#include <cstdint>
#include <chrono>
#include <functional>

#include "Externals/glm/detail/func_common.hpp"
#include "Externals/glm/glm.hpp"
#include "Externals/glm/gtc/epsilon.hpp"
#include "Externals/glm/gtx/intersect.hpp"
#include "Externals/cmdline/cmdline.h"

///////////////////////////////////////////////////////////////////////////////
// STB
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "Externals/stb_image_write.h"
///////////////////////////////////////////////////////////////////////////////

#include "Externals/zstr.hpp"

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
  typedef glm::vec2 vec2;
  #define REAL float
  #define SCANF_FORMAT "%f %f %f"
  #define REAL_MAX FLT_MAX  
#else
  typedef double real;
  typedef glm::dvec3 vec3;
  typedef glm::dvec2 vec2;
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
  int                       simplify_split_x;
  int                       simplify_split_y;
  bool                      simplify_split_use_ratio;
  bool                      always_negate;
  bool                      generate_normal;
  vec2                      resolution;
  AABB                      aabb;
  AABB                      aabb_limit;
  bool                      aabb_limit_valid;
  int                       point_count;
  real                      normal_z;
  std::vector<real>         bitmap;
  int                       bitmap_width;
  int                       bitmap_height;
  bool                      export_bbox;
  bool                      use_NODATA;
  double                    NODATA;

  std::vector<real>         previous_bitmap;
  int                       previous_bitmap_width;
  int                       previous_bitmap_height;

  void swap_bitmap() {
    this->previous_bitmap.swap(this->bitmap);
    std::swap(this->previous_bitmap_height, this->bitmap_height);
    std::swap(this->previous_bitmap_width,  this->bitmap_width);
  }
};

void print_aabb(const std::string& pre, const AABB& aabb) {
  std::cout << pre << std::fixed << std::setprecision(4);
  std::cout << " ("<< aabb.min.x << ", " << aabb.min.y;
  std::cout << ") -> ("; 
  std::cout << aabb.max.x << ", " << aabb.max.y << ")" << std::endl;
}

template<typename T>
void split_array(const std::string& string, char delimiter, std::function<T(const char*)> converter, std::vector<T>& output) {
  std::string       element;
  std::stringstream stream(string);
  while (std::getline(stream, element, delimiter)) {
    output.push_back(converter(element.c_str()));
  }
}

void parse_cmd_line(int ac, char** av, process_t& process) {
  cmdline::parser cmdparser;
  cmdparser.add<std::string>("simplify-split",
    'a',
    "Number of split used to generate the simplification grid",
    false,
    "256");

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

  cmdparser.add<std::string>("view",
    'v',
    "Limit export to AABB format : min_x:min_y:max_x:max_y",
    false
  );

  cmdparser.add<bool>("export-bbox",
    'e',
    "Export bounding box as JSON, format { min_x, min_y, min_z, max_x, max_y, max_z}",
    false,
    true
  );

  cmdparser.add<double>("no-data",
    'd',
    "Define NODATA value",
    false,
    0
  );

  cmdparser.parse_check(ac, av);

  process.inputs                   = cmdparser.rest();
  
  process.use_NODATA = cmdparser.exist("no-data");
  if (process.use_NODATA) {
    process.NODATA  = cmdparser.get<double>("no-data");
    std::cout << "Using NODATA " << process.NODATA << std::endl;
  }
  std::string      simplify_split       = cmdparser.get<std::string>("simplify-split");
  std::vector<int> splits;

  split_array<int>(simplify_split, ':', atoi, splits);

  if (splits.size() == 2) {
    process.simplify_split_x         = splits[0];
    process.simplify_split_y         = splits[1];
    process.simplify_split_use_ratio = false;
  } else {
    process.simplify_split_x         = splits[0];
    process.simplify_split_y         = process.simplify_split_x;
    process.simplify_split_use_ratio = cmdparser.get<bool>("simplify-split-use-ratio");    
  }

  if (process.simplify_split_use_ratio) {
    std::cout << "Using rationed output size x = " << process.simplify_split_x << std::endl;
  } else {
    std::cout << "Using output size x = " << process.simplify_split_x << " y = " << process.simplify_split_y << std::endl;
  }


  process.always_negate            = cmdparser.get<bool>("always-negate");
  process.generate_normal          = cmdparser.get<bool>("generate-normal");
  process.export_bbox              = cmdparser.get<bool>("export-bbox");

  std::string               view   = cmdparser.get<std::string>("view");
  std::vector<double>       aabb;
  std::string               element;

  split_array<double>(view, ':', atof, aabb);

  process.aabb_limit_valid = (aabb.size() >= 4);
  if (process.aabb_limit_valid) {
    int offset = (aabb.size() == 4) ? 0 : 1;
    process.aabb_limit.add(vec3(aabb[0], aabb[1], -1));
    process.aabb_limit.add(vec3(aabb[2 + offset], aabb[3 + offset], 1));
    print_aabb("AABB Limit", process.aabb_limit);
  }
}


void process_xyz_to_bin(process_t& process, const std::string& inputFileName, const std::string& outputFileName) {
  if (!file_exists(outputFileName)) {
    std::ifstream input(inputFileName);
    zstr::ofstream output(outputFileName, std::ios::binary | std::ios::out);

    for (std::string line; std::getline(input, line); ) {
      size_t comment = line.find("#");
      if (comment > 0) {
        if (comment != std::string::npos) {
          line = line.substr(0, comment);
        }
        vec3  p;
        const real epsilon = std::numeric_limits<real>::epsilon();
        if (std::sscanf(line.c_str(), SCANF_FORMAT, &p.x, &p.y, &p.z) == 3) {

          // maybe glm::epsilonNotEqual(process.NODATA, p.z, epsilon) ?
          if (!process.use_NODATA || process.NODATA != p.z) {
            output.write((char*)&p, sizeof(p));
          }
        }
      }
    }

  }
}

void process_bin_get_resolution(process_t& process, const std::string& input_as_bin) {    
  const std::string input_as_bin_res = input_as_bin + ".res";

  if (file_exists(input_as_bin_res)) {
    std::ifstream input(input_as_bin_res, std::ios::binary | std::ios::in);
    input.read((char*)&process.resolution, sizeof(process.resolution));
  } else {    
    zstr::ifstream i1(input_as_bin, std::ios::binary | std::ios::in);  
    vec3           p1, p2;
    process.resolution.x = REAL_MAX;
    process.resolution.y = REAL_MAX;

    const real epsilon = std::numeric_limits<real>::epsilon();
    
    // this is real slow, and we hope that the data is well formatted
    // ( so that the point are already sorted, basically, which is the case for 99.5% of our files
    // other wise, the user might have to input the correct resolution

    while (i1.read((char*)&p1, sizeof(p1)) && i1.read((char*)&p2, sizeof(p2))) {
      // zstr::ifstream i2(input_as_bin, std::ios::binary | std::ios::in);

      // while (i2.read((char*)&p2, sizeof(p2))) {
      real dx = glm::abs(p1.x - p2.x);
      real dy = glm::abs(p1.y - p2.y);
      
      
      if (dx > epsilon) {
        process.resolution.x = glm::min(process.resolution.x, dx);  
      }
      if (dy > epsilon) {
        process.resolution.y = glm::min(process.resolution.y, dy);
      }
        
      // }
    }

    std::ofstream output(input_as_bin_res, std::ios::binary | std::ios::out);
    output.write((char*)&process.resolution, sizeof(process.resolution));
  }
  std::cout << "Resolution: " << process.resolution.x << ", " << process.resolution.y << std::endl;
}

void process_bin_get_aabb(process_t& process, const std::string& input_as_bin) {
  zstr::ifstream input(input_as_bin, std::ios::binary | std::ios::in);
  vec3           p;
  process.point_count = 0;
  while (input.read((char*)&p, sizeof(p))) {
    if (process.aabb_limit_valid && !process.aabb_limit.contains2d(p)) {
      continue;
    }
    
    process.point_count++;
    process.aabb.add(p);
  }
  
  if (process.aabb.min.x > process.aabb_limit.min.x) {
    process.aabb.min.x -= process.resolution.x;
  }
  if (process.aabb.min.y > process.aabb_limit.min.y) {
    process.aabb.min.y -= process.resolution.y;
  }
  if (process.aabb.max.x < process.aabb_limit.max.x) {
    process.aabb.max.x += process.resolution.x;
  }
  if (process.aabb.max.y < process.aabb_limit.max.y) {
    process.aabb.max.y += process.resolution.y;
  }
}

real process_bin_to_bitmap(process_t& process, const std::string& input_as_bin) {
  vec3                len    = process.aabb.len();
  int                 width  = process.bitmap_width;
  int                 height = process.bitmap_height;
  vec3                ratio  = vec3(process.bitmap_width - 1, process.bitmap_height - 1, 1) / len;
  vec3                p;
  zstr::ifstream      input(input_as_bin, std::ios::binary | std::ios::in);
  std::vector<real>   count;
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
    if (process.aabb_limit_valid && !process.aabb.contains2d(p)) {
      continue;
    }
    real x = (p.x - process.aabb.min.x) * ratio.x;
    real y = (p.y - process.aabb.min.y) * ratio.y;
    real z = (p.z - process.aabb.min.z) * z_factor ;

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

  size_t  size  = 0;
  size_t  empty = 0;

  const real limits = 0;
  for (int y = 0; y < height; ++y) {
    int offset  = y * width;

    int start_x = 0;
    int end_x   = width - 1;

    for (; (start_x < width) && (count[offset + start_x] < limits); ++start_x);
    for (; (end_x > start_x) && (count[offset + end_x]   < limits); --end_x  );

    if (start_x < end_x) {
      for(int x = start_x; x <= end_x; ++x) {
        int  i = offset + x;
        real c = count[i];
        if (c > 0) {
          process.bitmap[i] = (bitmap[i] / count[i]) * z_len + process.aabb.min.z;
        }
        if (c < 1) {
          empty++;
        }
      }
      size += (end_x - start_x) + 1;
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
      real current = (v - process.aabb.min.z) / (process.aabb.max.z - process.aabb.min.z);
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
  if (process.always_negate) {
    for (auto& v : process.bitmap) {
      if (v != REAL_MAX) {
        v = (v > 0) ? -v : v;
      }
    }
    real minZ = (process.aabb.min.z > 0) ? -process.aabb.min.z : process.aabb.min.z;
    real maxZ = (process.aabb.max.z > 0) ? -process.aabb.max.z : process.aabb.max.z;
    process.aabb.min.z = glm::min(minZ, maxZ);
    process.aabb.max.z = glm::max(minZ, maxZ);
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

  write_binary<float>(output, process.aabb.min.x);
  write_binary<float>(output, process.aabb.min.y);
  write_binary<float>(output, min_z);

  write_binary<float>(output, process.aabb.max.x);
  write_binary<float>(output, process.aabb.max.y);
  write_binary<float>(output, max_z);

  for (auto& v : process.bitmap) {
    real d = (v == REAL_MAX) ? (process.aabb.min.z - 1) : v;
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


inline uint8 real_to_byte(real value) {
  return (uint8)  glm::max(
                    0.0,
                    glm::min(
                      (value + 1.0) * (255.0 / 2.0),
                      255.0
                    )
                  );
}

vec3 compute_normal_sobel(const process_t& process, const std::vector<real>& image, int w, int h, int u, int v)
{
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

  return glm::normalize(vec3(dX, dY, process.normal_z));
}

bool compute_normal(vec3& normal, const process_t& process, const std::vector<real>& image, real stepx, real stepy, int x, int y, int dx, int dy, int width) {
  int ax = x, ay = y;
  int bx = x + dx, by = y;
  int cx = x, cy = y + dy;

  real af = image[ax + ay * width];
  real bf = image[bx + by * width];
  real cf = image[cx + cy * width];

  if (af == REAL_MAX || bf == REAL_MAX || cf == REAL_MAX) {
    normal = vec3(0, 0, 1);
    return false;
  }

  vec3 a(ax * stepx, ay * stepy, af);
  vec3 b(bx * stepx, by * stepy, bf);
  vec3 c(cx * stepx, cy * stepy, cf);

  vec3 ab = b - a;
  vec3 ac = c - a;

  normal = glm::normalize(
    glm::cross(ab, ac)
  );
  return true;
}

void process_generate_normalmap(const process_t& process, const std::string& file_name) {
  const int width  = process.bitmap_width;
  const int height = process.bitmap_height;
  const int maxx   = width - 1;
  const int maxy   = height - 1;
  const real stepx = (process.aabb.max.x - process.aabb.min.x) / (real) width;
  const real stepy = (process.aabb.max.y - process.aabb.min.y) / (real) height;
  // Matrix filter = get_gaussian(5, 5, 10.0);

  // std::vector<double> bitmap = apply_filter(process.bitmap, width, height, filter);
  const auto& bitmap = process.bitmap;
  std::vector<uint8> data;
  data.resize(width * height * 4, 255);
/////////////////////////////////////////////////////////////////////////////
#define OUTPUT_NORMAL(x, y, dx, dy) {\
  vec3 n; \
  bool v = compute_normal(\
    n, \
    process, bitmap, stepx, stepy, \
    x, y, dx, dy, \
    width \
  ); \
  data[y*width*4 + x*4 + 0] = real_to_byte(n.x); \
  data[y*width*4 + x*4 + 1] = real_to_byte(n.y); \
  data[y*width*4 + x*4 + 2] = real_to_byte(n.z); \
  data[y*width*4 + x*4 + 3] = v ? 255 : 0; \
}\

/////////////////////////////////////////////////////////////////////////////
  int y;
  for (y = 0; y < height - 1; y++) {
    int x;
    for (x = 0; x < width - 1; x++) {
      OUTPUT_NORMAL(x, y, 1, 1);
    }
  }

  int x = width - 1;
  for (y = 0; y < height; y++) {
    data[y*width*4 + x*4 + 0] = data[y*width*4 + (x - 1)*4 + 0];
    data[y*width*4 + x*4 + 1] = data[y*width*4 + (x - 1)*4 + 1];
    data[y*width*4 + x*4 + 2] = data[y*width*4 + (x - 1)*4 + 2];
    data[y*width*4 + x*4 + 3] = data[y*width*4 + (x - 1)*4 + 3];
  }
  y = height - 1;
  for (x = 0; x < width; ++x) {
    data[y*width*4 + x*4 + 0] = data[(y-1)*width*4 + x *4 + 0];
    data[y*width*4 + x*4 + 1] = data[(y-1)*width*4 + x *4 + 1];
    data[y*width*4 + x*4 + 2] = data[(y-1)*width*4 + x *4 + 2];
    data[y*width*4 + x*4 + 3] = data[(y-1)*width*4 + x *4 + 3];
  }


/////////////////////////////////////////////////////////////////////////////
#undef OUTPUT_NORMAL
/////////////////////////////////////////////////////////////////////////////

  stbi_write_png(file_name.c_str(), width, height, 4, data.data(), 0);
}

void process_normalmap(process_t& process, const std::string& input, const std::string& input_as_bin) {
  int  split            = 128;
  vec3 len              = process.aabb.len();

  int  split_x       = process.simplify_split_x;
  int  split_y       = process.simplify_split_y;

  if (process.simplify_split_use_ratio) {
    real bitmap_ratio  = len.x / len.y;

    process.bitmap_width  = (int)glm::max<real>((split_x * bitmap_ratio), 4);
    process.bitmap_height = (int)glm::max<real>((split_y * (1 / bitmap_ratio)), 4);
  } else {
    process.bitmap_width  = (int)glm::max<real>(split_x, 4);
    process.bitmap_height = (int)glm::max<real>(split_y, 4);
  }

  process.bitmap.resize(process.bitmap_width * process.bitmap_height, REAL_MAX);

  while (true) {
    std::cout << "process_bin_to_bitmap " << split << std::endl;

    real coverage = process_bin_to_bitmap(process, input_as_bin);

    // we plotted less than 99% of the bitmap, bail out.
    std::cout << "computed coverage " << coverage << std::endl;
    if (coverage < 0.2) {
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
    split = split * 2;
    process_bitmap_double(process);
  }
  std::cout << "process_bitmap_negate" << std::endl;
  process_bitmap_negate(process);

  std::cout << "process_generate_normalmap" << std::endl;
  process_generate_normalmap(process, input + ".normal.png");
}

void process_export_bbox(process_t& process, const std::string& input, const std::string& input_as_bin) {
  std::string out_file = input + ".json";
  std::ofstream output(out_file.c_str(), std::ios_base::trunc);
  output << "{" << std::endl << std::fixed;
  output << "\t\"min_x\":" << process.aabb.min.x << "," << std::endl;
  output << "\t\"min_y\":" << process.aabb.min.y << "," << std::endl;
  output << "\t\"min_z\":" << process.aabb.min.z << "," << std::endl;
  output << "\t\"max_x\":" << process.aabb.max.x << "," << std::endl;
  output << "\t\"max_y\":" << process.aabb.max.y << "," << std::endl;
  output << "\t\"max_z\":" << process.aabb.max.z << std::endl;
  output << "}" << std::endl;
}

void process_xvb(process_t& process, const std::string& input, const std::string& input_as_bin) {
  vec3 len           = process.aabb.len();
  
  int  split_x       = process.simplify_split_x;
  int  split_y       = process.simplify_split_y;

  if (process.simplify_split_use_ratio) {
    real bitmap_ratio  = len.x / len.y;

    process.bitmap_width  = (int)glm::max<real>((split_x * bitmap_ratio), 4);
    process.bitmap_height = (int)glm::max<real>((split_y * (1 / bitmap_ratio)), 4);
  } else {
    process.bitmap_width  = (int)glm::max<real>(split_x, 4);
    process.bitmap_height = (int)glm::max<real>(split_y, 4);
  }


  process.bitmap.resize(process.bitmap_width * process.bitmap_height, REAL_MAX);

  std::cout << "process_bin_to_bitmap" << std::endl;
  process_bin_to_bitmap(process, input_as_bin);
  std::cout << "process_bitmap_negate" << std::endl;
  process_bitmap_negate(process);
  std::cout << "process_bitmap_to_png" << std::endl;
  process_bitmap_to_png(process, input + ".before.png");
  // process_fill_bitmap(process);
  std::cout << "process_bitmap_to_xvb" << std::endl;
  process_bitmap_to_xvb(process, input + ".xvb");
}

int main(int ac, char **av) {
  setlocale(LC_ALL, "C");

  process_t process;
  parse_cmd_line(ac, av, process);
  for (size_t arg = 0; arg < process.inputs.size(); arg++) {
    const std::string& input        = process.inputs[arg];
    const std::string& input_as_bin = input + ".bin";

    std::cout << "process_xyz_to_bin" << std::endl;
    process_xyz_to_bin(process, input, input_as_bin);
    std::cout << "process_bin_get_resolution" << std::endl;
    process_bin_get_resolution(process, input_as_bin);
    std::cout << "process_bin_get_aabb" << std::endl;

    process_bin_get_aabb(process, input_as_bin);

    print_aabb("AABB Dataset", process.aabb);

    if (process.point_count == 0) {
      std::cout << "Was not able to read any data." << std::endl;
      return -1;
    }

    if (process.generate_normal) {
      process_normalmap(process, input, input_as_bin);
    } else {
      process_xvb(process, input, input_as_bin);
    }

    if (process.export_bbox) {
      process_export_bbox(process, input, input_as_bin);
    }

  }
  return 0;
}
