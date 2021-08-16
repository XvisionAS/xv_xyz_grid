#pragma once

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

typedef double real;
typedef glm::dvec3 vec3;
typedef glm::dvec2 vec2;
#define REAL double
#define SCANF_FORMAT " %lf%*[ \t,;] %lf%*[ \t,;] %lf %*[\n]"
#define REAL_MAX DBL_MAX

#include "AABB.h"


template <typename OUT, typename T>
bool write_binary(std::ofstream& out, T& value) {
  OUT v = (OUT)value;
  return out.write((char*)&v, sizeof(v)).good();
}


inline bool file_exists(const std::string &Filename) {
  return access(Filename.c_str(), 0) == 0;
}