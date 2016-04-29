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

#if defined(SINGLE)
typedef float real;
typedef glm::vec3 vec3;
#define REAL float
#else
typedef double real;
typedef glm::dvec3 vec3;
#define REAL double
#endif

#define ANSI_DECLARATORS
#define VOID void

extern "C" {
#include <triangle/triangle.h>
}

struct AABB {

  AABB()
      : min(std::numeric_limits<real>::max()),
        max(-std::numeric_limits<real>::max()) {}

  void add(const vec3 &p) {
    min = glm::min(p, min);
    max = glm::max(p, max);
  }

  vec3 len() const { return max - min; }

  bool contains(const vec3 &p) {
    return p.x > min.x && p.x <= max.x && p.y > min.y && p.y <= max.y &&
           p.z > min.z && p.z <= max.z;
  }
  vec3 min;
  vec3 max;
};

typedef std::vector<vec3> tPoints;

void generate_grid(const tPoints &points, struct triangulateio *out) {
  struct triangulateio in;
  memset(&in, 0, sizeof(in));
  in.numberofpoints = points.size();
  in.numberofpointattributes = 1;

  in.pointlist = (REAL *)malloc(in.numberofpoints * 2 * sizeof(REAL));
  in.pointattributelist = (REAL *)malloc(
      in.numberofpoints * in.numberofpointattributes * sizeof(REAL));
  for (auto i = 0; i < in.numberofpoints; ++i) {
    auto indexDst = i * 2;
    const auto &p = points[i];
    in.pointlist[indexDst + 0] = p.x;
    in.pointlist[indexDst + 1] = p.y;
    // not sure if it should be put somewhere else.
    in.pointattributelist[i] = p.z;
  }

  in.numberofsegments = 0;
  in.numberofholes = 0;
  in.numberofregions = 0;
  memset(out, 0, sizeof(*out));
  triangulate("zQ", &in, out, NULL);
  free(in.pointattributelist);
  free(in.pointlist);
}

struct tNode {
  std::vector<uint32_t> indices;
  bool find_point(const tPoints &points, const vec3 &p) {
    for (auto index : indices) {
      const vec3 &o = points[index];
      if (p == o) {
        return true;
      }
    }
    return false;
  }
  bool find_point_epsilon(const tPoints &points, const vec3 &p, real epsilon) {
    for (auto index : indices) {
      const vec3 &o = points[index];
      if (glm::epsilonEqual(p.x, o.x, epsilon) &&
          glm::epsilonEqual(p.y, o.y, epsilon) &&
          glm::epsilonEqual(p.z, o.z, epsilon)) {
        return true;
      }
    }
    return false;
  }
};

struct tNodeGrid {
  tNodeGrid(std::int_fast32_t width, std::int_fast32_t height, const AABB &aabb)
      : mCoordStart(aabb.min), mCoordSize(aabb.len()) {
    mNodes.resize(width * height);

    mGridWidth = width;
    mGridHeight = height;
  }

  tNode &nodeAt(const vec3 &p) {
    int x = (int)glm::floor((p.x - mCoordStart.x) / mCoordSize.x *
                            (real)(mGridWidth - 1));
    int y = (int)glm::floor((p.y - mCoordStart.y) / mCoordSize.y *
                            (real)(mGridHeight - 1));
    return mNodes[x + y * mGridWidth];
  }

  std::vector<tNode> mNodes;

  std::int_fast32_t mGridWidth;
  std::int_fast32_t mGridHeight;

  vec3 mCoordStart;
  vec3 mCoordSize;
};

#define START(description)                                                     \
std::cout                                                                      \
      << description;                                                          \
start = std::chrono::steady_clock::now();

#define END()                                                                  \
  {                                                                            \
auto time = std::chrono::duration_cast<std::chrono::milliseconds>(             \
                std::chrono::steady_clock::now() - start)                      \
                .count();                                                      \
std::cout                                                                      \
        << " [took " << time << " milliseconds]" << std::endl;                   \
}

int main(int ac, char **av) {
  cmdline::parser cmdparser;

  cmdparser.add<bool>("simplify", 'i',
                      "Simplify output by generating a grid if given size",
                      false, true);
  cmdparser.add<int>("simplify-split", 'c', "Number of split used to generate "
                                            "the simplification grid (value x "
                                            "value)",
                     false, 256);
  cmdparser.add<bool>("simplify-generate-missing", 'g',
                      "When generating simplification grid, point that are not "
                      "part of the dataset will be generated using the lowest "
                      "value in the dataset",
                      false, true);
  cmdparser.add<bool>("output-obj", 'o', "Output an OBJ file", false, false);
  cmdparser.add<float>("scale", 's', "Multiply all points by this value", false,
                       1);

  cmdparser.parse_check(ac, av);
  setlocale(LC_ALL, "C");

  auto &inputs = cmdparser.rest();

  for (size_t arg = 0; arg < inputs.size(); arg++) {
    tPoints points_input;

    std::cout << "** processing " << inputs[arg] << std::endl;

    std::chrono::steady_clock::time_point start;

    START("1. loads points, this can take a bit of time ( and memory ) "
          "depending on the size of the datasets");
    start = std::chrono::steady_clock::now();
    std::ifstream input(inputs[arg]);
    real x, y, z;
    AABB aabb;
    while (input >> x >> y >> z) {
      vec3 p(x, y, z);

      points_input.push_back(p);
      aabb.add(p);
    }
    END();

    START("2. removing duplicate points");

    tPoints points_result;

    points_result.reserve(points_input.size());

    auto size = points_input.size();
    vec3 len = aabb.len();

    const int cNodeGridWidth = 250;
    const int cNodeGridHeight = 250;

    {
      tNodeGrid grid(cNodeGridWidth, cNodeGridHeight, aabb);
      for (auto &p : points_input) {
        tNode &node = grid.nodeAt(p);
        if (!node.find_point(points_result, p)) {
          node.indices.push_back(points_result.size());
          points_result.push_back(p);
        }
      }
    }
    END();

    // let's create a 100x100 grid from the aabb.
    points_input.clear();

    START("3. triangulation of the dataset");

    struct triangulateio mid;
    generate_grid(points_result, &mid);
    END();

    points_result.clear();
    if (cmdparser.get<bool>("simplify")) {
      START("4. generation of the simplified dataset");
      {
        tNodeGrid grid(cNodeGridWidth, cNodeGridHeight, aabb);

        for (auto i = 0; i < mid.numberoftriangles; i++) {
          for (auto j = 0; j < mid.numberofcorners; j++) {
            int index = mid.trianglelist[i * mid.numberofcorners + j];
            real *ptrxy = mid.pointlist + index * 2;

            tNode &node = grid.nodeAt(vec3(ptrxy[0], ptrxy[1], 0));
            node.indices.push_back(i);
          }
        }

        for (auto i = 0; i < mid.numberofpoints; ++i) {
          auto index = i * 2;
          points_input.push_back(vec3(mid.pointlist[index + 0],
                                      mid.pointlist[index + 1],
                                      mid.pointattributelist[i]));
        }

        int pt_count = cmdparser.get<int>("simplify-split");

        vec3 direction(0, 0, -1);
        vec3 start(0, 0, aabb.max.z * 2);
        real delta_x = len.x / (real)pt_count;
        real delta_y = len.y / (real)pt_count;
        bool generate_missing =
            cmdparser.get<bool>("simplify-generate-missing");
        for (int x = 0; x < (pt_count + 1); ++x) {
          for (int y = 0; y < (pt_count + 1); ++y) {
            start.x = aabb.min.x + x * delta_x;
            start.y = aabb.min.y + y * delta_y;

            tNode &node = grid.nodeAt(start);
            bool generated = false;
            for (auto index : node.indices) {
              index *= mid.numberofcorners;

              int i1 = mid.trianglelist[index + 0];
              int i2 = mid.trianglelist[index + 1];
              int i3 = mid.trianglelist[index + 2];

              vec3 res;
              if (glm::intersectRayTriangle(start, direction, points_input[i1],
                                            points_input[i2], points_input[i3],
                                            res)) {
                points_result.push_back(start + direction * res.z);
                generated = true;
                break;
              }
            }

            if (generate_missing && !generated) {
              points_result.push_back(vec3(start.x, start.y, aabb.min.z));
            }
          }
        }
      }
      END();

      free(mid.pointlist);
      free(mid.pointattributelist);
      free(mid.pointmarkerlist);
      free(mid.trianglelist);
      START("5. tesselation of the simplified dataset");
      generate_grid(points_result, &mid);
      END();
    }

    START("6. generating output OBJ");
		std::string outputFileName = inputs[arg] + ".obj";

    std::ofstream output(outputFileName.c_str());
    auto scale = cmdparser.get<float>("scale");
    std::cout << ", scale value " << scale;
    for (auto i = 0; i < mid.numberofpoints; ++i) {
      auto index = i * 2;
      output << std::fixed << "v " << (mid.pointlist[index + 0] * scale) << " "
             << (mid.pointlist[index + 1] * scale) << " "
             << (mid.pointattributelist[i] * scale) << std::endl;
    }

    for (auto i = 0; i < mid.numberoftriangles; i++) {
      output << "f ";
      for (auto j = 0; j < mid.numberofcorners; j++) {
        output << mid.trianglelist[i * mid.numberofcorners + j] + 1 << " ";
      }
      output << std::endl;
    }

    END();

    free(mid.pointlist);
    free(mid.pointattributelist);
    free(mid.pointmarkerlist);
    free(mid.trianglelist);
  }
  return 0;
}
