#include "common.h"
#include <ostream>

struct process_t {
  std::vector<std::string> inputs;

  vec2 resolution;
};

void process_bin_get_resolution(process_t &process,
                                const std::string &input_as_bin) {
  const std::string input_as_bin_res = input_as_bin + ".res";

  // if (file_exists(input_as_bin_res)) {
  //   std::ifstream input(input_as_bin_res, std::ios::binary | std::ios::in);
  //   input.read((char*)&process.resolution, sizeof(process.resolution));
  // } else
  {
    zstr::ifstream i1(input_as_bin, std::ios::binary | std::ios::in);
    vec3 p1, p2;
    process.resolution.x = REAL_MAX;
    process.resolution.y = REAL_MAX;

    const real epsilon = std::numeric_limits<real>::epsilon();

    // this is real slow, and we hope that the data is well formatted
    // ( so that the point are already sorted, basically, which is the case
    // for 99.5% of our files other wise, the user might have to input the
    // correct resolution
    while (i1.read((char *)&p1, sizeof(p1)) &&
           i1.read((char *)&p2, sizeof(p2))) {
      real dx = glm::abs(p1.x - p2.x);
      real dy = glm::abs(p1.y - p2.y);

      if (dx > epsilon) {
        if (process.resolution.x == REAL_MAX) {
          process.resolution.x = dx;
        } else {
          process.resolution.x = glm::max(process.resolution.x, dx);
        }
      }
      if (dy > epsilon) {
        if (process.resolution.y == REAL_MAX) {
          process.resolution.y = dx;
        } else {
          process.resolution.y = glm::max(process.resolution.x, dx);
        }
      }
    }

    std::ofstream output(input_as_bin_res, std::ios::binary | std::ios::out);
    output.write((char *)&process.resolution, sizeof(process.resolution));
  }
}

int main(int ac, char **av) {
  process_t process;

  setlocale(LC_ALL, "C");
  cmdline::parser cmdparser;

  cmdparser.parse_check(ac, av);

  process.inputs = cmdparser.rest();
  for (size_t arg = 0; arg < process.inputs.size(); arg++) {
    const std::string &input = process.inputs[arg];
    const std::string &input_as_bin = input + ".bin";

    std::cout << "process_bin_get_resolution " << input_as_bin << std::endl;
    
    process_bin_get_resolution(process, input_as_bin);

    std::cout << "-- computed resolution (x=" << process.resolution.x;
    std::cout << ",y=" << process.resolution.y << ")" << std::endl;
  }
}
