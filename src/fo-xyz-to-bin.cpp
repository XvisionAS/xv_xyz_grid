/*
 * This convert a XYZ to a binary blob that is way faster to load.
 * This is done because we are processing the same file multiple time over and over
 * and we are doing everything "out of core" ( so not in memory )
 */

#include "common.h"

struct process_t {
  std::vector<std::string>  inputs;

  bool use_NODATA;
  double NODATA;
};

void process_xyz_to_bin(process_t &process, const std::string &inputFileName, const std::string &outputFileName) {
  std::ifstream input(inputFileName);
  zstr::ofstream output(outputFileName, std::ios::binary | std::ios::out);

  for (std::string line; std::getline(input, line);) {
    size_t comment = line.find("#");
    if (comment > 0) {
      if (comment != std::string::npos) {
        line = line.substr(0, comment);
      }
      vec3 p;
      const real epsilon = std::numeric_limits<real>::epsilon();
      if (std::sscanf(line.c_str(), SCANF_FORMAT, &p.x, &p.y, &p.z) == 3) {

        // maybe glm::epsilonNotEqual(process.NODATA, p.z, epsilon) ?
        if (!process.use_NODATA || process.NODATA != p.z) {
          output.write((char *)&p, sizeof(p));
        }
      }
    }
  }
}

int main(int ac, char **av) {
  process_t process;

  setlocale(LC_ALL, "C");
  cmdline::parser cmdparser;

  cmdparser.add<double>("no-data",
    'd',
    "Define NODATA value",
    false,
    0
  );

  cmdparser.parse_check(ac, av);

  process.use_NODATA = cmdparser.exist("no-data");
  if (process.use_NODATA) {
    process.NODATA  = cmdparser.get<double>("no-data");
    std::cout << "Using NODATA " << process.NODATA << std::endl;
  }

  process.inputs                   = cmdparser.rest();

  for (size_t arg = 0; arg < process.inputs.size(); arg++) {
    const std::string& input        = process.inputs[arg];
    const std::string& input_as_bin = input + ".bin";

    std::cout << "process_xyz_to_bin " << input <<  " -> " << input_as_bin << std::endl;
    process_xyz_to_bin(process, input, input_as_bin);
  }
}
