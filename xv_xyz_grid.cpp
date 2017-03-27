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

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#ifdef _WIN32
#include <io.h>
#define access    _access_s
#else
#include <unistd.h>
#endif

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
	std::vector<std::string>	inputs;
	int							simplifySplit;
	bool						simplifySplitUseRatio;
	bool						alwaysNegate;
	AABB						aabb;
	int							point_count;

	std::vector<real>	bitmap;
	int								bitmap_width;
	int								bitmap_height;

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


	cmdparser.parse_check(ac, av);

	process.inputs = cmdparser.rest();
	process.simplifySplit = cmdparser.get<int>("simplify-split");
	process.simplifySplitUseRatio = cmdparser.get<bool>("simplify-split-use-ratio");
	process.alwaysNegate = cmdparser.get<bool>("always-negate");
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
				vec3	p;
				if (std::sscanf(line.c_str(), SCANF_FORMAT, &p.x, &p.y, &p.z) == 3) {
					output.write((char*)&p, sizeof(p));
					process.point_count++;
					process.aabb.add(p);
				}
			}
		}
	}
}

void process_bin_get_aabb(process_t& process, const std::string& inputAsBIN) {
	std::ifstream input(inputAsBIN, std::ios::binary | std::ios::in);
	vec3 p;
	while (input.read((char*)&p, sizeof(p))) {
		process.point_count++;
		process.aabb.add(p);
	}
}

void process_bin_to_bitmap(process_t& process, const std::string& inputAsBIN) {
	process.bitmap_width = process.simplifySplit;
	process.bitmap_height= process.simplifySplit;

	vec3 len = process.aabb.len();

	if (process.simplifySplitUseRatio) {
		real bitmap_ratio = len.x / len.y;
		if (bitmap_ratio > 0) {
			process.bitmap_width = (int)glm::max<real>((process.bitmap_width * bitmap_ratio), 4);
			process.bitmap_height= (int)glm::max<real>((process.bitmap_height * (1 / bitmap_ratio)), 4);
		}
	}

	process.bitmap.resize(process.bitmap_width * process.bitmap_height, 0);

	std::vector<real>		count;
	int									width = process.bitmap_width;
	int									height = process.bitmap_height;
	vec3								ratio = vec3(process.bitmap_width - 1, process.bitmap_height - 1, 1) / len;
	std::ifstream				input(inputAsBIN, std::ios::binary | std::ios::in);
	vec3								p;

	if (glm::isinf(ratio.x)) {
		ratio.x = 0;
	}
	if (glm::isinf(-ratio.y)) {
		ratio.y = 0;
	}
	if (glm::isinf(-ratio.z)) {
		ratio.z = 0;
	}
	
	count.resize(process.bitmap.size(), 0);

	while (input.read((char*)&p, sizeof(p))) {
		real x = (p.x - process.aabb.mMin.x) * ratio.x;
		real y = (p.y - process.aabb.mMin.y) * ratio.y;

		int		px = (int)x;
		int   py = (int)y;

		int   nx = px + 1;
		int   ny = py + 1;

		real nfx = (x - (real)px);
		real nfy = (y - (real)py);

		real pfx = 1 - nfx;
		real pfy = 1 - nfy;

		real s = pfx * pfy;

		count[px + py * width] += s;
		process.bitmap[px + py * width] += (p.z * s);

		s = nfx * pfy;
		if (nx < width) {
			count[nx + py * width] += s;
			process.bitmap[nx + py * width] += (p.z * s);
		}
		
		s = pfx * nfy;
		if (ny < height) {
			count[px + ny * width] += s;
			process.bitmap[px + ny * width] += (p.z * s);
		}

		s = nfx * nfy;
		if ( (ny < height) && (nx < width) ) {
			count[nx + ny * width] += s;
			process.bitmap[nx + ny * width] += (p.z * s);
		}
	}

	size_t  size = process.bitmap.size();
	int			empty = 0;
	for (size_t i = 0; i < size; ++i) {
		real c = count[i];
		if (c > 0) {
			process.bitmap[i] /= count[i];
		}
		else {
			process.bitmap[i] = REAL_MAX;
		}
	}
}

void process_fill_bitmap(process_t& process) {
	int  index = 0;
	for (int y = 0; y < process.bitmap_height; ++y) {
		for (int x = 0; x < process.bitmap_width; ++x) {
			if (process.bitmap[index] == REAL_MAX) {
				real values = 0.0f;
				int	 count = 0;

				for (int yy = -1; yy != 1; ++yy) {
					for (int xx = -1; xx != 1; ++xx) {
						int cy = yy + y;
						int cx = xx + x;
						if (cy >= 0 && cx >= 0 && cy < process.bitmap_height && cx < process.bitmap_width) {
							int index_sub = cy * process.bitmap_width + cx;
							real actual = process.bitmap[index_sub];
							if (actual != REAL_MAX) {
								values += actual;
								count++;
							}
						}
					}
				}
				if (count > 0) {
					process.bitmap[index] = values / count;
				}
			}
			index++;
		}
	}
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

	std::ofstream output(outputFile.c_str(), std::ios_base::binary | std::ios_base::trunc);

	write_binary<int>(output, process.bitmap_width);
	write_binary<int>(output, process.bitmap_height);

	write_binary<float>(output, process.aabb.mMin.x);
	write_binary<float>(output, process.aabb.mMin.y);
	write_binary<float>(output, process.aabb.mMin.z);

	write_binary<float>(output, process.aabb.mMax.x);
	write_binary<float>(output, process.aabb.mMax.y);
	write_binary<float>(output, process.aabb.mMax.z);

	for (auto& v : process.bitmap) {
		real d = (v == REAL_MAX) ? (process.aabb.mMin.z * 2.0f) : v;
		write_binary<float>(output, d);
	}
}

int main(int ac, char **av) {
	setlocale(LC_ALL, "C");

	process_t process;
	parse_cmd_line(ac, av, process);
	for (size_t arg = 0; arg < process.inputs.size(); arg++) {
		const std::string& input = process.inputs[arg];
		const std::string& inputAsBIN = input + ".bin";
		std::cout << "xv_xyz_grid: process_xyz_to_bin" << std::endl;
		process_xyz_to_bin(process, input, inputAsBIN);
		std::cout << "xv_xyz_grid: process_bin_get_aabb" << std::endl;
		process_bin_get_aabb(process, inputAsBIN);

		if (process.point_count == 0) {
			std::cout << "Was not able to read any data." << std::endl;
			return -1;
		}

		/*vec3 len = process.aabb.len();
		if (glm::length(len) < glm::epsilon<real>()) {
			std::cout << "Dataset is empty." << std::endl;
			return -2;
		}*/

		std::cout << "xv_xyz_grid: process_bin_to_bitmap" << std::endl;
		process_bin_to_bitmap(process, inputAsBIN);
		std::cout << "xv_xyz_grid: process_bitmap_negate" << std::endl;
		process_bitmap_negate(process);
		std::cout << "xv_xyz_grid: process_bitmap_to_png" << std::endl;
		process_bitmap_to_png(process, input + ".before.png");
		// process_fill_bitmap(process);
		std::cout << "xv_xyz_grid: process_bitmap_to_xvb" << std::endl;
		process_bitmap_to_xvb(process, input + ".xvb");
	}
  return 0;
}
