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

#define VOID void

#include <triangle/triangle.h>

#include "AABB.h"


typedef std::vector<vec3> tPoints;


void generate_grid(const tPoints &points, struct triangulateio *out) {
	struct triangulateio in;
	memset(&in, 0, sizeof(in));
	in.numberofpoints = (int) points.size();
	in.numberofpointattributes = 1;

	in.pointlist = (REAL *)malloc(in.numberofpoints * 2 * sizeof(REAL));
	in.pointattributelist = (REAL *)malloc(in.numberofpoints * in.numberofpointattributes * sizeof(REAL));

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
	triangulate((char*)"zQ", &in, out, NULL);
	
	free(in.pointattributelist);
	free(in.pointlist);
}

struct Timing {
	Timing(const char* description) {
		mStart = std::chrono::steady_clock::now();
		std::cout << description;
		std::cout.flush();

	}

	~Timing() {
		auto diff = std::chrono::steady_clock::now() - mStart;
		auto time = std::chrono::duration_cast<std::chrono::milliseconds>(diff);                                                      
		std::cout << " [took " << time.count() << " milliseconds]" << std::endl;
	}

	std::chrono::steady_clock::time_point mStart;
};

typedef std::vector<uint32_t> tIndices;

struct tNode {
	tIndices indices;
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
		: mCoordStart(aabb.mMin), mCoordSize(aabb.len()) {
		mNodes.resize((width + 1) * (height + 1));

		mGridWidth = width;
		mGridHeight = height;
	}

	tNode &nodeAt(const vec3 &p) {
		int x = (int)glm::floor(((p.x - mCoordStart.x) / mCoordSize.x) *
			(real)(mGridWidth));
		int y = (int)glm::floor(((p.y - mCoordStart.y) / mCoordSize.y) *
			(real)(mGridHeight));
		return mNodes[x + y * mGridWidth];
	}


	std::vector<tNode> mNodes;

	std::int_fast32_t mGridWidth;
	std::int_fast32_t mGridHeight;

	vec3 mCoordStart;
	vec3 mCoordSize;
};

template <typename OUT, typename T>
bool write_binary(std::ofstream& out, T& value) {
	OUT v = (OUT)value;
	return out.write((char*)&v, sizeof(v)).good();
}


bool file_exists(const std::string &Filename)
{
	return access(Filename.c_str(), 0) == 0;
}

int main(int ac, char **av) {
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
		false);
	cmdparser.add<bool>("simplify-generate-missing",
		'c',
		"Generate missing depth value in dataset, using previous value",
		false,
		true);

	cmdparser.add<bool>("simplify-use-min",
		'd',
		"Missing value will use minimum depth value",
		false,
		true);

	cmdparser.add<bool>("simplify-use-max",
		'e',
		"Missing value will use maximum depth value",
		false,
		false);

	cmdparser.add<bool>("always-negate",
		'f',
		"Depth value will always be set to negative value",
		false,
		true);

	cmdparser.add<bool>("generate-obj",
		'g',
		"Generate an OBJ file of the generated mesh",
		false,
		false);

	cmdparser.parse_check(ac, av);
	setlocale(LC_ALL, "C");

	auto&	inputs = cmdparser.rest();
	int		split = cmdparser.get<int>("simplify-split");

	
	for (size_t arg = 0; arg < inputs.size(); arg++) {
		tPoints points_input;
		tPoints points_result;
		AABB aabb;

		size_t point_count = 0;
		if (!file_exists(inputs[arg] + ".bin"))
		{
			Timing _("1. loads points, this can take a bit of time ( and memory ) depending on the size of the datasets, compute AABB");
			std::ifstream input(inputs[arg]);
			std::ofstream output(inputs[arg] + ".bin", std::ios::binary | std::ios::out);

			vec3 p;
			for (std::string line; std::getline(input, line); ) {
				size_t comment = line.find("#");
				if (comment > 0) {
					if (comment != std::string::npos) {
						line = line.substr(0, comment);
					}
					if (std::sscanf(line.c_str(), SCANF_FORMAT, &p.x, &p.y, &p.z) == 3) {
						point_count++;
						output.write((char*) &p, sizeof(p));
						aabb.add(p);
					}
				}
			}
		} else {
			Timing _("1. loads points ( from cached binary ), this can take a bit of time ( and memory ) depending on the size of the datasets, compute AABB");
			std::ifstream input(inputs[arg] + ".bin", std::ios::binary | std::ios::in);
			vec3 p;
			while (input.read((char*)&p, sizeof(p))) {
				aabb.add(p);
				point_count++;
			}
		}

		if (point_count == 0) {			
			std::cout << "Was not able to read any data." << std::endl;
			return -1;
		}

		vec3 len = aabb.len();
		
		if (glm::length(len) < glm::epsilon<real>()) {
			std::cout << "Dataset is empty." << std::endl;
			return -2;
		}

		int pt_count_x = cmdparser.get<int>("simplify-split");
		int pt_count_y = pt_count_x;

		if (cmdparser.get<bool>("simplify-split-use-ratio")) {
			real ratio = len.x / len.y;
			pt_count_x = (int)glm::max<real>((pt_count_y * ratio), 4);
			pt_count_x = glm::min<int>(pt_count_x, pt_count_y);
		}

		std::vector<real>	value;
		

		value.resize(pt_count_x * pt_count_y, 0);
		
		{
			Timing _("2. converts point to a bitmap");
			
			std::vector<real>	count;
			count.resize(value.size(), 0);
			vec3 ratio = vec3(pt_count_x - 1, pt_count_y - 1, 1) / len;		
			std::ifstream input(inputs[arg] + ".bin", std::ios::binary |  std::ios::in);
			vec3 p;
			
			while (input.read((char*)&p, sizeof(p))) {
				real x = (p.x - aabb.mMin.x) * ratio.x;
				real y = (p.y - aabb.mMin.y) * ratio.y;

				int		px = (int)x;
				int   py = (int)y;
				
				int   nx = px + 1;
				int   ny = py + 1;

				real nfx = (x - (real) px);
				real nfy = (y - (real) py);

				real pfx = 1 - nfx;
				real pfy = 1 - nfy;

				real s = pfx * pfy;
				if (nx < pt_count_x) {
					count[px + py * pt_count_x] += s;
					value[px + py * pt_count_x] += (p.z * s);
				}
				
				s = nfx * pfy;
				if (nx < pt_count_x) {
					count[nx + py * pt_count_x] += s;
					value[nx + py * pt_count_x] += (p.z * s);
				}
				
				s = nfx * nfy;
				if (ny < pt_count_y) {
					count[nx + ny * pt_count_x] += s;
					value[nx + ny * pt_count_x] += (p.z * s);
				}
				s = pfx * nfy;
				if (ny < pt_count_y) {
					count[px + ny * pt_count_x] += s;
					value[px + ny * pt_count_x] += (p.z * s);
				}
			}
			size_t  size = value.size();
			int			empty = 0;
			for (size_t i = 0; i < size; ++i) {
				real c = count[i];
				if (c > 0) {
					value[i] /= count[i];
				} else {
					value[i] = REAL_MAX;
				}
			}			
		}
		

		struct triangulateio mid;
		{
			Timing _("3. generating grid from bitmap points");
			size_t size = value.size();
			points_input.clear();
			points_input.reserve(size);

			vec3 ratio = len / vec3(pt_count_x - 1, pt_count_y - 1, 1);
			for (int i = 0; i < size; ++i) {
				real v = value[i];
				int y = i / pt_count_x;
				int x = i - (y * pt_count_x);

				if (v != REAL_MAX) {
					points_input.push_back(vec3(x * ratio.x + aabb.mMin.x, y * ratio.y + aabb.mMin.y, value[i]));
				}
				else {
					points_input.push_back(vec3(x * ratio.x + aabb.mMin.x, y * ratio.y + aabb.mMin.y, aabb.mMin.z));
				}
			}
			
			generate_grid(points_input, &mid);
			points_input.clear();					
		}
		
		{
			const int cNodeGridWidth = 16;
			const int cNodeGridHeight = 16;

			Timing _("4. generation of the simplified dataset");
			tNodeGrid grid(cNodeGridWidth, cNodeGridHeight, aabb);

			for (auto i = 0; i < mid.numberofpoints; ++i) {
				auto index = i * 2;
				points_input.push_back(vec3(mid.pointlist[index + 0], mid.pointlist[index + 1], mid.pointattributelist[i]));
			}

			for (auto i = 0; i < mid.numberoftriangles; i++) {
				auto index = i * mid.numberofcorners;

				auto i1 = mid.trianglelist[index + 0];
				auto i2 = mid.trianglelist[index + 1];
				auto i3 = mid.trianglelist[index + 2];

				AABB tri_bbox;
				tri_bbox.add(points_input[i1]);
				tri_bbox.add(points_input[i2]); 
				tri_bbox.add(points_input[i3]);

				vec3 grid_delta = grid.mCoordSize / (vec3(grid.mGridWidth, grid.mGridHeight, 1));
				for (int y = 0; y < grid.mGridHeight; ++y) {
					for (int x = 0; x < grid.mGridWidth; ++x) {
						AABB aabb(grid.mCoordStart, grid.mCoordStart + grid_delta);
						aabb.offset(grid_delta * vec3(x, y, 0));
						if (aabb.contains(tri_bbox)) {
							auto node_index = x + y * grid.mGridWidth;
							grid.mNodes[node_index].indices.push_back(index);
						}
					}
				}
			}
			vec3 direction(0, 0, -1);
			vec3 start(0, 0, glm::abs(aabb.mMax.z) * 2);
			
			real delta_x = len.x / (real)(pt_count_x - 1);
			real delta_y = len.y / (real)(pt_count_y - 1);
			bool generate_missing = cmdparser.get<bool>("simplify-generate-missing");
			real minZ = glm::abs(aabb.mMin.z);
			real maxZ = glm::abs(aabb.mMax.z);

			real missing = cmdparser.get<bool>("simplify-use-max") ? glm::max(minZ, maxZ) : glm::min(minZ, maxZ);
			
			points_result.reserve(pt_count_x * pt_count_y);
			for (int y = 0; y < pt_count_y; ++y) {
				for (int x = 0; x < pt_count_x; ++x) {
					int index = y * pt_count_x + x;
					start.x = aabb.mMin.x + x * delta_x;
					start.y = aabb.mMin.y + y * delta_y;
					
					if (value[index] == REAL_MAX) {
						tNode &node = grid.nodeAt(start);
						bool generated = false;
						for (auto index : node.indices) {
							auto i1 = mid.trianglelist[index + 0];
							auto i2 = mid.trianglelist[index + 1];
							auto i3 = mid.trianglelist[index + 2];

							vec3 res;
							if (glm::intersectRayTriangle(start, direction, points_input[i1], points_input[i2], points_input[i3], res)) {
								points_result.push_back(start + direction * res.z);
								generated = true;
								break;
							}
						}

						if (generate_missing && !generated) {
							points_result.push_back(vec3(start.x, start.y, missing));
						}
					}
					else {
						points_result.push_back(vec3(start.x, start.y, value[index]));
					}
				}
			}
		}

		{
			Timing _("6. generating output XVB");

			std::string		output_file_name = inputs[arg] + ".xvb";
			std::ofstream output(output_file_name.c_str(), std::ios_base::binary | std::ios_base::trunc);


			if (cmdparser.get<bool>("always-negate")) {
				aabb.reset();
				for (auto& v : points_result) {
					if (v.z > 0) {
						v.z = -v.z;
					}
					aabb.add(v);
				}
			}
			
			write_binary<int>(output, pt_count_x);
			write_binary<int>(output, pt_count_y);

			write_binary<float>(output, aabb.mMin.x);
			write_binary<float>(output, aabb.mMin.y);
			write_binary<float>(output, aabb.mMin.z);

			write_binary<float>(output, aabb.mMax.x);
			write_binary<float>(output, aabb.mMax.y);
			write_binary<float>(output, aabb.mMax.z);

			for (auto& v : points_result) {
				write_binary<float>(output, v.z);
			}
		}

		if (cmdparser.get<bool>("generate-obj")) {
			Timing _("7. generating output OBJ");
			std::string output_file_name = inputs[arg] + ".obj";
			std::ofstream output(output_file_name.c_str());
			for (auto i = 0; i < mid.numberofpoints; ++i) {
				auto index = i * 2;
				output << std::fixed
					<< "v " << (mid.pointlist[index + 0])
					<< "  " << (mid.pointlist[index + 1])
					<< "  " << (mid.pointattributelist[i]) << std::endl;
			}

			for (auto i = 0; i < mid.numberoftriangles; i++) {
				output << "f ";
				for (auto j = 0; j < mid.numberofcorners; j++) {
					output << mid.trianglelist[i * mid.numberofcorners + j] + 1 << " ";
				}
				output << std::endl;
			}
		}

		free(mid.pointlist);
		free(mid.pointattributelist);
		free(mid.pointmarkerlist);
		free(mid.trianglelist);		
	}	
  return 0;
}
