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
#include <s_hull_pro/s_hull_pro.h>

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
// #include "utils.h"

struct AABB {

  AABB()
      : mMin(std::numeric_limits<real>::max()),
        mMax(-std::numeric_limits<real>::max()) 
	{}

	AABB(const vec3& mi, const vec3& ma) : mMin(mi), mMax(ma) {

	}

	void offset(const vec3& o) {
		mMin += o;
		mMax += o;
	}

  void add(const vec3 &p) {
    mMin = glm::min(p, mMin);
    mMax = glm::max(p, mMax);
  }

  vec3 len() const { return mMax - mMin; }

  bool contains(const vec3 &p) {
    return p.x > mMin.x && p.x <= mMax.x && p.y > mMin.y && p.y <= mMax.y &&
           p.z > mMin.z && p.z <= mMax.z;
  }

	bool contains(const AABB& bb) {
		if (bb.mMin.x > mMax.x || bb.mMax.x < mMin.x)
			return false;
		else if (bb.mMin.y > mMax.y || bb.mMax.y < mMin.y)
			return false;
		else if (bb.mMin.z > mMax.z || bb.mMax.z < mMin.z)
			return false;

		return true;
	}
  vec3 mMin;
  vec3 mMax;
};

typedef std::vector<vec3> tPoints;
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

bool pointSortPredicate(const Shx& a, const Shx& b)
{
	if (a.r < b.r)
		return true;
	else if (a.r > b.r)
		return false;
	else if (a.c < b.c)
		return true;
	else
		return false;
};

void generate_grid_s_hull(const tPoints &points, struct triangulateio *out) {
	size_t size = points.size();
	std::vector<Shx> pts;
	
	pts.resize(size);

	
	
	for (size_t i = 0; i < size; ++i) {	
		const auto &p = points[i];		
		pts[i] = Shx(p.x, p.y);
		pts[i].id = i;
	}

	// std::sort(pts.begin(), pts.end(), pointSortPredicate);
	


	memset(out, 0, sizeof(*out));
	//triangulate("zQ", &in, out, NULL);
	std::vector<Triad> triads;

	s_hull_pro(pts, triads);

	out->numberofpoints = pts.size();
	out->numberofpointattributes = 1;
	out->pointlist = (REAL *)malloc(out->numberofpoints * 2 * sizeof(REAL));
	out->pointattributelist = (REAL *)malloc(out->numberofpoints * out->numberofpointattributes * sizeof(REAL));
	
	for (size_t i = 0; i < out->numberofpoints; ++i) {
		auto indexDst = i * 2;
		auto id = pts[i].id;
		const auto &p = points[i];
		out->pointlist[indexDst + 0] = p.x;
		out->pointlist[indexDst + 1] = p.y;
		out->pointattributelist[i]   = p.z;
	}
	out->numberofcorners = 3;
	out->numberoftriangles = triads.size();
	out->trianglelist = (int*)malloc(sizeof(int) * out->numberofcorners * out->numberoftriangles);
	for (size_t i = 0; i < out->numberoftriangles; ++i) {
		auto index = i * out->numberofcorners;
		auto tri = triads[i];
		
		out->trianglelist[index + 0] = tri.a;
		out->trianglelist[index + 1] = tri.b;
		out->trianglelist[index + 2] = tri.c;
	}

}

void generate_grid(const tPoints &points, struct triangulateio *out) {
	generate_grid_s_hull(points, out);
	
	/*struct triangulateio in;
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
	triangulate("zQ", &in, out, NULL);
	
	free(in.pointattributelist);
	free(in.pointlist);*/
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

void generate_test_grid() {
	std::ofstream output("test.xyz");

	vec3 min(-100, -100, 100);
	vec3 max(100, 100, 200);
	const int grid_width = 8;
	const int grid_height = 8;
	vec3 delta = (max - min) / vec3(grid_width, grid_height, 1);

	for (int y = 0; y < grid_height; ++y) {
		for (int x = 0; x < grid_width; ++x) {
			vec3 pt = min + delta * vec3(x, y, 1);
			pt.z = min.z + ((1 + glm::cos(((real)(x + y) / (real)grid_width) * glm::two_pi<real>()) ) * (real) 0.5) * delta.z;
			output << pt.x << " " << pt.y << " " << pt.z << std::endl;
		}
	}
}

void output_grid(int* triangles, tPoints& points, tNodeGrid& nodes) {
	std::ofstream output("debug.obj");
	int index = 0;
	output << "g " << std::endl;
	for (auto& v : points) {
		output << "v " << v.x << " " << v.y << " " << v.z << std::endl;
	}

	for (auto& node : nodes.mNodes) {
		if (!node.indices.empty()) {
			output << "g node" << index++ << std::endl;

			for (auto& i : node.indices) {
				auto i1 = triangles[i] + 1;
				auto i2 = triangles[i + 1] + 1;
				auto i3 = triangles[i + 2] + 1;

				output << "f " << i1 << " " << i2 << " " << i3 << std::endl;
			}
		}
	}
}

template <typename OUT, typename T>
bool write_binary(std::ofstream& out, T& value) {
	OUT v = (OUT)value;
	return out.write((char*)&v, sizeof(v)).good();
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

	cmdparser.parse_check(ac, av);
	setlocale(LC_ALL, "C");

	auto &inputs = cmdparser.rest();

	
	for (size_t arg = 0; arg < inputs.size(); arg++) {
		tPoints points_input;
		tPoints points_result;
		AABB aabb;


		std::cout << "** processing " << inputs[arg] << std::endl;

		{
			Timing _("1. loads points, this can take a bit of time ( and memory ) depending on the size of the datasets");
			std::ifstream input(inputs[arg]);
			vec3 p;
			for (std::string line; std::getline(input, line); ) {
				size_t comment = line.find("#");
				if (comment > 0) {
					if (comment != std::string::npos) {
						line = line.substr(0, comment);
					}
					if (std::sscanf(line.c_str(), "%f %f %f", &p.x, &p.y, &p.z) == 3) {
						points_input.push_back(p);
						aabb.add(p);
					}
				}
			}
		}
		vec3 len = aabb.len();


		{
			Timing _("2. removing duplicate points");

			const int cNodeGridWidth = 256;
			const int cNodeGridHeight = 256;

			points_result.reserve(points_input.size());
			auto size = points_input.size();
			tNodeGrid grid(cNodeGridWidth, cNodeGridHeight, aabb);
			for (auto &p : points_input) {
				tNode &node = grid.nodeAt(p);
				if (!node.find_point(points_result, p)) {
					node.indices.push_back((uint32_t) points_result.size());
					points_result.push_back(p);
				}
			}
		}


		points_input.clear();
		struct triangulateio mid;

		{
			Timing _("3. triangulation of the dataset");
			generate_grid(points_result, &mid);
		}


		points_result.clear();

		int pt_count_x = cmdparser.get<int>("simplify-split");
		int pt_count_y = pt_count_x;

		if (cmdparser.get<bool>("simplify-split-use-ratio")) {
			real ratio = len.x / len.y;
			pt_count_x = (int)glm::max<real>((pt_count_y * ratio), 4);
			pt_count_x = glm::min<int>(pt_count_x, pt_count_y);
		}

		std::cout << "** simplification grid size " << pt_count_x << " * " << pt_count_y << std::endl;

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


			// output_grid(mid.trianglelist, points_input, grid);


			vec3 direction(0, 0, -1);
			vec3 start(0, 0, glm::abs(aabb.mMax.z) * 2);
			real delta_x = len.x / (real)(pt_count_x + 1);
			real delta_y = len.y / (real)(pt_count_y + 1);
			bool generate_missing = cmdparser.get<bool>("simplify-generate-missing");
			real previous = cmdparser.get<bool>("simplify-use-min") ? glm::min(aabb.mMin.z, aabb.mMax.z) : glm::max(aabb.mMin.z, aabb.mMax.z);
			bool use_previous = !(cmdparser.get<bool>("simplify-use-min") | cmdparser.get<bool>("simplify-use-max"));

			for (int y = 0; y < pt_count_y; ++y) {
				for (int x = 0; x < pt_count_x; ++x) {
					start.x = aabb.mMin.x + x * delta_x;
					start.y = aabb.mMin.y + y * delta_y;

					tNode &node = grid.nodeAt(start);
					bool generated = false;
					for (auto index : node.indices) {
						auto i1 = mid.trianglelist[index + 0];
						auto i2 = mid.trianglelist[index + 1];
						auto i3 = mid.trianglelist[index + 2];

						vec3 res;
						if (glm::intersectRayTriangle(start, direction, points_input[i1], points_input[i2], points_input[i3], res)
							|| glm::intersectRayTriangle(start, direction, points_input[i1], points_input[i3], points_input[i2], res)) {
							points_result.push_back(start + direction * res.z);
							generated = true;
							break;
						}
					}

					if (generate_missing && !generated) {
						if (use_previous) {
							size_t index = y * pt_count_x + x;
							if (y > 0) {
								previous = points_result[index - pt_count_x].z;
							}
							else if (x > 0) {
								previous = points_result[index - 1].z;
							}
						}
						points_result.push_back(vec3(start.x, start.y, previous));
					}
				}
			}
		}

		free(mid.pointlist);
		free(mid.pointattributelist);
		free(mid.pointmarkerlist);
		free(mid.trianglelist);
	
		{
			Timing _("5. tesselation of the simplified dataset");
			generate_grid(points_result, &mid);
		}

		{
			Timing _("6. generating output XVB");

			std::string		output_file_name = inputs[arg] + ".xvb";
			std::ofstream output(output_file_name.c_str(), std::ios_base::binary | std::ios_base::trunc);
				

			write_binary<int>(output, pt_count_x);
			write_binary<int>(output, pt_count_y);

			float min = aabb.mMin.x;
			write_binary<float>(output, aabb.mMin.x);
			write_binary<float>(output, aabb.mMin.y);
			write_binary<float>(output, aabb.mMin.z);

			write_binary<float>(output, aabb.mMax.x);
			write_binary<float>(output, aabb.mMax.y);
			write_binary<float>(output, aabb.mMax.z);

			for (auto& v : points_result) {
				float z = -v.z;
				write_binary<float>(output, z);
			}
		}

		{
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
