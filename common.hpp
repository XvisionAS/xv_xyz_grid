#ifndef __COMMON_H__
#define __COMMON_H__

#include <glm/glm.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtx/intersect.hpp>

#if defined(SINGLE)
typedef float real;
typedef glm::vec3 vec3;
#	define REAL float
# define EPSILON FLT_EPSILON
#else
typedef double real;
typedef glm::dvec3 vec3;
#	define REAL double
# define EPSILON DBL_EPSILON
#endif


#endif
