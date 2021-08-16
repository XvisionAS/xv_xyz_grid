#ifndef __AABB_H__
#define __AABB_H__

struct AABB {

  AABB()
      : min(std::numeric_limits<real>::max()),
        max(-std::numeric_limits<real>::max()) {}

  AABB(const vec3 &mi, const vec3 &ma) : min(mi), max(ma) {}
  void reset() {
    min = vec3(std::numeric_limits<real>::max());
    max = -min;
  }

  void offset(const vec3 &o) {
    min += o;
    max += o;
  }

  void add(const vec3 &p) {
    min = glm::min(p, min);
    max = glm::max(p, max);
  }

  vec3 len() const { return max - min; }

  bool contains(const vec3 &p) {
    return p.x > min.x && p.x <= max.x && p.y > min.y && p.y <= max.y &&
           p.z > min.z && p.z <= max.z;
  }

  bool contains2d(const vec3 &p) {
    return p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y;
  }

  bool contains(const AABB &bb) {
    if (bb.min.x > max.x || bb.max.x < min.x)
      return false;
    else if (bb.min.y > max.y || bb.max.y < min.y)
      return false;
    else if (bb.min.z > max.z || bb.max.z < min.z)
      return false;

    return true;
  }
  vec3 min;
  vec3 max;
};

#endif