#ifndef __AABB_H__
#define __AABB_H__

struct tAABB {
  tAABB()
      : mMin(std::numeric_limits<real>::max()),
        mMax(-std::numeric_limits<real>::max())
	{}
	tAABB(const vec3& mi, const vec3& ma) : mMin(mi), mMax(ma) {

	}
	void offset(const vec3& o) {
		mMin += o;
		mMax += o;
	}
  void add(const vec3 &p) {
    mMin = glm::min(p, mMin);
    mMax = glm::max(p, mMax);
  }
  vec3 size() const {
		return mMax - mMin;
	}
  bool contains(const vec3 &p) {
    return p.x > mMin.x && p.x <= mMax.x
				&& p.y > mMin.y && p.y <= mMax.y
				&& p.z > mMin.z && p.z <= mMax.z;
  }
	bool contains(const tAABB& bb) {
		return !(bb.mMin.x > mMax.x || bb.mMax.x < mMin.x
		|| bb.mMin.y > mMax.y || bb.mMax.y < mMin.y
		|| bb.mMin.z > mMax.z || bb.mMax.z < mMin.z);
	}
  vec3 mMin;
  vec3 mMax;
};

#endif
