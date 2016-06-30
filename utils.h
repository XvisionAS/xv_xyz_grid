#pragma once 

namespace glm {
	template <typename T>
	GLM_FUNC_QUALIFIER void operator=
		(
			tvec3<T> & v1,
			tvec3<T> const & v2
			)
	{
		v1 = v2;
	}
}
