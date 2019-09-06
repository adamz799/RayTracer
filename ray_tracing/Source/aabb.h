#pragma once


#include <stdlib.h>
#include <math.h>
//#include "vec.h"
#include "ray.h"
#include "utils.h"

#ifndef _AABB_
#define _AABB_

class aabb {
public:
	vec4 _min, _max;

	aabb() {}
	aabb(const vec4 &a, const vec4 &b) { _min = a; _max = b; }
#ifdef SIMD
	bool hit(const ray &r, float tmin, float tmax)const {
		float flag[4], t0[4], t1[4];
		__m128 temp_min = _mm_sub_ps(_mm_loadu_ps(_min.e), _mm_loadu_ps(r.origin().e));
		__m128 temp_max = _mm_sub_ps(_mm_loadu_ps(_max.e), _mm_loadu_ps(r.origin().e));

		__m128 invD = _mm_rcp_ps(_mm_load_ps(r.direction().e));

		temp_min = _mm_mul_ps(temp_min, invD);
		temp_max = _mm_mul_ps(temp_max, invD);

		__m128 t = _mm_cmpgt_ps(_mm_setzero_ps(), invD);
		
		_mm_storeu_ps(flag, t);
		_mm_storeu_ps(t0, temp_min);
		_mm_storeu_ps(t1, temp_max);

		
		for (int i = 0; i < 3; ++i) {
			if (flag[i]) {
				std::swap(t0[i], t1[i]);
			}
			tmin = max(t0[i], tmin);
			tmax = min(t1[i], tmax);
		}
		return (tmax > tmin);
}
#else
	bool hit(const ray &r, float tmin, float tmax)const {
		vec4 temp_min = _min - r.origin();
		vec4 temp_max = _max - r.origin();
		for (int a = 0; a < 3; ++a) {
			float invD = 1.0f / r.direction()[a];
			float t0 = temp_min[a] * invD;
			float t1 = temp_max[a] * invD;
			if (invD < 0.0f) {
				std::swap(t0, t1);
			}
			tmin = max(t0, tmin);
			tmax = min(t1, tmax);
			if (tmax <= tmin) { return false; }
		}
		return true;
}
#endif // SIMD

	
};

aabb surrounding_box(aabb box0, aabb box1) {
	vec4 smaller(min(box0._min.x(), box1._min.x()),
		min(box0._min.y(), box1._min.y()),
		min(box0._min.z(), box1._min.z())
	);
	vec4 bigger(max(box0._max.x(), box1._max.x()),
		max(box0._max.y(), box1._max.y()),
		max(box0._max.z(), box1._max.z())
	);
	return aabb(smaller, bigger);
}

#endif // !_AABB_
