#pragma once

#include "vec.h"

#ifndef _RAY_
#define _RAY_

class ray {
public:
	vec4 ori, dir;
	float _time;
	ray() {}
	ray(const vec4 &a, const vec4 &b, float ti = 0.0f) :ori(a), dir(b), _time(ti) { }
	vec4 origin() const { return ori; }
	vec4 direction() const { return dir; }
	vec4 hit_point(const float t) { return ori + t * dir; }
	float time() const { return _time; }
};



#endif // !_RAY_