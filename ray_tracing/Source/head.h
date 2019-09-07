#pragma once

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <windows.h>
#include <WinUser.h>
#include <memory>
#include "vec.h"
#include "ray.h"

#ifndef _HEAD_
#define _HEAD_



struct Color {
	float R, G, B;
};


struct coord {
	int x, y;
};



class Material;
struct hit_record {
	float t;
	vec4 p;
	vec4 normal;
	std::shared_ptr<Material> mat_ptr;
	float u, v;
};

class PDF;
struct scatter_record {
	ray specular_ray;
	bool is_specular;
	vec4 attenuation;
	std::shared_ptr<PDF> pdf_ptr;
	//PDF *pdf_ptr;
};

class onb {
public:
	vec3 axis[3];

	onb() {}
	vec3 u() const { return axis[0]; }
	vec3 v() const { return axis[1]; }
	vec3 w() const { return axis[2]; }
	vec3 local(const vec3 &a)const {
		vec3 result(0.0f);
		for (int i = 0; i < 3; ++i) {
			result += a[i] * axis[i];
		}
		return result;
	}
	void build(const vec4& normal) {
		axis[2] = normal.toVec3().unit();
		vec3 a;
		if (fabs(w()[0]) > 0.9) {
			a = vec3(0.0f, 1.0f, 0.0f);
		}
		else {
			a = vec3(1.0f, 0.0f, 0.0f);
		}
		axis[1] = cross(w(), a).unit();
		axis[0] = cross(w(), v());
	}
};

#endif // !_HEAD_