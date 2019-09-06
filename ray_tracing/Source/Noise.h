#pragma once

#include <math.h>
#include "vec.h"
#include "utils.h"

#ifndef _NOISE_
#define _NOISE_

//Functions.

inline float hermite(float f) {
	return f * f * (3.f - 2.f * f);
}

static vec4* perlin_generate() {
	vec4* p = new vec4[256];
	for (int i = 0; i < 256; ++i) {
		p[i] = vec4(2.f * vec3(randf(), randf(), randf()) - 1.f).unit();
	}
	return p;
}

void permute(int* p, int n);

static int* perlin_generate_perm() {
	int* p = new int[256];
	for (int i = 0; i < 256; ++i) {
		p[i] = i;
	}
	permute(p, 256);
	return p;
}

float trilinear_interp(float c[2][2][2], float u, float v, float w);

float perlin_interp(const vec4 c[2][2][2], float u, float v, float w);


class Perlin
{
public:
	static vec4* randVec4;
	static int* perm_x, * perm_y, * perm_z;
	Perlin()=default;
	float noise(const vec4& p)const;
	float turb(const vec4& p, int depth = 7)const;

};



#endif //!_NOISE_
