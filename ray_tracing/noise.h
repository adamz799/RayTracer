#pragma once

#ifndef NOISE
#define NOISE

#include <math.h>
#include "vec.h"


inline float hermite(float f) {
	return f * f*(3.f - 2.f * f);
}

static vec4* perlin_generate() {
	vec4 *p = new vec4[256];
	for (int i = 0; i < 256; ++i) {
		p[i] = vec4(2.f *vec3(randf(), randf(), randf()) - 1.f).unit();
	}
	return p;
}

void permute(int *p, int n) {
	for (int i = n - 1; i > 0; --i) {
		int target = int(randf()*(i + 1));
		int tmp = p[i];
		p[i] = p[target];
		p[target] = tmp;
	}
	return;
}

static int * perlin_generate_perm() {
	int *p = new int[256];
	for (int i = 0; i < 256; ++i) {
		p[i] = i;
	}
	permute(p, 256);
	return p;
}

inline float trilinear_interp(float c[2][2][2], float u, float v, float w) {
	float accum = 0;
	for (int i = 0; i < 2; ++i) {
		for (int k = 0; k < 2; ++k) {
			for (int j = 0; j < 2; ++j) {
				accum += (i*u + (1.f - i)*(1.f - u))*
					(j*v + (1.f - j)*(1.f - v))*
					(k*w + (1.f - k)*(1.f - w))*
					c[i][j][k];
			}
		}
	}
	return accum;
}

inline float perlin_interp(const vec4 c[2][2][2], float u, float v, float w) {
	float uu = hermite(u);
	float vv = hermite(v);
	float ww = hermite(w);
	float accum = 0;
	for (int i = 0; i < 2; ++i) {
		for (int k = 0; k < 2; ++k) {
			for (int j = 0; j < 2; ++j) {
				vec4 weight_v(u - i, v - j, w - k);
				accum += (i*uu + (1.f - i)*(1.f - uu))
					*(j*vv + (1.f - j)*(1.f - vv))
					*(k*ww + (1.f - k)*(1.f - ww))
					*dot(c[i][j][k], weight_v);
			}
		}
	}
	return accum;
}


class Perlin {
public:
	static vec4 *randVec4;
	static int *perm_x, *perm_y, *perm_z;
	Perlin() {
		
	}
	float noise(const vec4 &p)const {
		float u = p.x() - floor(p.x());
		float v = p.y() - floor(p.y());
		float w = p.z() - floor(p.z());
		int i = floor(p.x());
		int j = floor(p.y());
		int k = floor(p.z());
		vec4 c[2][2][2];
		for (int di = 0; di < 2; ++di) {
			for (int dk = 0; dk < 2; ++dk) {
				for (int dj = 0; dj < 2; ++dj) {
					c[di][dj][dk] = randVec4[
						perm_x[(i+di)&255] 
						^ perm_y[(j+dj)&255] 
						^ perm_z[(k+dk)&255]];
				}
			}
		}
		return perlin_interp(c,u,v,w);
	}
	

	float turb(const vec4 &p, int depth = 7)const {
		float accum = 0;
		vec4 temp_p = p;
		float weight = 1.0f;
		for (int i = 0; i < depth; ++i) {
			accum += weight * noise(temp_p);
			weight *= 0.5;
			temp_p *= 2;
		}
		return fabs(accum);
	}

};

vec4 * Perlin::randVec4 = perlin_generate();
int * Perlin::perm_x = perlin_generate_perm();
int * Perlin::perm_y = perlin_generate_perm();
int * Perlin::perm_z = perlin_generate_perm();

#endif // !NOISE

