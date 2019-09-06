#pragma once

#include <windows.h>
#include <WinUser.h>
#include "vec.h"
#include <stdlib.h>
#include <random>

#ifndef _UTILS_
#define _UTILS_

LRESULT CALLBACK WinProc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam);

inline float randf()
{
	return static_cast<float>(rand() / float(RAND_MAX));
	//return dist(gen);
}

vec4 random_in_unit_disk();

vec4 random_in_unit_sphere();

vec4 random_cosine_direction();

vec4 reflect(const vec4& in, const vec4& n);

bool refract(const vec4& v, const vec4& n, float ni_over_nt, vec4& refracted);

float schlick(float cosine, float ref_idx);

vec4 random_to_sphere(float radius, float distance_quared);

vec4 de_nan(const vec4& c);


#endif // !_UTILS_
