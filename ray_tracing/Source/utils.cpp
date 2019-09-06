#pragma once

#include "utils.h"


LRESULT CALLBACK WinProc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{

	switch (msg)
	{
	case WM_CLOSE: {
		PostQuitMessage(0);
		break;
	}
	}
	return DefWindowProc(hwnd, msg, wparam, lparam);
}


//std::minstd_rand gen(std::random_device{}());
//std::uniform_real_distribution<float> dist(0, 1);

vec4 random_in_unit_disk()
{
	vec4 p;
	do {
		p = 2.f * vec4(randf(), randf(), 0.f, 0.f) - vec4(1.f, 1.f, 0.f, 0.f);
	} while (dot(p, p) > 1.f);
	return p;
}

vec4 random_in_unit_sphere()
{
	vec4 p;
	do {
		p = vec4(randf(), randf(), randf()) * 2.f - vec4(1.f, 1.f, 1.f);
	} while (p.square_length() > 1.0f);
	return p.unit();
}

vec4 random_cosine_direction() {
	float r1 = randf();
	float r2 = randf();
	float srqt_r2 = 2.f * sqrt(r2);
	float z = sqrt(1 - r2);
	float phi = 2.f * M_PI * r1;
	float x = cos(phi) * srqt_r2;
	float y = sin(phi) * srqt_r2;
	return vec4(x, y, z);
}


vec4 reflect(const vec4& in, const vec4& n)
{
	vec4 normal = n.unit();
	return in - 2.f * dot(in, normal) * normal;
}

bool refract(const vec4& v, const vec4& n, float ni_over_nt, vec4& refracted)
{
	vec4 unit_v = v.toVec3().unit();
	float dt = dot(unit_v, n);
	float dis = 1.f - ni_over_nt * ni_over_nt * (1.f - dt * dt);
	if (dis > 0.f)
	{
		refracted = ni_over_nt * (unit_v - n * dt) - n * sqrt(dis);
		return true;
	}
	else
	{
		return false;
	}
}

float schlick(float cosine, float ref_idx) {
	float r0 = (1.f - ref_idx) / (1.f + ref_idx);
	r0 *= r0;
	return r0 + (1.f - r0) * pow((1.f - cosine), 5);
}

vec4 random_to_sphere(float radius, float distance_quared) {
	float r1 = randf(), r2 = randf();
	float z = 1 + r2 * (sqrt(1.f - radius * radius / distance_quared) - 1.f);
	float phi = 2 * M_PI * r1;
	float temp = sqrt(1.f - z * z);
	float x = cos(phi) * temp;
	float y = sin(phi) * temp;
	return vec4(x, y, z);
}


vec4 de_nan(const vec4& c) {
	vec4 temp(c);
	for (int i = 0; i < 3; ++i) {
		if (!(temp[i] == temp[i])) { temp[i] = 0; }
	}
	return temp;
}

