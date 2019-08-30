#pragma once

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "vec.h"
#include <windows.h>
#include <WinUser.h>
//#include <random>


LRESULT CALLBACK WinProc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{

	switch (msg)
	{
	case WM_CLOSE: {
		PostQuitMessage(0);
		break;
	}
				   //case WM_KEYDOWN: {
				   //	switch (wparam) {
				   //	case 'S': {
				   //		vec3 dir = camera.N;
				   //		camera.origin -= dir * 0.1f;
				   //		camera.look_at -= dir * 0.1f;
				   //		camera.update();
				   //		break;
				   //	}
				   //	case 'W': {
				   //		vec3 dir = camera.N;
				   //		camera.origin += dir * 0.1f;
				   //		camera.look_at += dir * 0.1f;
				   //		camera.update();
				   //		break;
				   //	}
				   //	case 'A': {
				   //		vec3 dir = camera.U;
				   //		camera.origin -= dir * 0.1f;
				   //		camera.look_at -= dir * 0.1f;
				   //		camera.update();
				   //		break;
				   //	}
				   //	case 'D': {
				   //		vec3 dir = camera.U;
				   //		camera.origin += dir * 0.1f;
				   //		camera.look_at += dir * 0.1f;
				   //		camera.update();
				   //		break;
				   //	}
				   //	default: {
				   //		break;
				   //	}
				   //	}
				   //	break;
				   //}
				   //case WM_MOUSEMOVE: {
				   //	//printf("move: %d\n", leftMouseDown);
				   //	if (leftMouseDown) {
				   //		int x = (int)LOWORD(lparam);
				   //		int y = (int)HIWORD(lparam);
				   //		if (lastX != -1 && lastY != -1) {
				   //			if (x == lastX && y == lastY) { break; }
				   //			float sensitivity = 0.03f;
				   //			float xoffset = (x - lastX) * sensitivity;
				   //			float yoffset = (lastY - y) * sensitivity;
				   //			float rotateAmount = sqrt(xoffset*xoffset + yoffset * yoffset);
				   //			vec4 rotateAxis = vec4(-yoffset / rotateAmount, xoffset / rotateAmount, 0.f, 0.f).unit();
				   //			/*printf("offset: %f,%f", xoffset / rotateAmount, yoffset / rotateAmount);
				   //			printf("axis: %f,%f,%f\n", rotateAxis[0], rotateAxis[1], rotateAxis[2]);*/
				   //			Matrix4x4 view = calViewMatrix(camera);
				   //			rotateAxis = view.inverse().mul_vec(rotateAxis);
				   //			Matrix4x4 rotateMat = identityMatrix4x4();
				   //			vec3 axis = rotateAxis.toVec3().unit();
				   //			matRotate(rotateMat, rotateAmount, axis);
				   //			//printf("axis: %f,%f,%f\n", axis[0], axis[1], axis[2]);
				   //			vec4 cameraPos((camera.origin - camera.look_at), 1.f);
				   //			cameraPos = rotateMat.mul_vec(cameraPos);
				   //			camera.origin = cameraPos.toVec3() + camera.look_at;
				   //			camera.update();
				   //		}
				   //		lastX = x;
				   //		lastY = y;
				   //	}
				   //	break;
				   //}
				   //case WM_LBUTTONDOWN: {
				   //	leftMouseDown = true;
				   //	//printf("down: %d\n", leftMouseDown);
				   //	break;
				   //}
				   //case WM_LBUTTONUP: {
				   //	leftMouseDown = false;
				   //	lastX = -1;
				   //	lastY = -1;
				   //	//printf("up: %d\n", leftMouseDown);
				   //	break;
				   //}
				   //}
				   
	}
	return DefWindowProc(hwnd, msg, wparam, lparam);
}


struct Color {
	float R, G, B;
};


struct coord {
	int x, y;
};


//std::minstd_rand gen(std::random_device{}());
//std::uniform_real_distribution<float> dist(0, 1);

inline float randf()
{
	return static_cast<float>(rand() / float(RAND_MAX));
	//return dist(gen);
}

inline vec4 random_in_unit_disk()
{
	vec4 p;
	do {
		p = 2.f*vec4(randf(), randf(), 0.f, 0.f) - vec4(1.f, 1.f, 0.f, 0.f);
	} while (dot(p, p) > 1.f);
	return p;
}

inline vec4 random_in_unit_sphere()
{
	vec4 p;
	do {
		p = vec4(randf(), randf(), randf())*2.f - vec4(1.f, 1.f, 1.f);
	} while (p.square_length() > 1.0f);
	return p.unit();
}

inline vec4 random_cosine_direction() {
	float r1 = randf();
	float r2 = randf();
	float srqt_r2 = 2 * sqrt(r2);
	float z = sqrt(1 - r2);
	float phi = 2 * M_PI*r1;
	float x = cos(phi) * srqt_r2;
	float y = sin(phi) * srqt_r2;
	return vec4(x, y, z);
}

#ifndef RAY
#define RAY
class ray{
public:
	vec4 ori, dir;
	float _time;
	ray(){}
	ray(const vec4 &a, const vec4 &b, float ti = 0.0f):ori(a), dir(b), _time(ti) { }
	vec4 origin() const { return ori; }
	vec4 direction() const { return dir; }
	vec4 hit_point(const float t) { return ori + t * dir; }
	float time() const { return _time; }
};



#endif // !Ray

class Material;
struct hit_record {
	float t;
	vec4 p;
	vec4 normal;
	Material *mat_ptr;
	float u, v;
};

inline vec4 reflect(const vec4 &in, const vec4 &n)
{
	vec4 normal = n.unit();
	return in - 2.f*dot(in, normal)*normal;
}

bool refract(const vec4 &v, const vec4 &n, float ni_over_nt, vec4 & refracted)
{
	vec4 unit_v = v.toVec3().unit();
	float dt = dot(unit_v, n);
	float dis = 1.f - ni_over_nt * ni_over_nt*(1.f - dt * dt);
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
	return r0 + (1.f - r0)*pow((1.f - cosine), 5);
}

inline vec4 random_to_sphere(float radius, float distance_quared) {
	float r1 = randf(), r2 = randf();
	float z = 1 + r2 * (sqrt(1.f - radius * radius / distance_quared) - 1.f);
	float phi = 2 * M_PI*r1;
	float temp = sqrt(1.f - z * z);
	float x = cos(phi)*temp;
	float y = sin(phi)*temp;
	return vec4(x, y, z);
}


class onb {
public:
	vec3 axis[3];

	onb(){}
	vec3 u() const { return axis[0]; }
	vec3 v() const { return axis[1]; }
	vec3 w() const { return axis[2]; }
	vec3 local(const vec3 &a)const{
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

inline vec4 de_nan(const vec4 &c) {
	vec4 temp(c);
	for (int i = 0; i < 3; ++i) {
		if (!(temp[i] == temp[i])) { temp[i] = 0; }
	}
	return temp;
}