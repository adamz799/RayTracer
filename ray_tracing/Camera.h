#pragma once

#include "vec.h"
#include "ray.h"
#include "utils.h"

#ifndef _CAMERA_
#define _CAMERA_


class Camera {
public:
	vec3 origin;
	vec3 look_at;
	vec3 screen_center;
	vec3 horizontal;
	vec3 vertical;
	vec3 U;
	vec3 V;
	vec3 W;
	vec3 up;

	Camera(){}
	Camera(vec3 ori, vec3 look_at_, vec3 vup,  float v_fov, float aspect){//aspect = width/height;
		origin = ori;
		look_at = look_at_;
		up = vup;
		W = (look_at - origin).unit();
		U = cross(W, up).unit();
		V = cross(U, W);
		
		float theta = v_fov * M_PI / 180.;
		float half_height = tan(theta / 2.);
		float half_width = aspect * half_height;
		screen_center = look_at_;
		horizontal = half_width * U;
		vertical = half_height * V;
	}

	virtual ray get_ray(float x, float y)
	{
		return ray(vec4(origin), vec4(screen_center + x * horizontal + y * vertical - origin));
	}
};

class FocusCamera : public Camera {
public:
	float lens_radius;
	FocusCamera(vec3 ori, vec3 look_at_, vec3 vup, float v_fov, float aspect, float aperture, float focus_dist) {
		origin = ori;
		look_at = look_at_;
		up = vup;
		lens_radius = aperture / 2.;
		W = (look_at - origin).unit();
		U = cross(W, up).unit();
		V = cross(U, W);

		float theta = v_fov * M_PI / 180.;
		float half_height = tan(theta / 2.);
		float half_width = aspect * half_height;
		screen_center = origin + focus_dist * W;
		horizontal = half_width * U;
		vertical = half_height * V;	
	}

	virtual ray get_ray(float x, float y)
	{
		vec4 rd = lens_radius * random_in_unit_disk();
		vec4 offset = vec4(U * rd.x() + V * rd.y());
		return ray(vec4(origin) + offset,vec4( screen_center + x * horizontal + y * vertical - origin) - offset);
	}
};

class MovingCamera : public Camera {
public:
	float lens_radius;
	float time0, time1;
	MovingCamera(vec3 ori, vec3 look_at_, vec3 vup, float v_fov, float aspect, float aperture, float focus_dist, float t0, float t1) {
		origin = ori;
		look_at = look_at_;
		up = vup;
		time0 = t0;
		time1 = t1;
		lens_radius = aperture / 2.;
		W = (look_at - origin).unit();
		U = cross(W, up).unit();
		V = cross(U, W);

		float theta = v_fov * M_PI / 180.;
		float half_height = tan(theta / 2.);
		float half_width = aspect * half_height;
		screen_center = origin + focus_dist * W;
		horizontal = half_width * U;
		vertical = half_height * V;
	}

	virtual ray get_ray(float x, float y)
	{
		vec4 rd = lens_radius * random_in_unit_disk();
		vec4 offset = vec4(U * rd.x() + V * rd.y());
		float time = time0 + randf()*(time1 - time0);
		return ray(vec4(origin) + offset, vec4(screen_center + x * horizontal + y * vertical - origin) - offset, time);
	}
};

#endif // !_CAMERA_

