#pragma once

#include "utils.h"
#include "Material.h"
#include "aabb.h"

#ifndef _OBJECT_
#define _OBJECT_

aabb surrounding_box(aabb box0, aabb box1);

class HitableObj {
public:
	virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const = 0;
	virtual bool bounding_box(float t0, float t1, aabb &box) = 0;
	virtual float pdf_value(const vec4 &o, const vec4 &v)const = 0;
	virtual vec4 random(const vec4 &o)const = 0;
};



class HitableList : public HitableObj{
public:
	HitableObj **list;
	int list_size;
	HitableList() {};
	HitableList(HitableObj **l, int n) :list(l), list_size(n) {};

	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec)const;
	virtual bool bounding_box(float t0, float t1, aabb& box);	
	virtual float pdf_value(const vec4& o, const vec4& v)const;
	virtual vec4 random(const vec4& o)const;
};


int box_x_compare(const void* a, const void* b);

int box_y_compare(const void* a, const void* b);

int box_z_compare(const void* a, const void* b);


class BVHNode : public HitableObj {
public:
	HitableObj * left, *right;
	aabb box;

	BVHNode(): left(nullptr), right(nullptr){}
	BVHNode(HitableObj** l, int n, float time0, float time1);

	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& b);
	virtual float pdf_value(const vec4& o, const vec4& v)const;
	virtual vec4 random(const vec4& o)const;
};


class Sphere : public HitableObj {
public:
	vec4 center;
	float radius;
	Material *mat_ptr = NULL;
	aabb aabb_box;
	bool dirty;
	Sphere(){}
	Sphere(vec4 c, float r) : center(c), radius(r),dirty(true) { };
	Sphere(vec4 c, float r, Material *m_p) : center(c), radius(r), mat_ptr(m_p), dirty(true) { };
	
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box);
	void get_sphere_uv(const vec4& p, float& u, float& v) const;
	virtual float pdf_value(const vec4& o, const vec4& v)const;
	virtual vec4 random(const vec4& o)const;
};

class MovingSphere :public HitableObj {
public:
	vec4 center0, center1;
	float radius;
	float time0, time1;
	Material *mat_ptr;
	aabb aabb_box;
	bool dirty;

	MovingSphere() {};
	MovingSphere(const vec4& cen0, const vec4& cen1, float t0, float t1, float r, Material* m);
	
	vec4 center(float time) const;

	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec);
	virtual bool bounding_box(float t0, float t1, aabb& box);
};


class Parallelogram : public HitableObj {
public:
	vec4 ori, u, v, normal;
	Material *mp;
	bool dirty;
	aabb aabb_box;

	Parallelogram() {}
	Parallelogram(const vec4& _ori, const vec4& _u, const vec4& _v, Material* mat);

	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box);
	virtual float pdf_value(const vec4& o, const vec4& v)const;
	virtual vec4 random(const vec4& o)const;
};

class ConstantMedium : public HitableObj {
public:
	HitableObj *boundry;
	float density;
	Material *phase_function;

	ConstantMedium(HitableObj* h, float d, Texture* p) ;

	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box);
	virtual float pdf_value(const vec4& o, const vec4& v)const;
	virtual vec4 random(const vec4& o)const;
};


#endif // !_OBJECT_