#pragma once

#include "utils.h"
#include "Material.h"
#include "aabb.h"
#include "Object.h"


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


bool HitableList::hit(const ray& r, float t_min, float t_max, hit_record& rec)const
{
	hit_record temp_rec;
	bool hit_anything = false;
	double closest_so_far = t_max;
	for (int i = 0; i < list_size; ++i)
	{
		if (list[i]->hit(r, t_min, closest_so_far, temp_rec))
		{
			hit_anything = true;
			closest_so_far = temp_rec.t;
			rec = temp_rec;
		}
	}
	return hit_anything;
}
bool HitableList::bounding_box(float t0, float t1, aabb& box) {
	if (list_size < 1)return false;
	aabb temp_box;
	bool first_true = list[0]->bounding_box(t0, t1, temp_box);
	if (!first_true) { return false; }
	else { box = temp_box; }
	for (int i = 0; i < list_size; ++i) {
		if (list[0]->bounding_box(t0, t1, temp_box)) {
			box = surrounding_box(box, temp_box);
		}
		else { return false; }
	}
	return true;
}
float HitableList::pdf_value(const vec4& o, const vec4& v)const {
	float weight = 1.f / list_size;
	float sum = 0;
	for (int i = 0; i < list_size; ++i) {
		sum += weight * list[i]->pdf_value(o, v);
	}
	return sum;
}
vec4 HitableList::random(const vec4& o)const {
	int index = int(randf() * list_size);
	return list[index]->random(o);
}


int box_x_compare(const void* a, const void* b) {
	aabb box_left, box_right;
	HitableObj* ah = *(HitableObj * *)a;
	HitableObj* bh = *(HitableObj * *)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right)) {
		std::cerr << "no bounding box in BVHNode constructor(in compare function)\n";
	}
	if (box_left._min.x() - box_right._min.x() < 0.0f) {
		return -1;
	}
	else { return 1; }
}
int box_y_compare(const void* a, const void* b) {
	aabb box_left, box_right;
	HitableObj* ah = *(HitableObj * *)a;
	HitableObj* bh = *(HitableObj * *)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right)) {
		std::cerr << "no bounding box in BVHNode constructor(in compare function)\n";
	}
	if (box_left._min.y() - box_right._min.y() < 0.0f) {
		return -1;
	}
	else { return 1; }
}
int box_z_compare(const void* a, const void* b) {
	aabb box_left, box_right;
	HitableObj* ah = *(HitableObj * *)a;
	HitableObj* bh = *(HitableObj * *)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right)) {
		std::cerr << "no bounding box in BVHNode constructor(in compare function)\n";
	}
	if (box_left._min.z() - box_right._min.z() < 0.0f) {
		return -1;
	}
	else { return 1; }
}


BVHNode::BVHNode(HitableObj** l, int n, float time0, float time1) {
	int axis = int(3 * randf());
	if (axis == 0) {
		qsort(l, n, sizeof(HitableObj*), box_x_compare);
	}
	else if (axis == 1) {
		qsort(l, n, sizeof(HitableObj*), box_y_compare);
	}
	else {
		qsort(l, n, sizeof(HitableObj*), box_z_compare);
	}

	if (n == 1) { left = right = l[0]; }
	else if (n == 2) { left = l[0]; right = l[1]; }
	else {
		left = new BVHNode(l, n / 2, time0, time1);
		right = new BVHNode(l + n / 2, n - n / 2, time0, time1);
	}
	aabb box_left, box_right;
	if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right)) {
		std::cerr << "no bounding box in BHVNode constructor\n";
	}
	box = surrounding_box(box_left, box_right);

}
bool BVHNode::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	if (box.hit(r, t_min, t_max)) {
		hit_record left_record, right_record;
		bool hit_left = left->hit(r, t_min, t_max, left_record);
		bool hit_right = right->hit(r, t_min, t_max, right_record);
		if (hit_left && hit_right) {
			if (left_record.t < right_record.t) { rec = left_record; }
			else { rec = right_record; }
			return true;
		}
		else if (hit_left) {
			rec = left_record;
			return true;
		}
		else if (hit_right) {
			rec = right_record;
			return true;
		}
		else { return false; }
	}
	else { return false; }
}
bool BVHNode::bounding_box(float t0, float t1, aabb& b) {
	b = box;
	return true;
}
float BVHNode::pdf_value(const vec4& o, const vec4& v)const {
	return 0.5f * left->pdf_value(o, v) + 0.5f * right->pdf_value(o, v);
}
vec4 BVHNode::random(const vec4& o)const {
	if (randf() > 0.5f) { return left->random(o); }
	else { return right->random(o); }
}



bool Sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const
{
	vec4 oc = r.origin() - center;
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc, oc) - radius * radius;
	float delta = b * b - a * c;

	if (delta > 0)
	{
		delta = sqrtf(delta);
		float t1 = (-b + delta) / a, t2 = (-b - delta) / a;
		if (t1 > t2) { float temp = t1; t1 = t2; t2 = temp; }
		if (t1<t_max && t1>t_min)
		{
			rec.t = t1;
			rec.p = r.origin() + t1 * r.direction();
			rec.normal = (rec.p - center) / radius;
			get_sphere_uv(rec.normal, rec.u, rec.v);
			rec.mat_ptr = mat_ptr;
			return true;
		}
		t1 = t2;
		if (t1<t_max && t1>t_min)
		{
			rec.t = t1;
			rec.p = r.origin() + t1 * r.direction();
			rec.normal = (rec.p - center) / radius;
			get_sphere_uv(rec.normal, rec.u, rec.v);
			rec.mat_ptr = mat_ptr;
			return true;
		}

	}
	else
	{
		return false;
	}
	return false;
}
bool Sphere::bounding_box(float t0, float t1, aabb& box) {
	if (dirty) {
		aabb_box = aabb(center - vec4(radius), center + vec4(radius));
		dirty = false;
		box = aabb_box;
	}
	else {
		box = aabb_box;
	}
	return true;
}
void Sphere::get_sphere_uv(const vec4& p, float& u, float& v) const {
	float phi = atan2(p.z(), p.x());
	float theta = asin(p.y());
	u = 1 - (phi + M_PI) / (2 * M_PI);
	v = (theta + M_PI / 2) / M_PI;
}
float Sphere::pdf_value(const vec4& o, const vec4& v)const {
	hit_record rec;
	if (this->hit(ray(o, v), 1e-6f, 1e6f, rec)) {
		float cos_theta_max = sqrt(1.f - radius * radius / (center - o).square_length());
		float solid_angle = 2.f * M_PI * (1.f - cos_theta_max);
		return 1 / solid_angle;
	}
	else { return 0; }
}
vec4 Sphere::random(const vec4& o)const {
	vec4  dir = center - o;
	float distance_suqared = dir.square_length();
	onb uvw;
	uvw.build(dir);
	return uvw.local(random_to_sphere(radius, distance_suqared).toVec3());
}



MovingSphere::MovingSphere(const vec4& cen0, const vec4& cen1, float t0, float t1, float r, std::shared_ptr<Material> m)
{
	center0 = cen0;
	center1 = cen1;
	time0 = t0;
	time1 = t1;
	radius = r;
	mat_ptr = m;
	dirty = true;
}
vec4 MovingSphere::center(float time) const {
	return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
}
bool MovingSphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) {
	vec4 temp_center = center(r.time());
	vec4 oc = r.origin() - temp_center;
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc, oc) - radius * radius;
	float delta = b * b - a * c;

	if (delta > 0)
	{
		delta = sqrtf(delta);
		float t1 = (-b + delta) / a, t2 = (-b - delta) / a;
		if (t1 > t2) { float temp = t1; t1 = t2; t2 = temp; }
		if (t1<t_max && t1>t_min)
		{
			rec.t = t1;
			rec.p = r.origin() + t1 * r.direction();
			rec.normal = (rec.p - temp_center) / radius;
			rec.mat_ptr = mat_ptr;
			return true;
		}
		t1 = t2;
		if (t1<t_max && t1>t_min)
		{
			rec.t = t1;
			rec.p = r.origin() + t1 * r.direction();
			rec.normal = (rec.p - temp_center) / radius;
			rec.mat_ptr = mat_ptr;
			return true;
		}
		return false;
	}
	else
	{
		return false;
	}
}
bool MovingSphere::bounding_box(float t0, float t1, aabb& box) {
	if (dirty) {
		aabb_box = aabb(center0 - vec4(radius), center1 + vec4(radius));
		dirty = false;
	}
	box = aabb_box;
	return true;
}



Parallelogram::Parallelogram(const vec4& _ori, const vec4& _u, const vec4& _v, const std::shared_ptr<Material> m) {
	ori = _ori;
	u = _u; v = _v;
	normal = vec4(cross(_u.toVec3(), _v.toVec3()).unit());
	mat_ptr = m;
	dirty = true;
}
bool Parallelogram::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	float dot1 = dot(r.dir, normal);
	if (abs(dot1) < 1e-6) { return false; }
	float t = dot(ori - r.origin(), normal) / dot1;
	if (t<t_min || t>t_max) { return false; }
	vec4 hit_point = r.origin() + t * r.direction();

	float uu = dot(u, u), vv = dot(v, v), uv = dot(u, v);
	float delta = uu * vv - uv * uv;
	if (fabs(delta) < 1e-6) { return false; }

	delta = 1.f / delta;
	vec4 l = hit_point - ori;
	float lu = dot(l, u), lv = dot(l, v);
	float u_val = delta * (vv * lu - uv * lv);
	if (u_val < 0.f || u_val>1.f) { return false; }
	float v_val = delta * (uu * lv - uv * lu);
	if (v_val < 0.f || v_val>1.f) { return false; }

	rec.u = u_val;
	rec.v = v_val;
	rec.t = t;
	rec.mat_ptr = mat_ptr;
	rec.p = hit_point;
	rec.normal = normal;
	return true;
}
bool Parallelogram::bounding_box(float t0, float t1, aabb& box) {
	if (dirty) {
		vec4 temp = ori + u + v;
		vec4 o = ori;
		for (int i = 0; i < 3; ++i) {
			if (o[i] > temp[i]) { std::swap(o[i], temp[i]); }
			if (temp[i] - o[i] < 1e-6) { temp[i] += 1e-5; }
		}
		aabb_box = aabb(o, temp);
		box = aabb_box;
		dirty = false;
	}
	else {
		box = aabb_box;
	}

	return true;
}
float Parallelogram::pdf_value(const vec4& o, const vec4& v)const {
	return 0.0f;
}
vec4 Parallelogram::random(const vec4& o)const {
	return vec4(1.f, 0.f, 0.f);
}


ConstantMedium::ConstantMedium(HitableObj* h, float d, std::shared_ptr<Texture> p) :boundry(h), density(d) {
	phase_function = std::make_shared<Isotropic>(p);
}
bool ConstantMedium::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	//bool db = (randf() < 1e-5);
	bool db = false;
	hit_record rec1, rec2;
	if (boundry->hit(r, -FLT_MAX, FLT_MAX, rec1)) {
		if (boundry->hit(r, rec1.t + 1e-4, FLT_MAX, rec2)) {
			if (db) {
				std::cerr << "t0, t1 " << rec.t << " " << rec2.t << "\n";
			}
			if (rec1.t < t_min) { rec1.t = t_min; }
			if (rec2.t > t_max) { rec2.t = t_max; }
			if (rec1.t >= rec2.t) { return false; }
			if (rec1.t < 0.0f) { rec1.t = 0.0f; }
			float distance_inside_boundary = (rec2.t - rec1.t) * r.dir.length();
			float hit_distance = -(1.f / density) * log(randf());
			if (hit_distance < distance_inside_boundary) {
				if (db) {
					std::cerr << "hit_distance = " << hit_distance << "\n";
				}
				rec.t = rec1.t + hit_distance / r.direction().length();
				if (db) {
					std::cerr << "rec.t = " << rec.t << "\n";
				}
				rec.p = r.origin() + rec.t * r.direction();
				if (db) {
					std::cerr << "rec.p = " << rec.p << "\n";
				}
				rec.normal = vec4(1.0f, 0.0f, 0.0f);
				rec.mat_ptr = phase_function;
				return true;
			}
		}
	}
	return false;
}
bool ConstantMedium::bounding_box(float t0, float t1, aabb& box) {
	return boundry->bounding_box(t0, t1, box);
}
float ConstantMedium::pdf_value(const vec4& o, const vec4& v)const {
	return boundry->pdf_value(o, v);
}
vec4 ConstantMedium::random(const vec4& o)const {
	return boundry->random(o);
}
