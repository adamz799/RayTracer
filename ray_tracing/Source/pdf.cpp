#pragma once

#include "pdf.h"
#include "vec.h"
#include "Object.h"
#include "utils.h"
#include "head.h"


HitablePDF::HitablePDF(HitableObj* p, const vec4& o):ptr(p),ori(o){}

float HitablePDF::value(const vec4& dir) const
{
	return ptr->pdf_value(ori, dir);
}

vec4 HitablePDF::generate() const
{
	return ptr->random(ori);
}



CosinePDF::CosinePDF(const vec4& normal) { uvw.build(normal); }

float CosinePDF::value(const vec4& dir)const {
	float cosine = dot(dir.unit(), vec4(uvw.w()));
	if (cosine > 0) { return cosine / M_PI; }
	else { return 0; }
}

vec4 CosinePDF::generate()const {
	return uvw.local(random_cosine_direction().toVec3());
}




MixPDF::MixPDF(PDF* p0, PDF* p1) { p[0] = p0, p[1] = p1; }

float MixPDF::value(const vec4& dir)const {
	return 0.5f * p[0]->value(dir) + 0.5f * p[1]->value(dir);
}

vec4 MixPDF::generate() const {
	if (randf() < 0.5) {
		return p[0]->generate();
	}
	else {
		return p[1]->generate();
	}
}


