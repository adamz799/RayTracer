#pragma once

#ifndef _PDF
#define _PDF

#include "Object.h"

class PDF {
public:
	virtual float value(const vec4 &dir)const = 0;
	virtual vec4 generate() const = 0;
};


//class HitablePDF :public PDF {
//public:
//	vec4 ori;
//	HitableObj *ptr;
//
//	HitablePDF() {}
//
//	HitablePDF(HitableObj *p, const vec4 &o) {
//		ptr = p;
//		ori = o;
//	}
//
//	virtual float value(const vec4 &dir) const {
//		return ptr->pdf_value(ori, dir);
//	}
//	virtual vec4 generate() const {
//		return ptr->random(ori);
//	}
//
//};

class CosinePDF :public PDF {
public:
	onb uvw;
	CosinePDF(const vec4 &normal) { uvw.build(normal); }
	virtual float value(const vec4 &dir)const  {
		float cosine = dot(dir.unit(), vec4(uvw.w()));
		if (cosine > 0) { return cosine / M_PI; }
		else { return 0; }
	}

	virtual vec4 generate()const {
		return uvw.local(random_cosine_direction().toVec3());
	}
};




class MixPDF :public PDF {
public:
	PDF *p[2];
	MixPDF(PDF* p0, PDF *p1) { p[0] = p0, p[1] = p1; }

	virtual float value(const vec4 &dir)const {
		return 0.5f*p[0]->value(dir) + 0.5*p[1]->value(dir);
	}

	virtual vec4 generate() const {
		if (randf() < 0.5) {
			return p[0]->generate();
		}
		else {
			return p[1]->generate();
		}
	}
};

#endif // !_PDF

struct scatter_record {
	ray specular_ray;
	bool is_specular;
	vec4 attenuation;
	PDF *pdf_ptr;
};
