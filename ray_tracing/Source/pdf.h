#pragma once

#include "vec.h"
#include "Object.h"
#include "utils.h"

#ifndef _PDF_
#define _PDF_

class PDF {
public:
	virtual float value(const vec4 &dir)const = 0;
	virtual vec4 generate() const = 0;
};


class HitablePDF : public PDF {
public:
	vec4 ori;
	HitableObj *ptr;

	HitablePDF(HitableObj* p, const vec4& o);

	virtual float value(const vec4& dir) const;
	virtual vec4 generate() const;

};

class CosinePDF :public PDF {
public:
	onb uvw;
	CosinePDF(const vec4& normal);

	virtual float value(const vec4& dir)const;
	virtual vec4 generate()const;
};




class MixPDF :public PDF {
public:
	PDF* p[2];
	MixPDF(PDF* p0, PDF* p1);

	virtual float value(const vec4& dir)const;
	virtual vec4 generate() const;
};

#endif // !_PDF_


