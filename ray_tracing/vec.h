#pragma once

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <intrin.h>
#include <iostream>

#define SIMD
constexpr auto M_PI = 3.142592657;


#ifndef _VEC2_
#define _VEC2_

class vec2 {
public:
	float e[2];

	vec2() {}
	vec2(float a, float b) { e[0] = a; e[1] = b; }
	vec2(float a) { e[0] = a; e[1] = a; }
	float x() const { return e[0]; }
	float y() const { return e[1]; }
	float r() const { return e[0]; }
	float g() const { return e[1]; }
	int size() { return 2; }
	float operator[](int i) const { return e[i]; }
	float& operator[](int i) { return e[i]; }
	const vec2& operator+() const { return *this; }
	vec2 operator-() const { return vec2(-e[0], -e[1]); }
	vec2& operator+=(const vec2 &v2) { e[0] += v2[0]; e[1] += v2[1];  return *this; }
	vec2& operator+=(const float &v) { for (int i = 0; i < 2; ++i) { e[i] += v; } return *this; }
	vec2& operator-=(const vec2 &v2) { e[0] -= v2[0]; e[1] -= v2[1];  return *this; }
	vec2& operator-=(const float &v) { for (int i = 0; i < 2; ++i) { e[i] -= v; } return *this; }
	vec2& operator*=(const vec2 &v2) { e[0] *= v2[0]; e[1] *= v2[1];   return *this; }
	vec2& operator*=(const float &v) { for (int i = 0; i < 2; ++i) { e[i] *= v; } return *this; }
	vec2& operator/=(const vec2 &v2) { e[0] /= v2[0]; e[1] /= v2[1];  return *this; }
	vec2& operator/=(const float &v) { for (int i = 0; i < 2; ++i) { e[i] /= v; } return *this; }

	float length() const
	{
		return sqrt(square_length());
	}

	float square_length() const
	{
		return e[0] * e[0] + e[1] * e[1];
	}

	vec2 unit() const
	{
		float l = length();
		return vec2(e[0] / l, e[1] / l);
	}


};


inline std::istream& operator>>(std::istream &is, vec2 &t)
{
	is >> t.e[0] >> t.e[1];
	return is;
}

inline std::ostream& operator<<(std::ostream &os, vec2 &t)
{
	os << t.e[0] << t.e[1];
	return os;
}

inline vec2 operator+(vec2 const &e, float const v)
{
	return vec2(e.x() + v, e.y() + v);
}

inline vec2 operator+(vec2 const &e, vec2 const &f)
{
	return vec2(e.x() + f.x(), e.y() + f.y());
}

inline vec2 operator-(vec2 const &e, float const v)
{
	return vec2(e.x() - v, e.y() - v);
}

inline vec2 operator-(vec2 const &e, vec2 const &f)
{
	return vec2(e.x() - f.x(), e.y() - f.y());
}

inline vec2 operator*(vec2 const &e, float const v)
{
	return vec2(e.x() * v, e.y() * v);
}

inline vec2 operator*(vec2 const &e, vec2 const &f)
{
	return vec2(e.x() * f.x(), e.y() * f.y());
}

inline vec2 operator*(float const v, vec2 const &e)
{
	return vec2(e.x() * v, e.y() * v);
}

inline vec2 operator/(vec2 const &e, float const v)
{
	if (abs(v) > 1e-6) return (vec2(e[0] / v, e[1] / v));
}

inline float dot(vec2 const &e1, vec2 const &e2)
{
	return (e1.e[0] * e2.e[0] + e1.e[1] * e2.e[1]);
}

inline float cross(vec2 const &e1, vec2 const &e2)
{
	return e1[0] * e2[1] - e1[1] * e2[0];
}


#endif // !_VEC2_

#ifndef _VEC3_
#define _VEC3_

class vec3 {
public:
	vec3() {}
	float e[3];
	vec3(float e0) {
		e[0] = e0; e[1] = e0;
		e[2] = e0;
	}
	vec3(float e0, float e1, float e2) {
		e[0] = e0; e[1] = e1;
		e[2] = e2;
	}
	/*vec3(const vec4 &in) {
	for(int i = 0; i < 3; ++i) {
	e[i] = in[i];
	}
	}*/
	float x() const { return e[0]; }
	float y() const { return e[1]; }
	float z() const { return e[2]; }
	float r() const { return x(); }
	float g() const { return y(); }
	float b() const { return z(); }
	int size() { return 3; }
	const vec3& operator+() const { return *this; }
	vec3 operator-() const { return vec3(-x(), -y(), -z()); }

	float operator[](int i) const { return e[i]; }
	float& operator[](int i) { return e[i]; }

	vec3& operator+=(const vec3 &v2) { e[0] += v2.x(); e[1] += v2.y();  e[2] += v2.z(); return *this; }
	vec3& operator+=(const float &v) { for (int i = 0; i < 3; ++i) { e[i] += v; } return *this; }
	vec3& operator-=(const vec3 &v2) { e[0] -= v2.x(); e[1] -= v2.y();  e[2] -= v2.z(); return *this; }
	vec3& operator-=(const float &v) { for (int i = 0; i < 3; ++i) { e[i] -= v; } return *this; }
	vec3& operator*=(const vec3 &v2) { e[0] *= v2.x(); e[1] *= v2.y();  e[2] *= v2.z(); return *this; }
	vec3& operator*=(const float &v) { for (int i = 0; i < 3; ++i) { e[i] *= v; } return *this; }
	vec3& operator/=(const vec3 &v2) { e[0] /= v2.x(); e[1] /= v2.y();  e[2] /= v2.z(); return *this; }
	vec3& operator/=(const float &v) { for (int i = 0; i < 3; ++i) { e[i] /= v; } return *this; }

	float length() const
	{
		return sqrt(square_length());
	}

	float square_length() const
	{
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	vec3 unit() const
	{
		float l = length();
		return vec3(e[0] / l, e[1] / l, e[2] / l);
	}


};

inline std::istream& operator>>(std::istream &is, vec3 &t)
{
	is >> t.e[0] >> t.e[1] >> t.e[2];
	return is;
}

inline std::ostream& operator<<(std::ostream &os, vec3 &t)
{
	os << t.e[0] << t.e[1] << t.e[2];
	return os;
}

inline vec3 operator+(vec3 const &e, float const v)
{
	return vec3(e.x() + v, e.y() + v, e.z() + v);
}

inline vec3 operator+(vec3 const &e, vec3 const &f)
{
	return vec3(e.x() + f.x(), e.y() + f.y(), e.z() + f.z());
}

inline vec3 operator-(vec3 const &e, float const v)
{
	return vec3(e.x() - v, e.y() - v, e.z() - v);
}

inline vec3 operator-(vec3 const &e, vec3 const &f)
{
	return vec3(e.x() - f.x(), e.y() - f.y(), e.z() - f.z());
}

inline vec3 operator*(vec3 const &e, float const v)
{
	return vec3(e.x() * v, e.y() * v, e.z() * v);
}

inline vec3 operator*(vec3 const &e, vec3 const &f)
{
	return vec3(e.x() * f.x(), e.y() * f.y(), e.z() * f.z());
}

inline vec3 operator*(float const v, vec3 const &e)
{
	return vec3(e.x() * v, e.y() * v, e.z() * v);
}

inline vec3 operator/(vec3 const &e, float const v)
{
	if (abs(v) > 1e-6) return (vec3(e[0] / v, e[1] / v, e[2] / v));
}

inline float dot(vec3 const &e1, vec3 const &e2)
{
	return (e1.e[0] * e2.e[0] + e1.e[1] * e2.e[1] + e1.e[2] * e2.e[2]);
}

inline vec3 cross(vec3 const &e1, vec3 const &e2)
{
	//right-hand
	return vec3(
		(e1.e[1] * e2.e[2] - e1.e[2] * e2.e[1]),
		-(e1.e[0] * e2.e[2] - e1.e[2] * e2.e[0]),
		(e1.e[0] * e2.e[1] - e1.e[1] * e2.e[0])
	);
}

#endif // !_VEC2_


#ifndef _VEC4_
#define _VEC4_

class vec4 {
public:
	float e[4];
	vec4() { _mm_storeu_ps(e, _mm_setzero_ps()); };
#ifdef SIMD
	vec4(float a) {
		_mm_storeu_ps(e, _mm_load_ps1(&a));
	}
	vec4(float a, float b, float c, float d = 0.0f) {
		_mm_storeu_ps(e, { a,b,c,d });
	}
#else
	vec4(float a) { e[0] = a; e[1] = a; e[2] = a; e[3] = a; }
	vec4(float a, float b, float c, float d = 0.0f) { e[0] = a; e[1] = b; e[2] = c; e[3] = d; }

#endif // SIMD


	vec4(const vec3 &a, float b = 0.0f) { e[0] = a[0]; e[1] = a[1]; e[2] = a[2]; e[3] = b; }

	void set(float a, float b, float c, float d) { e[0] = a; e[1] = b; e[2] = c; e[3] = d; }
	float x() const { return e[0]; }
	float y() const { return e[1]; }
	float z() const { return e[2]; }
	float w() const { return e[3]; }

	float r() const { return e[0]; }
	float g() const { return e[1]; }
	float b() const { return e[2]; }
	float a() const { return e[3]; }
	int size() { return 4; }
	const vec4& operator+() const { return *this; }
	vec4 operator-() const { return vec4(-e[0], -e[1], -e[2], -e[3]); }

	float operator[](int i) const { return e[i]; }
	float& operator[](int i) { return e[i]; }

#ifdef SIMD
	vec4& operator+=(const vec4 &v2) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_loadu_ps(v2.e);
		_mm_storeu_ps(this->e, _mm_add_ps(v0, v1));
		return *this;
	}
	vec4& operator+=(const float &v) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_load_ps1(&v);
		_mm_storeu_ps(this->e, _mm_add_ps(v0, v1));
		return *this;
	}
	vec4& operator-=(const vec4 &v2) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_loadu_ps(v2.e);
		_mm_storeu_ps(this->e, _mm_sub_ps(v0, v1));
		return *this;
	}
	vec4& operator-=(const float &v) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_load_ps1(&v);
		_mm_storeu_ps(this->e, _mm_sub_ps(v0, v1));
		return *this;
	}
	vec4& operator*=(const vec4 &v2) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_loadu_ps(v2.e);
		_mm_storeu_ps(this->e, _mm_mul_ps(v0, v1));
		return *this;
	}
	vec4& operator*=(const float &v) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_load_ps1(&v);
		_mm_storeu_ps(this->e, _mm_mul_ps(v0, v1));
		return *this;
	}
	vec4& operator/=(const vec4 &v2) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_loadu_ps(v2.e);
		_mm_storeu_ps(this->e, _mm_div_ps(v0, v1));
		return *this;
	}
	vec4& operator/=(const float &v) {
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_load_ps1(&v);
		_mm_storeu_ps(this->e, _mm_div_ps(v0, v1));
		return *this;
	}

	float square_length() const
	{
		__m128 sums = _mm_dp_ps(_mm_loadu_ps(e), _mm_loadu_ps(e), 0xff);
		return _mm_cvtss_f32(sums);
	}

	vec4 unit() const
	{	
		vec4 result;
		float l = length();	
		__m128 v0 = _mm_loadu_ps(e);
		__m128 v1 = _mm_load_ps1(&l);
		_mm_storeu_ps(result.e, _mm_div_ps(v0, v1));
		return std::move(result);
		/*float l = length();
		if (abs(l) > 1e-6) {
			vec4 result;
			__m128 v0 = _mm_loadu_ps(e);
			__m128 v1 = _mm_load_ps1(&l);			
			_mm_storeu_ps(result.e, _mm_div_ps(v0, v1));
			return std::move(result);
		}
		else {
			return vec4(0.0f);
		}*/
	}
#else
	vec4& operator+=(const vec4 &v2) { e[0] += v2.x(); e[1] += v2.y();  e[2] += v2.z(); e[3] += v2.w(); return *this; }
	vec4& operator+=(const float &v) { for (int i = 0; i < 4; ++i) { e[i] += v; } return *this; }
	vec4& operator-=(const vec4 &v2) { e[0] -= v2.x(); e[1] -= v2.y();  e[2] -= v2.z(); e[3] -= v2.w(); return *this; }
	vec4& operator-=(const float &v) { for (int i = 0; i < 4; ++i) { e[i] -= v; } return *this; }
	vec4& operator*=(const vec4 &v2) { e[0] *= v2.x(); e[1] *= v2.y();  e[2] *= v2.z(); e[3] *= v2.w(); return *this; }
	vec4& operator*=(const float &v) { for (int i = 0; i < 4; ++i) { e[i] *= v; } return *this; }
	vec4& operator/=(const vec4 &v2) { e[0] /= v2.x(); e[1] /= v2.y();  e[2] /= v2.z(); e[3] /= v2.w(); return *this; }
	vec4& operator/=(const float &v) { for (int i = 0; i < 4; ++i) { e[i] /= v; } return *this; }

	float square_length() const
	{
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2] + e[3] * e[3];
	}

	vec4 unit() const
	{
		float l = length();
		if (abs(l) > 1e-6) {
			return vec4(this->e[0] / l, this->e[1] / l, this->e[2] / l, this->e[3] / l);
		}
		else {
			return vec4(0.0f);
		}
	}

#endif // SIMD


	vec3 toVec3()const {
		vec3 temp;
		for (int i = 0; i < 3; ++i) {
			temp[i] = e[i];
		}
		return std::move(temp);
	}

	float length() const
	{
		return sqrt(square_length());
	}

};

inline std::istream& operator>>(std::istream &is, vec4 &t)
{
	is >> t.e[0] >> t.e[1] >> t.e[2] >> t.e[3];
	return is;
}

inline std::ostream& operator<<(std::ostream &os, vec4 &t)
{
	os << t.e[0] << t.e[1] << t.e[2] << t.e[3];
	return os;
}

#ifdef SIMD
inline vec4 operator+(vec4 const &e, float const v)
{
	vec4 result;
	__m128 v0 = _mm_loadu_ps(e.e);
	__m128 v1 = _mm_load_ps1(&v);
	_mm_storeu_ps(result.e, _mm_add_ps(v0, v1));
	return std::move(result);
}

inline vec4 operator+(vec4 const &e, vec4 const &f)
{
	vec4 result;
	__m128 v0 = _mm_loadu_ps(e.e);
	__m128 v1 = _mm_loadu_ps(f.e);
	_mm_storeu_ps(result.e, _mm_add_ps(v0, v1));
	return std::move(result);
}

inline vec4 operator-(vec4 const &e, float const v)
{
	vec4 result;
	__m128 v0 = _mm_loadu_ps(e.e);
	__m128 v1 = _mm_load_ps1(&v);
	_mm_storeu_ps(result.e, _mm_sub_ps(v0, v1));
	return std::move(result);
}

inline vec4 operator-(vec4 const &e, vec4 const &f)
{
	vec4 result;
	__m128 v0 = _mm_loadu_ps(e.e);
	__m128 v1 = _mm_loadu_ps(f.e);
	_mm_storeu_ps(result.e, _mm_sub_ps(v0, v1));
	return std::move(result);
}

inline vec4 operator*(vec4 const &e, float const v)
{
	vec4 result;
	__m128 v0 = _mm_loadu_ps(e.e);
	__m128 v1 = _mm_load_ps1(&v);
	_mm_storeu_ps(result.e, _mm_mul_ps(v0, v1));
	return std::move(result);
}

inline vec4 operator*(vec4 const &e, vec4 const &f)
{
	vec4 result;
	__m128 v0 = _mm_loadu_ps(e.e);
	__m128 v1 = _mm_loadu_ps(f.e);
	_mm_storeu_ps(result.e, _mm_mul_ps(v0, v1));
	return std::move(result);
}

inline vec4 operator*(float const v, vec4 const &e)
{
	vec4 result;
	__m128 v0 = _mm_loadu_ps(e.e);
	__m128 v1 = _mm_load_ps1(&v);
	_mm_storeu_ps(result.e, _mm_mul_ps(v0, v1));
	return std::move(result);
}

inline vec4 operator/(vec4 const &e, float const v)
{
	if (abs(v) > 1e-6) {
		vec4 result;
		__m128 v0 = _mm_loadu_ps(e.e);
		__m128 v1 = _mm_load_ps1(&v);
		_mm_storeu_ps(result.e, _mm_div_ps(v0, v1));
		return std::move(result);
	}
	else {
		return vec4(0.0f);
	}
}

inline float dot(vec4 const &e1, vec4 const &e2)
{
	__m128 sums = _mm_dp_ps(_mm_loadu_ps(e1.e), _mm_loadu_ps(e2.e), 0xff);
	return _mm_cvtss_f32(sums);
}

#else
inline vec4 operator+(vec4 const &e, float const v)
{
	return vec4(e.x() + v, e.y() + v, e.z() + v, e.w() + v);
}

inline vec4 operator+(vec4 const &e, vec4 const &f)
{
	return vec4(e.x() + f.x(), e.y() + f.y(), e.z() + f.z(), e.w() + f.w());
}

inline vec4 operator-(vec4 const &e, float const v)
{
	return vec4(e.x() - v, e.y() - v, e.z() - v, e.w() - v);
}

inline vec4 operator-(vec4 const &e, vec4 const &f)
{
	return vec4(e.x() - f.x(), e.y() - f.y(), e.z() - f.z(), e.w() - f.w());
}

inline vec4 operator*(vec4 const &e, float const v)
{
	return vec4(e.x() * v, e.y() * v, e.z() * v, e.w() * v);
}

inline vec4 operator*(vec4 const &e, vec4 const &f)
{
	return vec4(e.x() * f.x(), e.y() * f.y(), e.z() * f.z(), e.w() * f.w());
}

inline vec4 operator*(float const v, vec4 const &e)
{
	return vec4(e.x() * v, e.y() * v, e.z() * v, e.w()*v);
}

inline vec4 operator/(vec4 const &e, float const v)
{
	if (abs(v) > 1e-6) return (vec4(e[0] / v, e[1] / v, e[2] / v, e[3] / v));
}

inline float dot(vec4 const &e1, vec4 const &e2)
{
	return (e1.e[0] * e2.e[0] + e1.e[1] * e2.e[1] + e1.e[2] * e2.e[2] + e1.e[3] * e2.e[3]);
}
#endif // SIMD


#endif // !_VEC2_


#ifndef _MATRIX_
#define _MATRIX_

bool gauss_jordan(float* a, int n); 


class Matrix3x3 {
public:
	float e[3][3];
	int width = 3;
	Matrix3x3() {}
	Matrix3x3(float **in) {
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < width; ++j) {
				*(*(e + j) + i) = *(*(in + j) + i);
			}
		}
	}
	Matrix3x3(const Matrix3x3 &in) {
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < width; ++j) {
				*(*(e + j) + i) = *(*(in.e + j) + i);
			}
		}
	}
	Matrix3x3 mat_mult(const Matrix3x3 &b) const {
		Matrix3x3 temp;
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < width; ++j) {
				temp.e[i][j] = 0;
				for (int k = 0; k < width; ++k) {
					temp.e[i][j] += e[i][k] * b.e[k][j];
				}
			}
		}
		return std::move(temp);
	}
	void LUdecomp(Matrix3x3 &L, Matrix3x3 &U) {
		//Matrix3x3 L, U;
		for (int r = 0; r < 3; ++r) {
			for (int i = r; i < 3; ++i) {
				U.e[r][i] = e[r][i];
				L.e[i][r] = e[i][r];
				L.e[r][i] = 0;
				U.e[i][r] = 0;
				for (int k = 0; k < r; ++k) {
					U.e[r][i] -= L.e[r][k] * U.e[k][i];
					L.e[i][r] -= L.e[i][k] * U.e[k][r];
				}
				L.e[i][r] /= U.e[r][r];
			}
		}
	}

	Matrix3x3 inverse() {
		Matrix3x3 L, U;
		LUdecomp(L, U);

		Matrix3x3 L_inv(L);
		L_inv.e[1][0] *= -1;
		L_inv.e[2][1] *= -1;
		L_inv.e[2][1] = L.e[1][0] * L.e[2][1] - L.e[2][1];

		Matrix3x3 U_inv(U);
		for (int i = 0; i < width; ++i) {
			for (int j = i; j < width; ++j) {
				U_inv.e[i][j] /= U.e[i][i];
			}
		}
		U_inv.e[0][1] *= -1;
		U_inv.e[1][2] *= -1;
		U_inv.e[0][2] = U_inv.e[0][1] * U_inv.e[1][2] - U_inv.e[0][2];

		return std::move(U_inv.mat_mult(L_inv));
	}
};

Matrix3x3 mat_mult(const Matrix3x3& a, const Matrix3x3& b);


class Matrix4x4 {
public:
	float e[4][4];
	static const int width = 4;


	Matrix4x4();
	Matrix4x4(float* in);
	Matrix4x4(const Matrix4x4& in);
	Matrix4x4(float val);
	Matrix4x4(const vec4& l0, const vec4& l1, const vec4& l2, const vec4& l3);

	vec4 mul_vec(const vec4& in);

	Matrix4x4 tranpose();
	Matrix4x4 mat_mult(const Matrix4x4& b);
	void LUdecomp(Matrix4x4& L, Matrix4x4& U);
	Matrix4x4 inverse();
};

Matrix4x4 operator+(const Matrix4x4& a, const Matrix4x4& b);

Matrix4x4 identityMatrix4x4(); 

void matTranslate(Matrix4x4& mat, const vec3& offset); 
void matRotate(Matrix4x4& mat, float theta, const vec3& axis); 
void matScale(Matrix4x4& mat, const vec3& scale); 

void applyTrans(Matrix4x4& model_matrix, std::vector<vec4>& pointArray); 

#endif // !_MATRIX_
