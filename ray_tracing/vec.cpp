#pragma once

#include "vec.h"


bool gauss_jordan(float* a, int n) {
	//d用于全选主元过程
	//t用于存储各种临时数据
	float d, t;
	//i用于遍历时的行号
	//j用于遍历时的列号
	//p，q用于存储标号
	//js用于存储全选主元元素的行号
	//k用于表示第几次操作
	int i, j, p, q, k;
	int* js = new int[n];
	int* is = new int[n];
	int two_width = n * 2;
	float* c = new float[two_width * n];
	//构造增广矩阵

	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			c[i * two_width + j] = a[i * n + j];
		}

		for (j = n; j < two_width; ++j) {
			if (j - n == i) { c[i * two_width + j] = 1; }
			else { c[i * two_width + j] = 0; }
		}
	}

	//消元过程
	for (k = 0; k < n; ++k) {
		//全选主元
		d = 0.0;
		for (i = k; i < n; ++i) {
			for (j = k; j < n; ++j) {
				t = c[i * two_width + j];
				if (fabs(t) > d) {
					d = fabs(t);
					js[k] = j;//记录列变换操作
					is[k] = i;//记录行号
				}
			}
		}

		//如果全0则奇异矩阵，无法计算		
		if (fabs(d) < 1e-6) {
			delete[] is;
			delete[] js;
			delete[] c;
			return false;
		}

		//行交换
		if (is[k] != k) {
			for (j = 0; j < two_width; ++j) {
				p = k * two_width + j;
				q = is[k] * two_width + j;
				t = c[p];
				c[p] = c[q]; c[q] = t;
			}
		}
		//列交换
		if (js[k] != k) {
			for (i = 0; i < n; ++i) {
				p = i * two_width + k;
				q = i * two_width + js[k];
				t = c[p]; c[p] = c[q]; c[q] = t;
			}
		}
		//归一化
		{
			d = c[k * two_width + k];
			for (j = k + 1; j < two_width; ++j) {
				p = k * two_width + j;
				c[p] = c[p] / d;
			}
		}
		//消元
		for (i = 0; i < n; ++i) {
			if (i != k) {
				for (j = k + 1; j < two_width; ++j) {
					p = i * two_width + j;
					q = i * two_width + k;
					c[p] = c[p] - c[q] * c[k * two_width + j];
				}
			}
		}
	}
	//恢复
	for (k = n - 1; k >= 0; --k) {
		for (j = 0; j < two_width; ++j) {//换列 : 其实是在换行
			t = c[k * two_width + j];
			c[k * two_width + j] = c[js[k] * two_width + j];
			c[js[k] * two_width + j] = t;
		}
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i * n + j] = c[i * two_width + j + n];
		}
	}
	//结束收工
	delete[] is;
	delete[] js;
	delete[] c;
	return true;

}


Matrix3x3 mat_mult(const Matrix3x3& a, const Matrix3x3& b) {
	Matrix3x3 temp;
	temp = a.mat_mult(b);
	return std::move(temp);
}



#ifdef SIMD
Matrix4x4::Matrix4x4() {
	for (int i = 0; i < 4; ++i) {
		_mm_storeu_ps(e[i], _mm_setzero_ps());
	}
}
Matrix4x4::Matrix4x4(float* in) {
	for (int i = 0; i < 4; ++i) {
		_mm_storeu_ps(e[i], _mm_loadu_ps(in + 4 * i));
	}
}
Matrix4x4::Matrix4x4(const Matrix4x4& in) {
	for (int i = 0; i < 4; ++i) {
		_mm_storeu_ps(e[i], _mm_loadu_ps(in.e[i]));
	}
}
Matrix4x4::Matrix4x4(float val) {
	__m128 data = _mm_load_ps1(&val);
	for (int i = 0; i < 4; ++i) {
		_mm_storeu_ps(e[i], data);
	}
}
Matrix4x4::Matrix4x4(const vec4& l0, const vec4& l1, const vec4& l2, const vec4& l3) {
	_mm_storeu_ps(e[0], _mm_loadu_ps(l0.e));
	_mm_storeu_ps(e[1], _mm_loadu_ps(l1.e));
	_mm_storeu_ps(e[2], _mm_loadu_ps(l2.e));
	_mm_storeu_ps(e[3], _mm_loadu_ps(l3.e));
}

vec4 Matrix4x4::mul_vec(const vec4& in) {
	vec4 temp;
	__m128 vec = _mm_loadu_ps(in.e);

	for (int i = 0; i < width; ++i) {
		__m128 sums = _mm_dp_ps(vec, _mm_loadu_ps(e[i]), 0xff);
		temp[i] = _mm_cvtss_f32(sums);
	}
	return std::move(temp);
}

#else
Matrix4x4::Matrix4x4() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			e[i][j] = 0.f;
		}
	}
}
Matrix4x4::Matrix4x4(float* in) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			*(*(e + j) + i) = *(in + 4 * j + i);
		}
	}
}
Matrix4x4::Matrix4x4(const Matrix4x4& in) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			*(*(e + j) + i) = *(*(in.e + j) + i);
		}
	}
}
Matrix4x4::Matrix4x4(float val) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			*(*(e + j) + i) = val;
		}
	}
}
Matrix4x4::Matrix4x4(const vec4& l0, const vec4& l1, const vec4& l2, const vec4& l3) {
	for (int i = 0; i < width; ++i) {
		e[0][i] = l0[i];
		e[1][i] = l1[i];
		e[2][i] = l2[i];
		e[3][i] = l3[i];
	}
}

vec4 Matrix4x4::mul_vec(const vec4& in) {
	vec4 temp;
	for (int i = 0; i < width; ++i) {
		temp.e[i] = 0;
		for (int m = 0; m < width; ++m) {
			temp.e[i] += e[i][m] * in.e[m];
		}
	}
	return std::move(temp);
}
#endif // SIMD

Matrix4x4 Matrix4x4::tranpose() {
	Matrix4x4 temp;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			*(*(temp.e + j) + i) = *(*(this->e + i) + j);
		}
	}
	return std::move(temp);
}

Matrix4x4 Matrix4x4::mat_mult(const Matrix4x4& b) {
	Matrix4x4 temp(0.0f);
	for (int i = 0; i < width; ++i) {
		for (int k = 0; k < width; ++k) {
			for (int j = 0; j < width; ++j) {
				temp.e[i][j] += e[i][k] * b.e[k][j];
			}
		}
	}
	return std::move(temp);
}

void Matrix4x4::LUdecomp(Matrix4x4& L, Matrix4x4& U) {
	//Matrix4x4 L, U;
	for (int r = 0; r < width; ++r) {
		for (int i = r; i < width; ++i) {
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

Matrix4x4 Matrix4x4::inverse() {//待改
					 //Matrix4x4 t(*this);
	float temp[16];
	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < width; ++j) {
			temp[i * 4 + j] = this->e[i][j];
		}
	}
	Matrix4x4 t;
	if (gauss_jordan(temp, width)) {
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < width; ++j) {
				t.e[i][j] = temp[i * 4 + j];
			}
		}
	}
	else {
		t = *this;
	}


	return t;
}


#ifdef SIMD
Matrix4x4 operator+(const Matrix4x4& a, const Matrix4x4& b) {
	Matrix4x4 t(a);
	for (int i = 0; i < 4; ++i) {
		__m128 m0 = _mm_loadu_ps(t.e[i]);
		__m128 m1 = _mm_loadu_ps(b.e[i]);
		_mm_storeu_ps(t.e[i], _mm_add_ps(m0, m1));
	}
	return std::move(t);
}

#else
Matrix4x4 operator+(const Matrix4x4& a, const Matrix4x4& b) {
	Matrix4x4 t(a);
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			*(*(t.e + i) + j) += *(*(b.e + i) + j);
		}
	}
	return std::move(t);
}
#endif // SIMD



Matrix4x4 identityMatrix4x4() {
	Matrix4x4 temp(0.f);
	for (int i = 0; i < 4; ++i) {
		temp.e[i][i] = 1.f;
	}
	return temp;
}

void matTranslate(Matrix4x4& mat, const vec3& offset) {
	Matrix4x4 t = identityMatrix4x4();
	for (int i = 0; i < 3; ++i) {
		t.e[i][3] += offset[i];
	}
	mat = t.mat_mult(mat);
}

void matRotate(Matrix4x4& mat, float theta, const vec3& axis) {
	theta = theta * M_PI / 180.f;
	float cos_theta = cos(theta), sin_theta = sin(theta);
	float one_minus_cos = 1.f - cos_theta;
	vec3 u = axis.unit();
	float u_pow[3] = { u.x() * u.x(), u.y() * u.y(), u.z() * u.z() };
	float u_xy = u.x() * u.y(), u_xz = u.x() * u.z(), u_yz = u.y() * u.z();

	Matrix4x4 t = identityMatrix4x4();
	for (int i = 0; i < 3; ++i) {
		t.e[i][i] = u_pow[i] * one_minus_cos + cos_theta;
	}
	t.e[0][1] = u_xy * one_minus_cos - u.z() * sin_theta; t.e[0][2] = u_xz * one_minus_cos + u.y() * sin_theta;
	t.e[1][0] = u_xy * one_minus_cos + u.z() * sin_theta; t.e[1][2] = u_yz * one_minus_cos - u.x() * sin_theta;
	t.e[2][0] = u_xz * one_minus_cos - u.y() * sin_theta; t.e[2][1] = u_yz * one_minus_cos + u.x() * sin_theta;

	Matrix4x4 traslate = identityMatrix4x4(), traslate_inv = identityMatrix4x4();
	matTranslate(traslate, axis);
	matTranslate(traslate_inv, -axis);
	mat = traslate.mat_mult(mat);
	mat = t.mat_mult(mat);
	mat = traslate_inv.mat_mult(mat);
}

void matScale(Matrix4x4& mat, const vec3& scale) {
	Matrix4x4 t = identityMatrix4x4();
	for (int i = 0; i < 3; ++i) {
		t.e[i][i] *= scale[i];
	}
	mat = t.mat_mult(mat);
}


#ifdef SIMD
void applyTrans(Matrix4x4& model_matrix, std::vector<vec4>& pointArray) {
	__m128 data[4];
	for (int j = 0; j < 4; ++j) {
		data[j] = _mm_loadu_ps(model_matrix.e[j]);
	}

	for (int i = 0; i < pointArray.size(); ++i) {
		__m128 vec = _mm_loadu_ps(pointArray[i].e);
		for (int j = 0; j < 4; ++j) {
			__m128 sums = _mm_dp_ps(data[j], vec, 0xff);
			pointArray[i][j] = _mm_cvtss_f32(sums);
		}
	}
}
#else
void applyTrans(Matrix4x4& model_matrix, std::vector<vec4>& pointArray) {
	for (int i = 0; i < pointArray.size(); ++i) {
		pointArray[i] = model_matrix.mul_vec(pointArray[i]);
	}
}
#endif // SIMD
