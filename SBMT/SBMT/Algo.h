#ifndef ALGORITHM
#define ALGORITHM
/*!
*      \File Algo.h
*      \brief general algorithms
*	   \author Wei Feng
*      \date 11/21/2024
*/

#include <vector>
#include<cassert>
#include<algorithm>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#define EIGEN_STRONG_INLINE inline
#include <Eigen/Dense>
#include <set>
using namespace Eigen;
using namespace cv;
using namespace std;
#define PI 3.14159265358979323846264338327950288419716939937510
#define MAXV 1048576
#define TOLERANCE 1e-5
#define EPSILON 1e-9
struct Pt {
	double x, y;//the coordinates on the integral image
	uchar val;
	double ori, oci;

	Pt() : x(0), y(0), val(0), ori(0), oci(0) {}
	Pt(float x, float y, uchar val, float ori, float oci)
		: x(x), y(y), val(val), ori(ori), oci(oci) {}

	double norm()
	{
		return sqrt(pow(x, 2) + pow(y, 2));
	}

	// Vector subtraction
	Pt operator-(const Pt& other) const {
		return Pt(x - other.x
			, y - other.y
			, val - other.val
			, ori - other.ori
			, oci - other.oci);
	}

	// Dot product
	float dot(const Pt& other) const {
		return x * other.x + y * other.y;
	}

	float euclidDist(const Pt& other) const {
		return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
	}

	float BMPDist(const Pt& other) const {
		return sqrt(pow(oci - other.oci, 2) + pow(ori - other.ori, 2));
	}
};

// 自定义比较函数，确保 Point2d 在 set 中唯一
struct ComparePoints {
	bool operator()(const Pt& p1, const Pt& p2) const {
		// VS2013 可能不支持直接比较浮点数，所以用 epsilon 进行浮点比较
		const float epsilon = 1e-6f;
		if (fabs(p1.x - p2.x) > epsilon) {
			return p1.x < p2.x;
		}
		return p1.y < p2.y;
	}
};

struct FootResult {
	cv::Point2d adjustedPoint;
	bool adjusted;
};

template<class T>
void MatrixIdentity(T* matrix, int nrows)
{ //  matrix -  nrows*nrows
	memset(matrix, 0, nrows*nrows * sizeof(T));
	for (int i = 0; i<nrows; ++i)
		matrix[i*nrows + i] = 1;
}

template<class TT1, class TT2, class TT3>//>
void Doolittle(TT1* aa, TT2* bb, TT3* xx, int rows)
{// aa * xx = bb        root - xx[rows][rows]
	int k, i, j, t, ik;
	int* M = new int[rows];
	double  *s, *l, *u, *a, *b;
	double temp, smax = 0, *y, *x;
	s = new double[rows];
	l = new double[rows*rows];
	u = new double[rows*rows];
	a = new double[rows*rows];
	b = new double[rows];
	y = new double[rows];
	x = new double[rows];
	//  QA  =  LU
	for (i = 0; i<rows; ++i)
	{
		M[i] = 0;
		for (j = 0; j<rows; ++j)
		{
			a[i*rows + j] = aa[i*rows + j];
		}
	}
	for (k = 0; k<rows; ++k)
	{
		for (i = k; i<rows; ++i)
		{
			s[i] = a[i*rows + k];
			for (t = 0; t < k; ++t)
			{
				s[i] -= l[i*rows + t] * u[t*rows + k];
			}

			if (i == k)
			{
				smax = s[i];
				ik = i;
			}
			if (fabs(smax)<fabs(s[i]))
			{
				smax = s[i];
				ik = i;
			}
		}
		M[k] = ik;
		if (ik != k)
		{
			for (t = 0; t<k; ++t)
			{
				temp = l[k*rows + t];
				l[k*rows + t] = l[ik*rows + t];
				l[ik*rows + t] = temp;
			}
			for (t = k; t<rows; ++t)
			{
				temp = a[k*rows + t];
				a[k*rows + t] = a[ik*rows + t];
				a[ik*rows + t] = temp;
			}
			temp = s[k];
			s[k] = s[ik];
			s[ik] = temp;
		}
		u[k*rows + k] = s[k];
		if (k<rows - 1)
		{
			for (j = k + 1; j<rows; ++j)
			{
				u[k*rows + j] = a[k*rows + j];
				for (t = 0; t < k; ++t)
				{
					u[k*rows + j] -= l[k*rows + t] * u[t*rows + j];
				}

			}
			for (i = k + 1; i < rows; ++i)
			{
				l[i*rows + k] = s[i] / (u[k*rows + k] + 0.00001);
			}

		}
	}
	//Qb  =  Ly   AND   Ux  =   y
	for (j = 0; j<rows; ++j)
	{
		for (i = 0; i < rows; ++i)
		{
			b[i] = bb[i*rows + j];
		}

		for (k = 0; k<rows - 1; ++k)
		{
			t = M[k];
			temp = b[k];
			b[k] = b[t];
			b[t] = temp;
		}
		y[0] = b[0];
		for (i = 1; i<rows; ++i)
		{
			y[i] = b[i];
			for (t = 0; t < i; ++t)
			{
				y[i] -= l[i*rows + t] * y[t];
			}

		}
		x[rows - 1] = y[rows - 1] / (u[rows*rows - 1] + 0.00001);
		for (i = rows - 2; i>-1; --i)
		{
			x[i] = y[i];
			for (t = i + 1; t < rows; ++t)
			{
				x[i] -= u[i*rows + t] * x[t];
			}

			x[i] /= (u[i*rows + i] + 0.00001);
		}
		for (i = 0; i<rows; ++i)
		{
			xx[i*rows + j] = x[i];
		}
	}
	delete[]M;
	delete[]s;
	delete[]l;
	delete[]u;
	delete[]a;
	delete[]b;
	delete[]y;
	delete[]x;
	M = NULL;
	s = NULL;
	l = NULL;
	u = NULL;
	a = NULL;
	b = NULL;
	y = NULL;
	x = NULL;
}


// invertible matrice
template<class TT1, class TT2>
void MatrixAnti(TT1* Matrix, TT2* MatrixA, int rows)
{//  Matrix * MatrixA = I          I = E
	double* E = new double[rows*rows];
	MatrixIdentity(E, rows);

	//Doolittle solution
	Doolittle(Matrix, E, MatrixA, rows);
	delete[]E;
	E = NULL;
}

template<class TT1, class TT2, class TT3>
void  __declspec(dllexport) MultMatrix(TT1* M, TT2* M1, TT3* M2, int rows, int soms, int cols)
{
	//M - rows*soms    M1 - soms*cols   M2 - rows*cols
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j<cols; ++j)
		{
			M2[i*cols + j] = 0.0;
			for (int k = 0; k < soms; ++k)
			{
				M2[i*cols + j] += M[i*soms + k] * M1[k*cols + j];
			}

		}
	}

}

// Function to check if point P is inside the triangle formed by A, B, C
inline bool isPointInTriangle(const Pt& A, const Pt& B, const Pt& C, const Pt& P) {
	// Calculate vectors
	Pt v0 = C - A;
	Pt v1 = B - A;
	Pt v2 = P - A;

	// Calculate dot products
	float dot00 = v0.dot(v0);
	float dot01 = v0.dot(v1);
	float dot02 = v0.dot(v2);
	float dot11 = v1.dot(v1);
	float dot12 = v1.dot(v2);

	// Calculate the inverse denominator
	float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);

	// Calculate barycentric coordinates u and v
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point P is inside the triangle
	return (u >= 0) && (v >= 0) && (u + v <= 1);
}

// Function to check if point p is inside quadrilateral
inline bool isPointInsideQuadrilateral(Pt p, const vector<Pt>& quad)
{
	if (quad.size() != 4)
	{
		std::cerr << "Error: Quadrilateral must have exactly 4 points.\n";
		return false;
	}

	// Step 1: Check bounding box
	double minX = std::min({ quad[0].x, quad[1].x, quad[2].x, quad[3].x });
	double maxX = std::max({ quad[0].x, quad[1].x, quad[2].x, quad[3].x });
	double minY = std::min({ quad[0].y, quad[1].y, quad[2].y, quad[3].y });
	double maxY = std::max({ quad[0].y, quad[1].y, quad[2].y, quad[3].y });

	if (p.x < minX || p.x > maxX || p.y < minY || p.y > maxY)
	{
		return false;
	}

	if (isPointInTriangle(quad[0], quad[1], quad[2], p))
	{
		return true;
	}
	else if (isPointInTriangle(quad[0], quad[2], quad[3], p))
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline int classify2(vector<int> &t1, double &cluster1, vector<int> &t2
	, double &cluster2, vector<Pt> & ps4)
{
	cluster1 = ps4[0].y;
	t1.push_back(0);
	int t_num = 1;
	for (int i = 1; i != 4; ++i)
	{
		if (abs(cluster1 - ps4[i].y) > 0.5)
		{
			cluster2 = ps4[i].y;
			t2.push_back(i);
			t_num = 2;
		}
		else
		{
			t1.push_back(i);
		}
	}

	return t_num;
}

inline void getRandomNum(Point2d &coors)
{
	coors.x = (std::rand() / (double)RAND_MAX) - 0.5;
	coors.y = (std::rand() / (double)RAND_MAX) - 0.5;

	coors.x /= 2;
	coors.y /= 2;
}

inline uchar BilineInterpola(uchar g00, uchar g01, uchar g10, uchar g11
	, double a, double b)
{
	return b*(a*g11 + (1 - a)*g01) + (1 - b)*(a*g10 + (1 - a)*g00);
}

inline uchar getGray(Mat &mat, Point2d p)
{
	Point2d pT(p.x - 0.5, p.y - 0.5);
	uchar g00 = mat.at<uchar>(floor(pT.y), floor(pT.x));
	uchar g01 = mat.at<uchar>(floor(pT.y), floor(pT.x) + 1);
	uchar g10 = mat.at<uchar>(floor(pT.y) + 1, floor(pT.x));
	uchar g11 = mat.at<uchar>(floor(pT.y) + 1, floor(pT.x) + 1);
	double a = pT.y - floor(pT.y);
	double b = pT.x - floor(pT.x);
	return BilineInterpola(g00, g01, g10, g11
		, a, b);
}

inline uchar getGray(Mat &mat, double x, double y)
{
	Point2d pT(x - 0.5, y - 0.5);
	uchar g00 = mat.at<uchar>(floor(pT.y), floor(pT.x));
	uchar g01 = mat.at<uchar>(floor(pT.y), floor(pT.x) + 1);
	uchar g10 = mat.at<uchar>(floor(pT.y) + 1, floor(pT.x));
	uchar g11 = mat.at<uchar>(floor(pT.y) + 1, floor(pT.x) + 1);
	double a = pT.y - floor(pT.y);
	double b = pT.x - floor(pT.x);
	return BilineInterpola(g00, g01, g10, g11
		, a, b);
}

inline void addTri(vector<vector<Point2d>> &tri_set
	, Point2d p1, Point2d p2, Point2d p3, Point2d p4, Point2d p5, Point2d p6)
{
	vector<Point2d> t1, t2, t3, t4;
	t1.push_back(p1);
	t1.push_back(p3);
	t1.push_back(p2);
	tri_set.push_back(t1);
	t2.push_back(p2);
	t2.push_back(p3);
	t2.push_back(p4);
	tri_set.push_back(t2);
	t3.push_back(p1);
	t3.push_back(p5);
	t3.push_back(p3);
	tri_set.push_back(t3);
	t4.push_back(p3);
	t4.push_back(p6);
	t4.push_back(p4);
	tri_set.push_back(t4);
}


inline void addTri(vector<vector<int>> &tri_set
	, int i1, int i2, int i3, int i4, int i5, int i6)
{
	vector<int> t1, t2, t3, t4;
	t1.push_back(i1);
	t1.push_back(i3);
	t1.push_back(i2);
	tri_set.push_back(t1);
	t2.push_back(i2);
	t2.push_back(i3);
	t2.push_back(i4);
	tri_set.push_back(t2);
	t3.push_back(i1);
	t3.push_back(i5);
	t3.push_back(i3);
	tri_set.push_back(t3);
	t4.push_back(i3);
	t4.push_back(i6);
	t4.push_back(i4);
	tri_set.push_back(t4);
}


inline void getReguTris(vector<Point2d> &ps, vector<uchar> &gs
	, vector<vector<Point2d>> &tri_set, Mat &mat)
{
	double mul = sqrt(0.9);
	double offset = 5.0;
	uchar g00, g01, g10, g11, gc;
	Point2d p1, p2, p3, p4, p5, p6;
	int i1, i2, i3, i4, i5, i6;

	Point2d p;
	int r = mat.rows;
	int c = mat.cols;

	double xStep = 1.0 * mul;
	double yStep = sqrt(3)*mul;

	int xN = floor(c / xStep) - 10.0;
	int yN = floor(r / yStep) - 10.0;

	for (int i = 0; i != yN; ++i)
	{
		for (int j = 0; j != xN; ++j)
		{
			if (0 == i && 0 == j)
			{
				p1.x = 0.5*mul + j*xStep + offset;
				p1.y = 0 * mul + i*yStep + offset;
				p2.x = 0 * mul + j*xStep + offset;
				p2.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p3.x = 0.5 * mul + j*xStep + offset;
				p3.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p1);
				gc = getGray(mat, p1);
				gs.push_back(gc);

				ps.push_back(p2);
				gc = getGray(mat, p2);
				gs.push_back(gc);

				ps.push_back(p3);
				gc = getGray(mat, p3);
				gs.push_back(gc);
			}
			else if (0 == i)
			{
				p4.x = 0.5*mul + j*xStep + offset;
				p4.y = 0 * mul + i*yStep + offset;
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p4);
				gc = getGray(mat, p4);
				gs.push_back(gc);
				i5 = ps.size() - 1;

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);
				i3 = ps.size() - 1;

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
				i6 = ps.size() - 1;

				i1 = i5 - 3;
				i2 = i3 - 3;
				i4 = i6 - 3;
				addTri(tri_set, ps[i1], ps[i2], ps[i3], ps[i4], ps[i5], ps[i6]);
			}
			else if (0 == j)
			{
				p2.x = 0 * mul + j*xStep + offset;
				p2.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p3.x = 0.5 * mul + j*xStep + offset;
				p3.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p2);
				gc = getGray(mat, p2);
				gs.push_back(gc);

				ps.push_back(p3);
				gc = getGray(mat, p3);
				gs.push_back(gc);
			}
			else if (1 == i)
			{
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);
				i3 = ps.size() - 1;

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
				i6 = ps.size() - 1;

				i1 = i3 - j * 2 - (xN - j) * 3 - 1;
				i2 = i3 - 2;
				i4 = i6 - 2;
				i5 = i3 - j * 2 - (xN - j - 1) * 3 - 1;
				addTri(tri_set, ps[i1], ps[i2], ps[i3], ps[i4], ps[i5], ps[i6]);
			}
			else
			{
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);
				i3 = ps.size() - 1;

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
				i6 = ps.size() - 1;

				i1 = i3 - j * 2 - (xN - j) * 2 - 1;
				i2 = i3 - 2;
				i4 = i6 - 2;
				i5 = i3 - j * 2 - (xN - j - 1) * 2 - 1;
				addTri(tri_set, ps[i1], ps[i2], ps[i3], ps[i4], ps[i5], ps[i6]);
			}

		}
	}
}


inline void getReguTris(vector<Point2d> &ps, vector<uchar> &gs
	, vector<vector<int>> &tri_set, Mat &mat)
{
	double mul = sqrt(0.45);
	double offset = 5.0;
	uchar g00, g01, g10, g11, gc;
	Point2d p1, p2, p3, p4, p5, p6;
	int i1, i2, i3, i4, i5, i6;

	Point2d p;
	int r = mat.rows;
	int c = mat.cols;

	double xStep = 1.0 * mul;
	double yStep = sqrt(3)*mul;

	int xN = floor(c / xStep) - 10.0;
	int yN = floor(r / yStep) - 10.0;

	for (int i = 0; i != yN; ++i)
	{
		for (int j = 0; j != xN; ++j)
		{
			if (0 == i && 0 == j)
			{
				p1.x = 0.5*mul + j*xStep + offset;
				p1.y = 0 * mul + i*yStep + offset;
				p2.x = 0 * mul + j*xStep + offset;
				p2.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p3.x = 0.5 * mul + j*xStep + offset;
				p3.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p1);
				gc = getGray(mat, p1);
				gs.push_back(gc);

				ps.push_back(p2);
				gc = getGray(mat, p2);
				gs.push_back(gc);

				ps.push_back(p3);
				gc = getGray(mat, p3);
				gs.push_back(gc);
			}
			else if (0 == i)
			{
				p4.x = 0.5*mul + j*xStep + offset;
				p4.y = 0 * mul + i*yStep + offset;
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p4);
				gc = getGray(mat, p4);
				gs.push_back(gc);
				i5 = ps.size() - 1;

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);
				i3 = ps.size() - 1;

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
				i6 = ps.size() - 1;

				i1 = i5 - 3;
				i2 = i3 - 3;
				i4 = i6 - 3;
				addTri(tri_set, i1, i2, i3, i4, i5, i6);
			}
			else if (0 == j)
			{
				p2.x = 0 * mul + j*xStep + offset;
				p2.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p3.x = 0.5 * mul + j*xStep + offset;
				p3.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p2);
				gc = getGray(mat, p2);
				gs.push_back(gc);

				ps.push_back(p3);
				gc = getGray(mat, p3);
				gs.push_back(gc);
			}
			else if (1 == i)
			{
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);
				i3 = ps.size() - 1;

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
				i6 = ps.size() - 1;

				i1 = i3 - j * 2 - (xN - j) * 3 - 1;
				i2 = i3 - 2;
				i4 = i6 - 2;
				i5 = i3 - j * 2 - (xN - j - 1) * 3 - 1;
				addTri(tri_set, i1, i2, i3, i4, i5, i6);
			}
			else
			{
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);
				i3 = ps.size() - 1;

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
				i6 = ps.size() - 1;

				i1 = i3 - j * 2 - (xN - j) * 2 - 1;
				i2 = i3 - 2;
				i4 = i6 - 2;
				i5 = i3 - j * 2 - (xN - j - 1) * 2 - 1;
				addTri(tri_set, i1, i2, i3, i4, i5, i6);
			}

		}
	}
}


inline void getReguPs(vector<Point2d> &ps, vector<uchar> &gs, Mat &mat)
{
	double mul = 1.0;
	double offset = 5.0;
	uchar g00, g01, g10, g11, gc;
	Point2d p1, p2, p3, p4, p5, p6;

	Point2d p;
	int r = mat.rows;
	int c = mat.cols;

	double xStep = 1.0 * mul;
	double yStep = sqrt(3)*mul;

	int xN = floor(c / xStep) - 10;
	int yN = floor(r / yStep) - 10;

	for (int i = 0; i != yN; ++i)
	{
		for (int j = 0; j != xN; ++j)
		{
			if (0 == i)
			{
				p4.x = 0.5*mul + j*xStep + offset;
				p4.y = 0 * mul + i*yStep + offset;
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p4);
				gc = getGray(mat, p4);
				gs.push_back(gc);

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
			}
			else
			{
				p5.x = 0 * mul + j*xStep + offset;
				p5.y = sqrt(3)*mul / 2 + i*yStep + offset;
				p6.x = 0.5 * mul + j*xStep + offset;
				p6.y = sqrt(3)*mul + i*yStep + offset;

				ps.push_back(p5);
				gc = getGray(mat, p5);
				gs.push_back(gc);

				ps.push_back(p6);
				gc = getGray(mat, p6);
				gs.push_back(gc);
			}

		}
	}
}

template<class T>
double crossProduct(const T& p1, const T& p2, const T& p3)
{
	return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

// 判断两个点是否相等（考虑浮点误差）
template<class T>
bool __declspec(dllexport) isSamePoint(const T &P1, const T &P2)
{
	return (abs(P1.x - P2.x) < EPSILON) && (abs(P1.y - P2.y) < EPSILON);
}

// 判断点是否在线段上
template<class T>
bool __declspec(dllexport) containsPoint(const T &lineP, const T &lineQ, const T &point)
{
	// 直接检查 point 是否是端点之一
	if (isSamePoint(point, lineP) || isSamePoint(point, lineQ))
	{
		return true; // 点是端点
	}

	// 计算叉积，判断点是否在线段所在直线上
	double crossPro = crossProduct(lineP, lineQ, point);
	if (std::fabs(crossPro) > EPSILON) {
		return false; // 叉积不为0，点不在线段所在直线上
	}

	// 计算点积，判断点是否在线段的范围内
	double dotProduct = (point.x - lineP.x) * (lineQ.x - lineP.x)
		+ (point.y - lineP.y) * (lineQ.y - lineP.y);
	if (dotProduct < 0.0f) {
		return false; // 点在延长线的反方向
	}

	// 计算线段的平方长度
	double squaredLengthBA = (lineQ.x - lineP.x) * (lineQ.x - lineP.x)
		+ (lineQ.y - lineP.y) * (lineQ.y - lineP.y);
	if (dotProduct > squaredLengthBA) {
		return false; // 点超出了线段的范围
	}

	return true; // 点在线段上
}


template<class T>
int __declspec(dllexport) pointInTriangle(const T& A, const T& B, const T& C, const T& p)
{
	double c1 = crossProduct(A, B, p);
	double c2 = crossProduct(B, C, p);
	double c3 = crossProduct(C, A, p);

	// 判断是否完全在内部（所有叉积同号）
	if ((c1 > 0 && c2 > 0 && c3 > 0) || (c1 < 0 && c2 < 0 && c3 < 0))
		return 1; // Inside triangle

	// 判断是否在线段上
	if (std::abs(c1) < EPSILON && containsPoint(A, B, p)) return 2;  // On edge AB
	if (std::abs(c2) < EPSILON && containsPoint(B, C, p)) return 3;  // On edge BC
	if (std::abs(c3) < EPSILON && containsPoint(C, A, p)) return 4;  // On edge CA

	return 0; // Outside triangle
}

// Check if the line segment p1p2 lies within the bounding box of the triangle ABC
inline bool isInsideBoundingBox(const Point2d& A, const Point2d& B, const Point2d& C
	, const Point2d& p1, const Point2d& p2)
{
	double minX = std::min({ A.x, B.x, C.x });
	double maxX = std::max({ A.x, B.x, C.x });
	double minY = std::min({ A.y, B.y, C.y });
	double maxY = std::max({ A.y, B.y, C.y });
	return !(p1.x < minX || p1.x > maxX || p1.y < minY || p1.y > maxY ||
		p2.x < minX || p2.x > maxX || p2.y < minY || p2.y > maxY);
}

// 检查两条线段是否有顶点重合
inline bool isVertexOverlap(const Point2d &A, const Point2d &B
	, const Point2d &C, const Point2d &D, Point2d &intersect)
{
	if (isSamePoint(A, C))
	{
		intersect = C;
		return true;
	}
	else if (isSamePoint(A, D))
	{
		intersect = D;
		return true;
	}
	else if (isSamePoint(B, C))
	{
		intersect = C;
		return true;
	}
	else if (isSamePoint(B, D))
	{
		intersect = D;
		return true;
	}
	else
	{
		return false;
	}
}

// 计算两条线段是否相交，并返回交点
inline int getIntersection(Point2d A, Point2d B, Point2d C, Point2d D, Point2d &intersec)
{
	if (isVertexOverlap(A, B, C, D, intersec))
	{
		return 2;
	}
	double a = (D.x - C.x) * (C.y - A.y) - (D.y - C.y) * (C.x - A.x);
	double b = (D.x - C.x) * (B.y - A.y) - (D.y - C.y) * (B.x - A.x);
	double c = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);

	// 处理共线情况
	if (fabs(a) < EPSILON && fabs(b) < EPSILON) {
		// 计算线段投影是否重叠
		if (max(A.x, B.x) < min(C.x, D.x) || max(C.x, D.x) < min(A.x, B.x) ||
			max(A.y, B.y) < min(C.y, D.y) || max(C.y, D.y) < min(A.y, B.y))
		{
			return 0; // 共线但不重叠
		}
		return 1; // 共线且重叠
	}

	if (fabs(b) < EPSILON) {
		return 0; // 平行不相交
	}

	double alpha = a / b;
	double beta = c / b;

	// 检查是否在线段范围内
	if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0) {
		return 0; // 交点在延长线上，不是线段相交
	}

	// 计算交点
	intersec.x = A.x + alpha * (B.x - A.x);
	intersec.y = A.y + alpha * (B.y - A.y);

	return 2; // 线段相交
}


inline bool doIntersect(const Point2d& p1, const Point2d& p2,
	const Point2d& p3, const Point2d& p4, Point2d& intersection)
{
	int res = getIntersection(p1, p2, p3, p4, intersection);

	if (0 == res)
	{
		return false;
	}
	else
	{
		if (1 == res)
		{
			// 处理共线情况：检查端点是否在线段上
			if (containsPoint(p1, p2, p3))
			{
				intersection = p3;
				return true;
			}
			if (containsPoint(p1, p2, p4))
			{
				intersection = p4;
				return true;
			}
			if (containsPoint(p3, p4, p1))
			{
				intersection = p1;
				return true;
			}
			if (containsPoint(p3, p4, p2))
			{
				intersection = p2;
				return true;
			}

			return false;
		}
		else
		{
			return true;
		}
	}
}

// 计算线段与三角形的交点，返回交点列表和 mask
inline vector<Point2d> triaLineSegIntersec(
	const Point2d& A, const Point2d& B, const Point2d& C,
	const Point2d& p1, const Point2d& p2, int &mask)
{
	std::vector<Point2d> intersections;
	mask = 0;

	// 分别计算 p1p2 与三角形三条边的交点
	Point2d intersection1, intersection2, intersection3;
	bool found1 = doIntersect(A, B, p1, p2, intersection1);
	bool found2 = doIntersect(B, C, p1, p2, intersection2);
	bool found3 = doIntersect(C, A, p1, p2, intersection3);

	if (found1) {
		intersections.push_back(intersection1);
		mask |= 2;  // AB
	}
	if (found2) {
		intersections.push_back(intersection2);
		mask |= 4;  // BC
	}
	if (found3) {
		intersections.push_back(intersection3);
		mask |= 8;  // CA
	}

	return intersections;
}

template<class T>
int __declspec(dllexport) isExi(std::vector<T>& es, const T& e) {
	auto it = std::find(es.begin(), es.end(), e);

	if (it != es.end()) {
		return static_cast<int>(std::distance(es.begin(), it)); // 返回索引
	}
	else {
		return -1; // 未找到
	}
}

// 比较两个 Point2d 是否近似相等
inline bool arePointsEqual(const cv::Point2d& p1
	, const cv::Point2d& p2)
{
	return (std::abs(p1.x - p2.x) < EPSILON)
		&& (std::abs(p1.y - p2.y) < EPSILON);
}

// 查找 Point2d 在 vector<Point2d> 中的索引，未找到返回 -1
inline int findPointIndex(const std::vector<cv::Point2d>& points
	, const cv::Point2d& target)
{
	auto it = std::find_if(points.begin()
		, points.end(), [&](const cv::Point2d& p)
	{
		return arePointsEqual(p, target);
	});

	if (it != points.end()) {
		return std::distance(points.begin(), it);
	}

	return -1; // 未找到
}


inline Point2d findCentroid(Point2d A, Point2d B, Point2d C)
{
	Point2d centroid;
	centroid.x = (A.x + B.x + C.x) / 3.0;
	centroid.y = (A.y + B.y + C.y) / 3.0;
	return centroid;
}

inline Point2d findMidpoint(Point2d A, Point2d B)
{
	Point2d midpoint;
	midpoint.x = (A.x + B.x) / 2.0;
	midpoint.y = (A.y + B.y) / 2.0;
	return midpoint;
}

// 检查两条线段 AB 和 CD 是否共线，排除点重合情况
inline bool areCollinear(const cv::Point2d& A
	, const cv::Point2d& B, const cv::Point2d& C, const cv::Point2d& D)
{
	// 排除 A 和 B 重合，或者 C 和 D 重合的情况
	if (arePointsEqual(A, B) || arePointsEqual(C, D)) {
		return false; // 如果线段退化为一个点，不可能共线
	}

	// 计算叉积，判断是否接近零（表示共线）
	double cp1 = crossProduct(A, B, C);
	double cp2 = crossProduct(A, B, D);

	return (std::abs(cp1) < TOLERANCE) && (std::abs(cp2) < TOLERANCE);
}


// 计算仿射变换矩阵
inline Matrix3d computeAffineTransformation(const MatrixXd& A, const MatrixXd& B) {
	int n = A.cols();
	if (A.rows() != 2 || B.rows() != 2 || A.cols() != B.cols() || n < 2) {
		cerr << "Error: A and B must be 2xN matrices with N >= 2 points" << endl;
		exit(EXIT_FAILURE);
	}

	// 计算 A 和 B 的质心
	Vector2d centroid_A = A.rowwise().mean();
	Vector2d centroid_B = B.rowwise().mean();

	// 去中心化
	MatrixXd A_prime = A.colwise() - centroid_A;
	MatrixXd B_prime = B.colwise() - centroid_B;

	// 计算最优变换矩阵 H
	Matrix2d AAT_inv = A_prime * A_prime.transpose();

	// 检查是否接近零矩阵，避免奇异
	if (AAT_inv.norm() < 1e-8) {
		cerr << "Warning: Matrix is near singular, applying regularization" << endl;
		AAT_inv += Matrix2d::Identity() * 1e-6;
	}

	// 使用更稳定的 SVD 计算
	BDCSVD<Matrix2d> svd(AAT_inv, ComputeFullU | ComputeFullV);

	// 计算伪逆
	Matrix2d AAT_pinv = svd.solve(Matrix2d::Identity());

	// 计算 H 矩阵
	Matrix2d H = B_prime * A_prime.transpose() * AAT_pinv;

	Vector2d t = centroid_B - H * centroid_A;

	// 组装 3x3 变换矩阵
	Matrix3d T = Matrix3d::Identity();
	T.block<2, 2>(0, 0) = H;
	T.block<2, 1>(0, 2) = t;

	return T;
}

// 使用仿射变换矩阵 T 变换一个 2D 点
inline Vector2d applyAffineTransformation(const Matrix3d& T, const Vector2d& point) {
	Vector3d point_homogeneous(point(0), point(1), 1.0); // 转换为齐次坐标
	Vector3d transformed_point_homogeneous = T * point_homogeneous; // 变换
	return transformed_point_homogeneous.head<2>(); // 提取 x', y'
}



inline bool findDuplicateIndices(const cv::Point2d points[4], int& index1, int& index2)
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = i + 1; j < 4; ++j)
		{
			if (points[i] == points[j])
			{
				index1 = i;
				index2 = j;
				return true;
			}
		}
	}
	return false;
}

inline void split1p1edge(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(C);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(C);
	t2.push_back(A);
	addTs.push_back(t2);
}

inline void split1p1interior1(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearA, Point2d intersecNearB)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersec);
	t1.push_back(intersecNearB);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersecNearB);
	t2.push_back(intersec);
	t2.push_back(intersecNearA);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(A);
	t3.push_back(intersecNearA);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(B);
	t4.push_back(C);
	t4.push_back(intersec);
	addTs.push_back(t4);

	vector<Point2d> t5;
	t5.push_back(C);
	t5.push_back(A);
	t5.push_back(intersec);
	addTs.push_back(t5);
}

inline void split1p1interior2(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearAB, Point2d intersecNearBC)
{
	vector<Point2d> t1;
	t1.push_back(A);
	t1.push_back(intersecNearAB);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersecNearAB);
	t2.push_back(B);
	t2.push_back(intersec);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(B);
	t3.push_back(intersecNearBC);
	t3.push_back(intersec);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(intersec);
	t4.push_back(intersecNearBC);
	t4.push_back(C);
	addTs.push_back(t4);

	vector<Point2d> t5;
	t5.push_back(A);
	t5.push_back(intersec);
	t5.push_back(C);
	addTs.push_back(t5);
}

inline void split1p2interior2(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearBC)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersec);
	t1.push_back(A);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(B);
	t2.push_back(intersecNearBC);
	t2.push_back(intersec);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(intersecNearBC);
	t3.push_back(C);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(A);
	t4.push_back(intersec);
	t4.push_back(C);
	addTs.push_back(t4);
}

inline void split1p2interior1_ca(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearCA)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersec);
	t1.push_back(A);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(B);
	t2.push_back(C);
	t2.push_back(intersec);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(C);
	t3.push_back(intersecNearCA);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(intersec);
	t4.push_back(intersecNearCA);
	t4.push_back(A);
	addTs.push_back(t4);
}

inline void split1p2interior1_ab(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearAB)
{
	vector<Point2d> t1;
	t1.push_back(intersecNearAB);
	t1.push_back(intersec);
	t1.push_back(A);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(B);
	t2.push_back(intersec);
	t2.push_back(intersecNearAB);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(B);
	t3.push_back(C);
	t3.push_back(intersec);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(A);
	t4.push_back(intersec);
	t4.push_back(C);
	addTs.push_back(t4);
}

inline void splitEdge1(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(C);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(C);
	t2.push_back(A);
	addTs.push_back(t2);
}

inline void split1p2Edge2_BC(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C
	, Point2d intersec, Point2d intersecNearBC)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersecNearBC);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(intersecNearBC);
	t2.push_back(C);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(C);
	t3.push_back(A);
	addTs.push_back(t3);
}

inline void split1p2Edge2_CA(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C
	, Point2d intersec, Point2d intersecNearCA)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(C);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(C);
	t2.push_back(intersecNearCA);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(intersecNearCA);
	t3.push_back(A);
	addTs.push_back(t3);
}

inline void split2p2interior(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec)
{
	vector<Point2d> t1;
	t1.push_back(A);
	t1.push_back(B);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(B);
	t2.push_back(C);
	t2.push_back(intersec);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(A);
	t3.push_back(intersec);
	t3.push_back(C);
	addTs.push_back(t3);
}

inline void split2p2outer1_r(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersecNearAB1
	, Point2d intersecNearAB2, Point2d intersecNearBC1
	, Point2d intersecNearBC2)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersecNearBC2);
	t1.push_back(intersecNearAB2);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersecNearBC2);
	t2.push_back(intersecNearBC1);
	t2.push_back(intersecNearAB2);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersecNearAB2);
	t3.push_back(intersecNearBC1);
	t3.push_back(intersecNearAB1);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(intersecNearAB1);
	t4.push_back(intersecNearBC1);
	t4.push_back(C);
	addTs.push_back(t4);

	vector<Point2d> t5;
	t5.push_back(A);
	t5.push_back(intersecNearAB1);
	t5.push_back(C);
	addTs.push_back(t5);
}

inline void split2p2outer1_l(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersecNearAB1
	, Point2d intersecNearAB2, Point2d intersecNearBC1
	, Point2d intersecNearBC2)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersecNearBC2);
	t1.push_back(intersecNearAB2);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersecNearBC2);
	t2.push_back(intersecNearBC1);
	t2.push_back(intersecNearAB2);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersecNearAB2);
	t3.push_back(intersecNearBC1);
	t3.push_back(intersecNearAB1);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(intersecNearAB1);
	t4.push_back(intersecNearBC1);
	t4.push_back(A);
	addTs.push_back(t4);

	vector<Point2d> t5;
	t5.push_back(A);
	t5.push_back(intersecNearBC1);
	t5.push_back(C);
	addTs.push_back(t5);
}

inline void split2p2outer2(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersecNearAB1
	, Point2d intersecNearAB2, Point2d intersecNearCA
	, Point2d intersecNearBC)
{
	vector<Point2d> t1;
	t1.push_back(A);
	t1.push_back(intersecNearAB1);
	t1.push_back(intersecNearCA);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersecNearAB1);
	t2.push_back(intersecNearAB2);
	t2.push_back(intersecNearCA);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersecNearAB2);
	t3.push_back(B);
	t3.push_back(intersecNearBC);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(intersecNearAB2);
	t4.push_back(intersecNearBC);
	t4.push_back(intersecNearCA);
	addTs.push_back(t4);

	vector<Point2d> t5;
	t5.push_back(intersecNearCA);
	t5.push_back(intersecNearBC);
	t5.push_back(C);
	addTs.push_back(t5);
}


inline void split2p2edge_noVt(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearCA, Point2d intersecNearBC)
{
	vector<Point2d> t1;
	t1.push_back(A);
	t1.push_back(intersec);
	t1.push_back(intersecNearCA);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(B);
	t2.push_back(intersecNearBC);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersecNearCA);
	t3.push_back(intersec);
	t3.push_back(intersecNearBC);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(intersecNearCA);
	t4.push_back(intersecNearBC);
	t4.push_back(C);
	addTs.push_back(t4);
}

inline void split2p2edge_Vt1_bc(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearBC)
{
	vector<Point2d> t1;
	t1.push_back(intersec);
	t1.push_back(B);
	t1.push_back(intersecNearBC);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(intersecNearBC);
	t2.push_back(C);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(C);
	t3.push_back(A);
	addTs.push_back(t3);
}

inline void split2p2edge2_bc(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearB, Point2d intersecNearC)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersecNearB);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(intersecNearB);
	t2.push_back(intersecNearC);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(A);
	t3.push_back(intersec);
	t3.push_back(intersecNearC);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(A);
	t4.push_back(intersecNearC);
	t4.push_back(C);
	addTs.push_back(t4);
}

inline void split2p2edge2_ca(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearC, Point2d intersecNearA)
{
	vector<Point2d> t1;
	t1.push_back(intersec);
	t1.push_back(B);
	t1.push_back(intersecNearC);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(B);
	t2.push_back(C);
	t2.push_back(intersecNearC);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(A);
	t3.push_back(intersec);
	t3.push_back(intersecNearA);
	addTs.push_back(t3);

	vector<Point2d> t4;
	t4.push_back(intersecNearA);
	t4.push_back(intersec);
	t4.push_back(intersecNearC);
	addTs.push_back(t4);
}

inline void split2p2edge_Vt1_ca(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearCA)
{
	vector<Point2d> t1;
	t1.push_back(intersec);
	t1.push_back(B);
	t1.push_back(C);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(C);
	t2.push_back(intersecNearCA);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(intersecNearCA);
	t3.push_back(A);
	addTs.push_back(t3);
}

inline void split2p3_bc(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearBC)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(intersecNearBC);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersecNearBC);
	t2.push_back(C);
	t2.push_back(intersec);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(C);
	t3.push_back(A);
	addTs.push_back(t3);
}

inline void split2p3_ca(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d intersec
	, Point2d intersecNearCA)
{
	vector<Point2d> t1;
	t1.push_back(B);
	t1.push_back(C);
	t1.push_back(intersec);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(intersec);
	t2.push_back(C);
	t2.push_back(intersecNearCA);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(intersec);
	t3.push_back(intersecNearCA);
	t3.push_back(A);
	addTs.push_back(t3);
}

inline void split2(vector<vector<Point2d>> &addTs
	, Point2d A, Point2d B, Point2d C, Point2d ab
	, Point2d bc)
{
	vector<Point2d> t1;
	t1.push_back(ab);
	t1.push_back(B);
	t1.push_back(bc);
	addTs.push_back(t1);

	vector<Point2d> t2;
	t2.push_back(ab);
	t2.push_back(bc);
	t2.push_back(C);
	addTs.push_back(t2);

	vector<Point2d> t3;
	t3.push_back(A);
	t3.push_back(ab);
	t3.push_back(C);
	addTs.push_back(t3);
}


inline bool isAlmostEqual(const Point2d& p1, const Point2d& p2)
{
	return (std::fabs(p1.x - p2.x) < TOLERANCE)
		&& (std::fabs(p1.y - p2.y) < TOLERANCE);
}


// 判断点 (px, py) 是否在线段 (x1, y1) - (x2, y2) 上
inline bool isPointOnSegment(Point2d p, Point2d a, Point2d b) {
	double crossProduct = (p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x);
	if (fabs(crossProduct) > 1e-6) return false;  // 不共线
	double dotProduct = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);
	if (dotProduct < 0) return false;  // 超出 a 端点
	double squaredLengthBA = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
	if (dotProduct > squaredLengthBA) return false;  // 超出 b 端点
	return true;  // 在线段上
}

// **Ray-Casting 射线法** 判断点是否在闭合曲线内
inline bool isPointInsidePolygon(Point2d p, const vector<Point2d>& polygon) {
	int crossings = 0;
	size_t n = polygon.size();

	for (size_t i = 0; i < n; i++) {
		Point2d a = polygon[i];
		Point2d b = polygon[(i + 1) % n];

		// 点在边上，直接返回 true
		if (isPointOnSegment(p, a, b)) return true;

		// 保证 a.y <= b.y
		if (a.y > b.y) std::swap(a, b);

		// 忽略水平边
		if (abs(a.y - b.y) < 1e-10) continue;

		// 判断 p.y 是否在 (a.y, b.y) 之间（开区间）
		if (p.y > a.y && p.y <= b.y) {
			double xIntersect = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);

			if (xIntersect > p.x) {
				crossings++;
			}
		}
	}

	return (crossings % 2 == 1);
}

// 计算三角形归一化
inline void normalizeTriangle(const vector<Point2d>& T, vector<Point2d>& T_norm, Point2d& center, double& scale) {
	center.x = (T[0].x + T[1].x + T[2].x) / 3.0;
	center.y = (T[0].y + T[1].y + T[2].y) / 3.0;

	double maxDist = 0;
	for (const auto& pt : T) {
		double dist = norm(pt - center);
		if (dist > maxDist) maxDist = dist;
	}

	scale = (maxDist > 0) ? (1.0 / maxDist) : 1.0;

	T_norm.clear();
	for (const auto& pt : T) {
		T_norm.push_back((pt - center) * scale);
	}
}

// 计算点 A 在 T1 中的重心坐标
inline Vec3d computeBarycentricCoords(const Point2d& A, const vector<Point2d>& T1) {
	Matrix3d M;
	M << T1[0].x, T1[1].x, T1[2].x,
		T1[0].y, T1[1].y, T1[2].y,
		1.0, 1.0, 1.0;

	Vector3d rhs(A.x, A.y, 1.0);
	Vector3d lambda = M.colPivHouseholderQr().solve(rhs);
	return Vec3d(lambda(0), lambda(1), lambda(2));
}

// 计算对应点在 T2 中的位置
inline Point2d computeCorrespondingPoint(const Vec3d& lambda, const vector<Point2d>& T2) {
	return lambda[0] * T2[0] + lambda[1] * T2[1] + lambda[2] * T2[2];
}


// 传入 T1、T2 三角形的三个点，以及 T1 中的查询点 A
inline Point2d findCorrespondingPoint_affine(const Point2d& A,
	const std::vector<Point2d>& T1,
	const std::vector<Point2d>& T2) {
	Eigen::Vector2d P0(T1[0].x, T1[0].y);
	Eigen::Vector2d P1(T1[1].x, T1[1].y);
	Eigen::Vector2d P2(T1[2].x, T1[2].y);
	Eigen::Vector2d Q0(T2[0].x, T2[0].y);
	Eigen::Vector2d Q1(T2[1].x, T2[1].y);
	Eigen::Vector2d Q2(T2[2].x, T2[2].y);
	Eigen::Vector2d P(A.x, A.y);

	Eigen::Matrix2d A_mat;
	A_mat.col(0) = P1 - P0;
	A_mat.col(1) = P2 - P0;

	Eigen::Vector2d coeffs = A_mat.inverse() * (P - P0);

	Eigen::Vector2d result = Q0 + coeffs[0] * (Q1 - Q0) + coeffs[1] * (Q2 - Q0);
	return Point2d(result.x(), result.y());
}



// 计算两点间欧几里得距离
template<class T>
double __declspec(dllexport) dist(const T& P1, const T& P2) {
	return sqrt((P1.x - P2.x) * (P1.x - P2.x) + (P1.y - P2.y) * (P1.y - P2.y));
}

// 计算三角形的面积
template<class PointT>
double __declspec(dllexport) triangleArea(const std::vector<PointT>& tri)
{
	assert(tri.size() == 3);
	return 0.5 * fabs(
		tri[0].x * (tri[1].y - tri[2].y) +
		tri[1].x * (tri[2].y - tri[0].y) +
		tri[2].x * (tri[0].y - tri[1].y));
}

inline double triangleArea(const Point2d t1, const Point2d t2, const Point2d t3)
{
	return 0.5 * fabs(
		t1.x * (t2.y - t3.y) +
		t2.x * (t3.y - t1.y) +
		t3.x * (t1.y - t2.y));
}

template<class T>
double __declspec(dllexport) triangleArea4ori(const vector<T>& T) {
	return fabs((T[0].oci * (T[1].ori - T[2].ori) + T[1].oci
		* (T[2].ori - T[0].ori) + T[2].oci * (T[0].ori - T[1].ori)) / 2.0);
}

inline double clamp(double x, double a, double b) {
	return std::max(a, std::min(b, x));
}

inline double angleBetween(const Point2d& u, const Point2d& v)
{
	double dot = u.ddot(v);
	double norm_u = norm(u);
	double norm_v = norm(v);
	if (norm_u == 0 || norm_v == 0) return 0.0; // 避免除0
	double cos_theta = dot / (norm_u * norm_v);
	cos_theta = clamp(cos_theta, -1.0, 1.0); // 兼容 VS2013
	return acos(cos_theta);
}

inline double triangleMinAngle(const Point2d& A, const Point2d& B, const Point2d& C)
{
	double angleA = angleBetween(B - A, C - A);
	double angleB = angleBetween(A - B, C - B);
	double angleC = angleBetween(A - C, B - C);
	return std::min(std::min(angleA, angleB), angleC);
}

// 计算三角形的正三角形评分
template<class T>
double __declspec(dllexport) triangleQuality(const vector<T>& Tr) {
	if (Tr.size() != 3) {
		cerr << "Error: Triangle must have exactly 3 points!" << endl;
		return numeric_limits<double>::max();
	}

	double L1 = dist(Tr[0], Tr[1]);
	double L2 = dist(Tr[1], Tr[2]);
	double L3 = dist(Tr[2], Tr[0]);

	// 计算边长比误差
	double maxL = max({ L1, L2, L3 });
	double minL = min({ L1, L2, L3 });
	double E_ratio = maxL / minL - 1.0;

	// 计算角度误差
	double thetaA = acos((L2 * L2 + L3 * L3 - L1 * L1) / (2 * L2 * L3)) * 180.0 / PI;
	double thetaB = acos((L1 * L1 + L3 * L3 - L2 * L2) / (2 * L1 * L3)) * 180.0 / PI;
	double thetaC = 180.0 - thetaA - thetaB;

	double E_angle = (fabs(thetaA - 60) + fabs(thetaB - 60) + fabs(thetaC - 60)) / 3.0;

	// 计算高宽比
	double s = (L1 + L2 + L3) / 2.0;
	double area = triangleArea(Tr);
	double R = (L1 * L2 * L3) / (4.0 * area);
	double r = area / s;
	double E_aspect = fabs(R / r - 2.0);

	// 综合评分
	double w1 = 1.0, w2 = 1.5, w3 = 1.2;
	return w1 * E_ratio + w2 * (E_angle / 60) + w3 * E_aspect;
}

// 找到最接近正三角形的三角形
template<class T>
vector<T> __declspec(dllexport) findBestTriangle(const vector<vector<T>>& triangles)
{
	vector<T> bestTriangle;
	double minScore = numeric_limits<double>::max();

	for (const auto& S : triangles)
	{
		double score = triangleQuality(S);
		if (score < minScore)
		{
			minScore = score;
			bestTriangle = S;
		}
	}
	return bestTriangle;
}


// 计算三角形边长的平均值
template<class T>
double __declspec(dllexport) averageEdgeLength(const vector<T>& tri) {
	if (tri.size() != 3) {
		cerr << "Error: Triangle must have exactly 3 points!" << endl;
		return -1.0;  // 返回无效值
	}

	double L1 = dist(tri[0], tri[1]);
	double L2 = dist(tri[1], tri[2]);
	double L3 = dist(tri[2], tri[0]);

	return (L1 + L2 + L3) / 3.0;
}

template<class T>
int __declspec(dllexport) nearestAngle(const vector<T>& tri, T&v) {
	if (tri.size() != 3) {
		cerr << "Error: Triangle must have exactly 3 points!" << endl;
		return -1.0;  // 返回无效值
	}

	double L1 = dist(tri[0], v);
	double L2 = dist(tri[1], v);
	double L3 = dist(tri[2], v);

	if (L1 <= L2 && L1 <= L3)
	{
		return 0;
	}
	else if (L2 <= L1 && L2 <= L3)
	{
		return 1;
	}
	else if (L3 <= L1 && L3 <= L2)
	{
		return 2;
	}
	else
	{
		printf("nearestAngle return 3\n");
		return 3;
	}
}

// 提取唯一的点集
inline vector<Pt> extractUniqueVertices(const vector<vector<Pt>>& triangles)
{
	set<Pt, ComparePoints> uniquePoints;

	// 遍历所有三角形的顶点，插入到 set 中
	for (size_t i = 0; i < triangles.size(); i++) {
		for (size_t j = 0; j < triangles[i].size(); j++) {
			uniquePoints.insert(triangles[i][j]);
		}
	}

	// 转换 set 为 vector 输出
	return vector<Pt>(uniquePoints.begin(), uniquePoints.end());
}

// 从点集M中选出所有包含P的三角形
template<class T>
vector<vector<T>> __declspec(dllexport) findTrianglesContainingP(const vector<T>& M, const T& P)
{
	vector<vector<T>> result;

	int N = static_cast<int>(M.size());
	if (N < 3) return result;

	for (int i = 0; i < N - 2; ++i)
	{
		for (int j = i + 1; j < N - 1; ++j)
		{
			for (int k = j + 1; k < N; ++k)
			{
				const T &A = M[i], &B = M[j], &C = M[k];

				if (pointInTriangle(A, B, C, P))
				{
					result.push_back({ A, B, C });
				}
			}
		}
	}
	return result;
}

// 点到直线(AB)距离函数
inline double ptLineDist(const Point2d& P, const Point2d& A, const Point2d& B)
{
	double numerator = fabs((B.x - A.x)*(A.y - P.y) - (A.x - P.x)*(B.y - A.y));
	double denominator = sqrt(pow(B.x - A.x, 2) + pow(B.y - A.y, 2));
	return numerator / denominator;
}

inline double squaredDistance(const Point2d& a, const Point2d& b)
{
	return (a - b).dot(a - b);  // 推荐，更高效
	// 或者：return cv::norm(a - b) * cv::norm(a - b); // 不推荐，重复开根号
}

inline double pointToSegmentDistance(const Point2d& p
	, const Point2d& a, const Point2d& b) {
	double abx = b.x - a.x;
	double aby = b.y - a.y;
	double apx = p.x - a.x;
	double apy = p.y - a.y;

	double abLenSq = abx * abx + aby * aby;
	if (abLenSq == 0.0)
		return std::numeric_limits<double>::max(); // a 和 b 重合

	double t = (abx * apx + aby * apy) / abLenSq;

	if (t < 0.0 || t > 1.0) {
		return std::numeric_limits<double>::max(); // 垂足不在线段上
	}

	cv::Point2d projection(a.x + t * abx, a.y + t * aby);
	return std::sqrt(squaredDistance(p, projection));
}


// 坐标插值：fA, fB, fC 是坐标，而非标量值
inline Vector2d meanValueCoordInterp2D(
	const Vector2d& P,
	const Vector2d& A, const Vector2d& fA,
	const Vector2d& B, const Vector2d& fB,
	const Vector2d& C, const Vector2d& fC)
{
	Vector2d v0 = A - P;
	Vector2d v1 = B - P;
	Vector2d v2 = C - P;

	double d0 = v0.norm();
	double d1 = v1.norm();
	double d2 = v2.norm();

	// 防止除零
	const double eps = 1e-12;
	d0 = std::max(d0, eps);
	d1 = std::max(d1, eps);
	d2 = std::max(d2, eps);

	// 夹角计算（用夹角公式 + clamp 防止 acos 精度误差）
	auto safe_acos = [](double x) {
		return std::acos(std::max(-1.0, std::min(1.0, x)));
	};

	double theta0 = safe_acos(v0.dot(v1) / (d0 * d1));
	double theta1 = safe_acos(v1.dot(v2) / (d1 * d2));
	double theta2 = safe_acos(v2.dot(v0) / (d2 * d0));

	double w0 = (std::tan(theta2 / 2.0) + std::tan(theta0 / 2.0)) / d0;
	double w1 = (std::tan(theta0 / 2.0) + std::tan(theta1 / 2.0)) / d1;
	double w2 = (std::tan(theta1 / 2.0) + std::tan(theta2 / 2.0)) / d2;

	double w_sum = w0 + w1 + w2;
	w0 /= w_sum;
	w1 /= w_sum;
	w2 /= w_sum;

	// 坐标插值：f(P) = w0 * fA + w1 * fB + w2 * fC
	return w0 * fA + w1 * fB + w2 * fC;
}

// 计算点 P 在三角形 ABC 中的重心坐标
inline Vec3d computeBarycentricCoords(const Point2d& P,
	const Point2d& A,
	const Point2d& B,
	const Point2d& C) {
	// 使用面积法或向量法都可以
	Vector2d v0(B.x - A.x, B.y - A.y);
	Vector2d v1(C.x - A.x, C.y - A.y);
	Vector2d v2(P.x - A.x, P.y - A.y);

	double d00 = v0.dot(v0);
	double d01 = v0.dot(v1);
	double d11 = v1.dot(v1);
	double d20 = v2.dot(v0);
	double d21 = v2.dot(v1);

	double denom = d00 * d11 - d01 * d01;
	if (std::abs(denom) < 1e-10)
		return Vec3d(1.0 / 3, 1.0 / 3, 1.0 / 3); // 防退化

	double v = (d11 * d20 - d01 * d21) / denom;
	double w = (d00 * d21 - d01 * d20) / denom;
	double u = 1.0 - v - w;

	return Vec3d(u, v, w);
}

// 三角形重心插值版本
inline Point2d findCorrespondingPoint_Barycentric(const Point2d& A,
	const vector<Point2d>& T1,
	const vector<Point2d>& T2) {
	// 获取重心坐标
	Vec3d lambda = computeBarycentricCoords(A, T1[0], T1[1], T1[2]);

	// 用重心坐标插值 T2 上的点
	Point2d B = lambda[0] * T2[0] + lambda[1] * T2[1] + lambda[2] * T2[2];
	return B;
}

// 在两个相似三角形之间找到对应点
inline Point2d findCorrespondingPoint(const Point2d& A, const vector<Point2d>& T1
	, const vector<Point2d>& T2) {
	vector<Point2d> T1_norm, T2_norm;
	Point2d center1, center2;
	double scale1, scale2;

	normalizeTriangle(T1, T1_norm, center1, scale1);
	normalizeTriangle(T2, T2_norm, center2, scale2);

	Point2d A_norm = (A - center1) * scale1;

	Vector2d A_norm_v, T1_a, T1_b, T1_c, T2_a, T2_b, T2_c;
	A_norm_v(0) = A_norm.x;
	A_norm_v(1) = A_norm.y;
	T1_a(0) = T1_norm[0].x;
	T1_a(1) = T1_norm[0].y;

	T1_b(0) = T1_norm[1].x;
	T1_b(1) = T1_norm[1].y;

	T1_c(0) = T1_norm[2].x;
	T1_c(1) = T1_norm[2].y;

	T2_a(0) = T2_norm[0].x;
	T2_a(1) = T2_norm[0].y;

	T2_b(0) = T2_norm[1].x;
	T2_b(1) = T2_norm[1].y;

	T2_c(0) = T2_norm[2].x;
	T2_c(1) = T2_norm[2].y;

	Vector2d B_norm = meanValueCoordInterp2D(A_norm_v, T1_a
		, T2_a, T1_b, T2_b, T1_c, T2_c);

	return Point2d(B_norm(0) / scale2 + center2.x, B_norm(1) / scale2 + center2.y);
}

template<class T>
void __declspec(dllexport) GetWu1(T Wu1[], const T *matrixB
	, const T *P, const T *l, int N, int U)
{
	T *transposedB = new T[N*U];
	T *transposedBP = new T[N*U];
	MatrixTranspose(matrixB, transposedB, N, U);
	MultMatrix(transposedB, P, transposedBP, U, N, N);
	MultMatrix(transposedBP, l, Wu1, U, N, 1);

	delete[] transposedB;
	delete[] transposedBP;
}

//matrix[3][4] change to  matrix0[4][3]
template<class T, class T0>
void __declspec(dllexport) MatrixTranspose(T* matrix, T0* matrix0, int nrows, int ncols)
{// matrix - nrows*ncols    matrix0 - ncols*nrows
	for (int i = 0; i < nrows; ++i)
	{
		for (int j = 0; j < ncols; ++j)
		{
			matrix0[j*nrows + i] = matrix[i*ncols + j];
		}
	}
}


template<class TT1>
void __declspec(dllexport) GetNBB(TT1 nbb[], const TT1 *matrixB, const TT1 *P, int N, int U)
{
	TT1 *transposedB = new TT1[N*U];
	TT1 *transposedBP = new TT1[N*U];
	MatrixTranspose(matrixB, transposedB, N, U);
	MultMatrix(transposedB, P, transposedBP, U, N, N);
	MultMatrix(transposedBP, matrixB, nbb, U, N, U);

	delete[] transposedB;
	delete[] transposedBP;
}


/*!
weighted indirect adjustment model
\param correction
\param matrixB
\param l
\param P
\param N
\param U
*/
template<class T>
void __declspec(dllexport) GetCorrection(T correction[], const T *matrixB
	, const T *l, const T *P, int N, int U)
{
	T *nbb = new T[U*U];
	T *inverseForNbb = new T[U*U];
	T *Wu1 = new T[U];

	GetNBB<T>(nbb, matrixB, P, N, U);
	MatrixAnti(nbb, inverseForNbb, U);
	GetWu1<T>(Wu1, matrixB, P, l, N, U);

	MultMatrix(inverseForNbb, Wu1, correction, U, U, 1);

	delete[] nbb;
	delete[] inverseForNbb;
	delete[] Wu1;
}

template<class T>
void getL(T l[], const double a1, const double a2
	, vector<Point2d> &p0, vector<Point2d> &p, T delta_x, T delta_y)
{
	assert(p0.size() == p.size());
	int N = p0.size();
	memset(l, 0, sizeof(T)*N * 2);

	for (int i = 0; i != N; ++i)
	{
		l[i * 2 + 0] = -1 * a1*p0[i].x + a2*p0[i].y - delta_x + p[i].x;
		l[i * 2 + 1] = -1 * a2*p0[i].x - a1*p0[i].y - delta_y + p[i].y;
	}
}

template<class T>
void getMatrixB(T matrixB[], vector<Point2d> &ps)
{
	int N = ps.size();
	memset(matrixB, 0, sizeof(T)*N * 2 * 4);

	for (int i = 0; i != N; ++i)
	{
		matrixB[i * 8 + 0] = ps[i].x;
		matrixB[i * 8 + 1] = -1 * ps[i].y;
		matrixB[i * 8 + 2] = 1;
		matrixB[i * 8 + 3] = 0;

		matrixB[i * 8 + 4 + 0] = ps[i].y;
		matrixB[i * 8 + 4 + 1] = ps[i].x;
		matrixB[i * 8 + 4 + 2] = 0;
		matrixB[i * 8 + 4 + 3] = 1;
	}
}

// 计算所有样本对之间的距离平方，并返回中位数
template<class T>
T __declspec(dllexport) compute_median_distance(const vector<Point2d>& data) {
	std::vector<T> distances;

	for (size_t i = 0; i < data.size(); ++i) {
		for (size_t j = i + 1; j < data.size(); ++j) {
			distances.push_back(norm(data[i] - data[j]));
		}
	}

	if (distances.empty())
		return 0.0;

	size_t n = distances.size();
	std::nth_element(distances.begin(), distances.begin() + n / 2, distances.end());
	T median;

	if (n % 2 == 0) {
		T a = *std::max_element(distances.begin(), distances.begin() + n / 2);
		T b = *std::min_element(distances.begin() + n / 2, distances.end());
		median = (a + b) / 2.0;
	}
	else {
		median = distances[n / 2];
	}

	return std::sqrt(median / 2.0);  // sigma
}

// 高斯核函数：计算 x 对中心点 c 的核权重
template<class T>
T __declspec(dllexport) gaussian_kernel(const Point2d& x, const Point2d& c,
	T sigma) {
	T dist2 = norm(x - c);
	return std::exp(-dist2 / (2.0 * sigma * sigma));
}


template<class T>
void __declspec(dllexport) getP4Curv(T *P, vector<Point2d> &ps, Point2d curr_point)
{
	T dist = 0;
	T w = 0;

	int num = ps.size();

	T sigma = compute_median_distance<T>(ps);

	for (int i = 0; i != num * 2; ++i)
	{
		for (int j = 0; j != num * 2; ++j)
		{
			if (i != j)
			{
				P[i*num * 2 + j] = 0;
			}
			else
			{
				w = gaussian_kernel(ps[i / 2], curr_point, sigma);
				P[i*num * 2 + j] = w;
			}
		}
	}
}

template<class T>
void __declspec(dllexport) getDeltaVector(T DeltaV[], const T a1, const T a2
	, vector<Point2d> &p0, vector<Point2d> &p, Point2d curr_point, T delta_x, T delta_y)
{
	assert(p0.size() == p.size());
	int N = p0.size();
	T *matrixB = new T[N * 2 * 4];
	T *l = new T[N * 2];
	T P[4096] = { 0 };
	getP4Curv(P, p0, curr_point);
	getMatrixB(matrixB, p0);
	getL(l, a1, a2, p0, p, delta_x, delta_y);
	GetCorrection(DeltaV, matrixB, l, P, N * 2, 4);

	delete[] matrixB;
	delete[] l;
}


inline void GetTransformParams(double &a1, double &a2, double &tran_x, double &tran_y
	, vector<Point2d> &p0, vector<Point2d> &p, Point2d current_pt)
{
	a1 = 0;
	a2 = 0;
	tran_x = 0;
	tran_y = 0;
	double deltaV[4] = { 0 };
	assert(p0.size() == p.size());

	do
	{
		getDeltaVector(deltaV, a1, a2, p0, p, current_pt, tran_x, tran_y);

		a1 += deltaV[0];
		a2 += deltaV[1];

		tran_x += deltaV[2];
		tran_y += deltaV[3];
	} while (sqrt(pow(deltaV[0], 2) + pow(deltaV[1], 2)
		+ pow(deltaV[2], 2) + pow(deltaV[3], 2)) < 0.00001);
}

inline Point2d applyAffineTransform(
	const Matx22d& A,
	const Point2d& t,
	const Point2d& P)
{
	// 矩阵乘以点 + 平移
	double x_new = A(0, 0) * P.x + A(0, 1) * P.y + t.x;
	double y_new = A(1, 0) * P.x + A(1, 1) * P.y + t.y;
	return Point2d(x_new, y_new);
}

// 输入两个三角形中 A→B 和 A′→B′ 的向量，输出旋转角（单位：弧度，范围：0~2π）
inline double estimateTriangleRotationRadian(
	const cv::Point2d& A, const cv::Point2d& B,
	const cv::Point2d& A_, const cv::Point2d& B_)
{
	cv::Point2d v1 = B - A;
	cv::Point2d v2 = B_ - A_;

	double angle1 = std::atan2(v1.y, v1.x);
	double angle2 = std::atan2(v2.y, v2.x);

	double angle_rad = angle2 - angle1;

	// 映射到 [0, 2π)
	if (angle_rad < 0)
		angle_rad += 2 * CV_PI;

	return angle_rad; // 单位：弧度
}

inline FootResult repelFromSegmentIfTooClose(
	const cv::Point2d& p,
	const cv::Point2d& a,
	const cv::Point2d& b,
	double min_dist)
{
	cv::Point2d ab = b - a;
	double ab_len_sq = ab.dot(ab);
	if (ab_len_sq == 0.0) {
		return{ p, false }; // a 和 b 重合，不处理
	}

	double ab_len = std::sqrt(ab_len_sq);
	cv::Point2d ab_unit;
	ab_unit.x = ab.x / ab_len;
	ab_unit.y = ab.y / ab_len;

	// 求垂足参数 t，判断是否在段上
	cv::Point2d ap = p - a;
	double t = ab.dot(ap) / ab_len_sq;
	bool onSegment = (t >= 0.0 && t <= 1.0);
	if (!onSegment) {
		return{ p, false };
	}

	// 计算垂足
	cv::Point2d foot = a + t * ab;

	// 距离过远，无需调整
	cv::Point2d vec_foot_to_p = p - foot;
	double dist_sq = vec_foot_to_p.dot(vec_foot_to_p);
	if (dist_sq >= min_dist * min_dist)
	{
		return{ p, false };
	}

	// 法向方向（左手法线）
	cv::Point2d normal(-ab_unit.y, ab_unit.x);

	// 判断 p 在 ab 的哪一侧（使用叉积）
	double cross = ab.x * ap.y - ab.y * ap.x;
	cv::Point2d direction = (cross >= 0.0) ? normal : -normal;

	// 沿远离 ab 的方向移动 (min_dist - dist)
	double dist = std::sqrt(dist_sq);
	double move_len = min_dist - dist;
	cv::Point2d adjustedP = p + direction * move_len;

	return{ adjustedP, true };
}

inline int func_nc8(int *b)
/* connectivity detection for each point */
{
	int n_odd[4] = { 1, 3, 5, 7 };  /* odd-number neighbors */
	int i, j, sum, d[10];           /* control variable */

	for (i = 0; i <= 9; i++)
	{
		j = i;
		if (i == 9)
			j = 1;
		if (abs(*(b + j)) == 1)
		{
			d[i] = 1;
		}
		else
		{
			d[i] = 0;
		}
	}
	sum = 0;
	for (i = 0; i < 4; i++)
	{
		j = n_odd[i];
		sum = sum + d[j] - d[j] * d[j + 1] * d[j + 2];
	}
	return (sum);
}

inline bool isEquilateralTriangle(const Point2d& A, const Point2d& B, const Point2d& C
	, double tol = 1e-6)
{
	double ab = norm(B - A);
	double bc = norm(C - B);
	double ca = norm(A - C);

	return (fabs(ab - bc) < tol) && (fabs(bc - ca) < tol) && (fabs(ca - ab) < tol);
}
#endif