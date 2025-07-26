#pragma once
/*!
*      \class ImgProc
*      \brief image analysis
*	   \author Wei Feng
*      \date 11/21/2024
*/
#include <opencv2/opencv.hpp>
#include "Algo.h"
#include "KdTree\kdtree.h"
#define MAX_VERTEX_NUM 1048576
#define MAX_TRIANGLE_NUM 524288
#define MINV 0.000001
#define SEARCHITEMNUM 16
#define BLACK 0
#define WHITE 255
#define GRAY 128
typedef double double2[2];
typedef int int3[3];

enum TOPOTYPE
{
	INNER, OUTER
};

struct LineSegIntersec
{
	vector<Point2d> ps;
	Point2d p1;
	Point2d p2;
	int res1;
	int res2;
	int mask;
};

struct StatisPara
{
	double minAngle;
	double minArea;
	int sliverNum;
	double areaVari;
	int equilateralCount;
};

class ImgProc
{
	friend class DataAnalysis;
public:
	ImgProc(const char *maskPath, const char *bdyPath, const char *oriPath);
	virtual ~ImgProc();
	void Analysis();
	void AvailModel4Output(double2 *ps_op, int &pN
		, uchar *gs_op, int &gN, int3 *tri_set, int &tN);
	void OuterModel4Output(double2 *ps_op, int &pN
		, uchar *gs_op, int &gN, int3 *tri_set, int &tN);
	void GeoQualEval4Avail(StatisPara &statisPara);
	void OutputPsBdy(double2 *innPs, int &pN, int *bdy, int &bN);
	void OutputBdy(double2 *ps, int *idx, int &pN);
	void OutputTris(int3 *tri_set, int &tN);
	void OutputBdyIds(int *ids, int &bN);
private:
	Mat ori_mat;
	Mat src_mat;
	Mat bdy_mat;
	Mat sobel_mat;
	Mat canny_mat;
	Mat mask_mat;
	int r;
	int c;
	vector<Point2d> ps;
	vector<uchar> gs;
	vector<int> boundary;
	vector<vector<int>> triangle_set;
	vector<vector<int>> avail_tria;
	vector<vector<int>> outer_tria;
	vector<Point2d> orderedBdy;

	vector<vector<Point2d>> ts_1;
	vector<vector<Point2d>> ts_1p1;
	vector<vector<Point2d>> ts_1p2;
	vector<vector<Point2d>> ts_2p2;
	vector<vector<Point2d>> ts_2p3;
	vector<vector<Point2d>> ts_1p3;

	// kd tree of model points
	kdtree::KDTree*     M_tree4Tri;
	kdtree::KDTreeArray M_data4Tri;

	kdtree::KDTree*     M_tree4P;
	kdtree::KDTreeArray M_data4P;

	kdtree::KDTree*     M_tree4Padd;
	kdtree::KDTreeArray M_data4Padd;
	double EuclidDist(double2 &v1, double2 &v2);
	void Sobel4Edge(Mat &sMat, Mat &srcMat);
	void GetBdy();
	uchar GetinterpoGray(Mat &mat, Point2d offset, int i, int j);
	void replaceBdy(vector<Point2d> &ps, vector<int> &untreatedIds, double thres);
	void replacePs(vector<Point2d> &ps, double thres);
	void DrawBdy(const char *filePath);
	bool GetFirstBdyPt(Mat &mat);
	void TriClassify();
	void ClipTri(vector<int> &intersectTridxSet
		, vector<vector<LineSegIntersec>> &lSegIntersecSet
		, vector<Point2d> &addPs, vector<Point2d> &edgeCtrs
		, vector<Point2d> &triCtrs, vector<int> &untreatedIds);
	void HandleNearEdgePt(vector<vector<Point2d>> &addTriS, vector<Point2d> &addPs
		, vector<int> &delTriS, int idx, double thres);
	void InitKdTree4Tri(vector<Point2d> &triCtrs);
	void Convert2IdxS(vector<vector<int>> &idxSet
		, vector<Point2d> &addPs, vector<vector<Point2d>> &addTs);
	void GetBdyCts(vector<Point2d> &edgeCtrs);
	void GetTriCts(vector<Point2d> &triCtrs);
	void Triangula4Edge(vector<LineSegIntersec> &lSegIntersec, int &tridx);

	void Triangula4Edge2(vector<LineSegIntersec> &lSegIntersec
		, Point2d A, Point2d B, Point2d C);

	void Triangula4Edge1(vector<LineSegIntersec> &lSegIntersec
		, Point2d A, Point2d B, Point2d C);

	int Triangula4Edge2_1p1(vector<vector<Point2d>> &ts_1p1
		, vector<LineSegIntersec> &lSegIntersec
		, Point2d A, Point2d B, Point2d C);

	int Triangula4Edge2_1p2(vector<vector<Point2d>> &addTs
		, vector<LineSegIntersec> &lSegIntersec
		, Point2d A, Point2d B, Point2d C);

	int Triangula4Edge2_2p2(vector<vector<Point2d>> &addTs
		, vector<LineSegIntersec> &lSegIntersec
		, Point2d A, Point2d B, Point2d C);

	int Triangula4Edge2_1p3(vector<vector<Point2d>> &addTs
		, vector<LineSegIntersec> &lSegIntersec
		, Point2d A, Point2d B, Point2d C);

	int Triangula4Edge2_2p3(vector<vector<Point2d>> &addTs
		, vector<LineSegIntersec> &lSegIntersec
		, Point2d A, Point2d B, Point2d C);

	int GetIntersecParams(Point2d &intersec, int &rest1, int &rest2
		, vector<LineSegIntersec> &lSegIntersec);

	void GetAdjacFarPs(vector<Point2d> &adjacentPoints, vector<Point2d> &farSet
		, vector<LineSegIntersec> &lSegIntersec, Point2d intersec);

	bool IsVertexCoincident(Point2d p, Point2d A, Point2d B, Point2d C);

	void DrawOrigMesh(Mat &mat, double zoom);
	void DrawTriMesh(Mat &mat, double zoom);
};

