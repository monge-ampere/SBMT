#include "stdafx.h"
#include "ImgProc.h"


ImgProc::ImgProc(const char *maskPath, const char *bdyPath, const char *oriPath)
{
	src_mat = imread(oriPath, IMREAD_GRAYSCALE);
	bdy_mat = imread(bdyPath, IMREAD_GRAYSCALE);
	mask_mat = imread(maskPath, IMREAD_GRAYSCALE);
	assert(src_mat.rows == mask_mat.rows);
	assert(src_mat.cols == mask_mat.cols);
	r = src_mat.rows;
	c = src_mat.cols;
	canny_mat.create(r, c, CV_8UC1);

	//ori_mat.create(r, c, CV_8UC3);

	//for (int i = 0; i != r; ++i)
	//{
	//	for (int j = 0; j != c; ++j)
	//	{
	//		uchar d = src_mat.at<uchar>(i, j);
	//		ori_mat.at<Vec3b>(i, j) = Vec3b(d, d, d);
	//	}
	//}
	ori_mat = imread(oriPath, IMREAD_COLOR);
	imwrite("F:\\gray.bmp", src_mat);
	Sobel4Edge(sobel_mat, src_mat);
	imwrite("F:\\sobel.bmp", sobel_mat);
	double lowThreshold = 80;
	double highThreshold = 120;
	Canny(src_mat, canny_mat, lowThreshold, highThreshold);

	imwrite("F:\\canny.bmp", canny_mat);

	GetBdy();
	getReguTris(ps, gs, triangle_set, src_mat);
	Point2d s1, s2, s3;
	s1.x = 163.392027;
	s1.y = 222.392027;
	s2.x = 163.391515;
	s2.y = 222.391515;
	s3.x = 163.038366;
	s3.y = 221.779842;
	for (int i = 0; i != triangle_set.size(); ++i)
	{
		Point2d A, B, C;
		A = ps[triangle_set[i][0]];
		B = ps[triangle_set[i][1]];
		C = ps[triangle_set[i][2]];

		if (norm(A - s1) < TOLERANCE)
		{
			printf("NO %d triangle's A is s1\n", i);
		}
		if (norm(A - s2) < TOLERANCE)
		{
			printf("NO %d triangle's A is s2\n", i);
		}
		if (norm(A - s3) < TOLERANCE)
		{
			printf("NO %d triangle's A is s3\n", i);
		}

		if (norm(B - s1) < TOLERANCE)
		{
			printf("NO %d triangle's B is s1\n", i);
		}
		if (norm(B - s2) < TOLERANCE)
		{
			printf("NO %d triangle's B is s2\n", i);
		}
		if (norm(B - s3) < TOLERANCE)
		{
			printf("NO %d triangle's B is s3\n", i);
		}

		if (norm(C - s1) < TOLERANCE)
		{
			printf("NO %d triangle's C is s1\n", i);
		}
		if (norm(C - s2) < TOLERANCE)
		{
			printf("NO %d triangle's C is s2\n", i);
		}
		if (norm(C - s3) < TOLERANCE)
		{
			printf("NO %d triangle's C is s3\n", i);
		}
	}

	for (int i = 0; i != orderedBdy.size(); ++i)
	{
		int oI = (orderedBdy.size() + i - 1) % orderedBdy.size();

		if (isPointOnSegment(s1, orderedBdy[oI], orderedBdy[i])
			&& isPointOnSegment(s2, orderedBdy[oI], orderedBdy[i]))
		{
			printf("s1s2 is on NO %d edge\n", i);
		}
	}
	Mat show_mat;
	show_mat.create(4096, 4096, CV_8UC3);
	show_mat.setTo(cv::Scalar(255, 255, 255));

	DrawOrigMesh(show_mat, 40);
	DrawTriMesh(show_mat, 40);
	imwrite("F:\\show_mat.bmp", show_mat);
}


ImgProc::~ImgProc()
{
}

bool ImgProc::GetFirstBdyPt(Mat &mat)
{
	for (int i = 0; i != r; ++i)
	{
		for (int j = 0; j != c; ++j)
		{
			uchar data = mat.at<uchar>(i, j);

			if (255 == data)
			{
				boundary.push_back(i*c + j);
				Point2d p;
				p.x = j + 0.5;
				p.y = i + 0.5;
				orderedBdy.push_back(p);
				mat.at<uchar>(i, j) = 0;
				return true;
			}
		}
	}

	return false;
}

void ImgProc::GetBdy()
{
	Mat bdy_cp = bdy_mat;
	GetFirstBdyPt(bdy_cp);
	int i, j;
	uchar d1, d2, d3, d4;
	Point2d p = orderedBdy[0];
	while (1)
	{
		i = p.y - 0.5;
		j = p.x - 0.5;

		d1 = bdy_cp.at<uchar>(i - 1, j);
		d2 = bdy_cp.at<uchar>(i + 1, j);
		d3 = bdy_cp.at<uchar>(i, j - 1);
		d4 = bdy_cp.at<uchar>(i, j + 1);

		if (255 == d1)
		{
			boundary.push_back((i - 1)*c + j);
			p.x = j + 0.5;
			p.y = i - 1 + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i - 1, j) = 0;
			continue;
		}
		else if (255 == d2)
		{
			boundary.push_back((i + 1)*c + j);
			p.x = j + 0.5;
			p.y = i + 1 + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i + 1, j) = 0;
			continue;
		}
		else if (255 == d3)
		{
			boundary.push_back(i*c + j - 1);
			p.x = j - 1 + 0.5;
			p.y = i + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i, j - 1) = 0;
			continue;
		}
		else if (255 == d4)
		{
			boundary.push_back(i*c + j + 1);
			p.x = j + 1 + 0.5;
			p.y = i + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i, j + 1) = 0;
			continue;
		}

		d1 = bdy_cp.at<uchar>(i - 1, j - 1);
		d2 = bdy_cp.at<uchar>(i - 1, j + 1);
		d3 = bdy_cp.at<uchar>(i + 1, j - 1);
		d4 = bdy_cp.at<uchar>(i + 1, j + 1);

		if (255 == d1)
		{
			boundary.push_back((i - 1)*c + j - 1);
			p.x = j - 1 + 0.5;
			p.y = i - 1 + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i - 1, j - 1) = 0;
			continue;
		}
		else if (255 == d2)
		{
			boundary.push_back((i - 1)*c + j + 1);
			p.x = j + 1 + 0.5;
			p.y = i - 1 + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i - 1, j + 1) = 0;
			continue;
		}
		else if (255 == d3)
		{
			boundary.push_back((i + 1)*c + j - 1);
			p.x = j - 1 + 0.5;
			p.y = i + 1 + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i + 1, j - 1) = 0;
			continue;
		}
		else if (255 == d4)
		{
			boundary.push_back((i + 1)*c + j + 1);
			p.x = j + 1 + 0.5;
			p.y = i + 1 + 0.5;
			orderedBdy.push_back(p);
			bdy_cp.at<uchar>(i + 1, j + 1) = 0;
			continue;
		}
		else
		{
			break;
		}
	}
}

void ImgProc::GetBdyCts(vector<Point2d> &edgeCtrs)
{
	int bN = orderedBdy.size();
	int pI;
	Point2d ct;
	for (int i = 0; i != bN; ++i)
	{
		pI = (bN + i - 1) % bN;
		ct = findMidpoint(orderedBdy[pI], orderedBdy[i]);
		edgeCtrs.push_back(ct);
	}
}

void ImgProc::GetTriCts(vector<Point2d> &triCtrs)
{
	int tN = triangle_set.size();
	Point2d ct, A, B, C;

	for (int i = 0; i != tN; ++i)
	{
		A = ps[triangle_set[i][0]];
		B = ps[triangle_set[i][1]];
		C = ps[triangle_set[i][2]];
		ct = findCentroid(A, B, C);
		triCtrs.push_back(ct);
	}
}

void ImgProc::Analysis()
{
	if (avail_tria.size()
		&& outer_tria.size())
	{
		return;
	}


	vector<int> untreatedIds;
	replaceBdy(ps, untreatedIds, 0.26);

	for (int i = 0; i != untreatedIds.size(); ++i)
	{
		ps.push_back(orderedBdy[untreatedIds[i]]);

		uchar data = src_mat.at<uchar>(
			floor(orderedBdy[untreatedIds[i]].y)
			, floor(orderedBdy[untreatedIds[i]].x));
		gs.push_back(data);
	}

	replacePs(ps, 0.183);

	vector<Point2d> addPs;
	vector<vector<Point2d>> addTs;
	vector<Point2d> edgeCtrs;
	vector<Point2d> triCtrs;
	vector<int> delTriS;
	vector<vector<int> > copy_triangle_set;
	vector<vector<int>> idxSet;
	for (int i = 0; i != orderedBdy.size(); ++i)
	{
		GetTriCts(triCtrs);
		InitKdTree4Tri(triCtrs);

		HandleNearEdgePt(addTs, addPs, delTriS, i, 0.125);

		for (int i = 0; i != triangle_set.size(); ++i)
		{
			if (-1 != isExi(delTriS, i))
			{
				continue;
			}
			copy_triangle_set.push_back(triangle_set[i]);
		}

		triangle_set = copy_triangle_set;


		Convert2IdxS(idxSet, addPs, addTs);

		for (int i = 0; i != addTs.size(); ++i)
		{
			triangle_set.push_back(idxSet[i]);
		}

		addPs.clear();
		addTs.clear();
		delTriS.clear();
		triCtrs.clear();
		copy_triangle_set.clear();
		idxSet.clear();
		delete M_tree4Tri;

		//printf("No %d is finished for handling the near edge point\n\n", i);
	}

	GetBdyCts(edgeCtrs);
	GetTriCts(triCtrs);

	vector<int> intersectTridxSet;
	vector<vector<LineSegIntersec>> lSegIntersecSet;
	ClipTri(intersectTridxSet, lSegIntersecSet
		, addPs, edgeCtrs, triCtrs, untreatedIds);

	assert(intersectTridxSet.size() == lSegIntersecSet.size());

	int SegIntersecNum = lSegIntersecSet.size();
	for (int i = 0; i != SegIntersecNum; ++i)
	{
		Triangula4Edge(lSegIntersecSet[i], intersectTridxSet[i]);
	}

	for (int i = 0; i != ts_1p1.size(); ++i)
	{
		addTs.push_back(ts_1p1[i]);
	}
	for (int i = 0; i != ts_1p2.size(); ++i)
	{
		addTs.push_back(ts_1p2[i]);
	}
	for (int i = 0; i != ts_2p2.size(); ++i)
	{
		addTs.push_back(ts_2p2[i]);
	}
	for (int i = 0; i != ts_1p3.size(); ++i)
	{
		addTs.push_back(ts_1p3[i]);
	}
	for (int i = 0; i != ts_2p3.size(); ++i)
	{
		addTs.push_back(ts_2p3[i]);
	}

	for (int i = 0; i != ts_1.size(); ++i)
	{
		addTs.push_back(ts_1[i]);
	}
	for (int i = 0; i != triangle_set.size(); ++i)
	{
		if (-1 != isExi(intersectTridxSet, i))
		{
			continue;
		}
		copy_triangle_set.push_back(triangle_set[i]);
	}

	triangle_set = copy_triangle_set;
	Convert2IdxS(idxSet, addPs, addTs);

	for (int i = 0; i != addTs.size(); ++i)
	{
		triangle_set.push_back(idxSet[i]);
	}

	Mat showMat;
	showMat.create(src_mat.rows * 15, src_mat.cols * 15, CV_8UC3);
	showMat.setTo(cv::Scalar(255, 255, 255));
	cv::Scalar color1(255, 0, 0);
	int thickness1 = 1;

	cv::Scalar color2(0, 0, 255);
	int thickness2 = 1;

	cv::Scalar color3(0, 0, 0);
	int thickness3 = 1;

	Point2d At, Bt, Ct;

	for (int i = src_mat.rows / 2; i != src_mat.rows - 1; ++i)
	{
		for (int j = 1/*src_mat.cols / 2*/; j != src_mat.cols - 1; ++j)
		{
			int ii, jj;
			ii = i - src_mat.rows / 2;
			jj = j /*- src_mat.cols / 2*/;
			cv::Point startPoint1((jj - 1) * 45, ii * 45);
			cv::Point endPoint1(jj * 45, ii * 45);

			if ((jj - 1) * 45 < src_mat.cols * 15
				&& jj * 45 < src_mat.cols * 15
				&& ii * 45 < src_mat.rows * 15)
			{
				line(showMat, startPoint1, endPoint1, color3, thickness3, cv::LINE_AA);
			}

			cv::Point startPoint2((jj + 1) * 45, ii * 45);
			cv::Point endPoint2(jj * 45, ii * 45);

			if ((jj + 1) * 45 < src_mat.cols * 15
				&& jj * 45 < src_mat.cols * 15
				&& ii * 45 < src_mat.rows * 15)
			{
				line(showMat, startPoint2, endPoint2, color3, thickness3, cv::LINE_AA);
			}

			cv::Point startPoint3(jj * 45, (ii - 1) * 45);
			cv::Point endPoint3(jj * 45, ii * 45);

			if ((ii - 1) * 45 < src_mat.rows * 15
				&& jj * 45 < src_mat.cols * 15
				&& ii * 45 < src_mat.rows * 15)
			{
				line(showMat, startPoint3, endPoint3, color3, thickness3, cv::LINE_AA);
			}

			cv::Point startPoint4(jj * 45, (ii + 1) * 45);
			cv::Point endPoint4(jj * 45, ii * 45);

			if ((ii + 1) * 45 < src_mat.rows * 15
				&& jj * 45 < src_mat.cols * 15
				&& ii * 45 < src_mat.rows * 15)
			{
				line(showMat, startPoint4, endPoint4, color3, thickness3, cv::LINE_AA);
			}
		}
	}

	for (int i = 0; i != triangle_set.size(); ++i)
	{
		Point2d A = ps[triangle_set[i][0]];
		Point2d B = ps[triangle_set[i][1]];
		Point2d C = ps[triangle_set[i][2]];

		At.x = (A.x /*- src_mat.cols / 2*/) * 45;
		At.y = (A.y - src_mat.rows / 2) * 45;

		Bt.x = (B.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (B.y - src_mat.rows / 2) * 45;

		Ct.x = (C.x /*- src_mat.cols / 2*/) * 45;
		Ct.y = (C.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x>src_mat.cols * 15
			|| At.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x>src_mat.cols * 15
			|| Bt.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Ct.x < 0
			|| Ct.y < 0
			|| Ct.x>src_mat.cols * 15
			|| Ct.y>src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		line(showMat, startPoint1, endPoint1, color1, thickness1, cv::LINE_AA);

		cv::Point startPoint2(Bt.x, Bt.y);
		cv::Point endPoint2(Ct.x, Ct.y);

		line(showMat, startPoint2, endPoint2, color1, thickness1, cv::LINE_AA);

		cv::Point startPoint3(Ct.x, Ct.y);
		cv::Point endPoint3(At.x, At.y);

		line(showMat, startPoint3, endPoint3, color1, thickness1, cv::LINE_AA);
	}

	for (int i = 0; i != ts_1.size(); ++i)
	{
		Point2d A = ts_1[i][0];
		Point2d B = ts_1[i][1];
		Point2d C = ts_1[i][2];

		At.x = (A.x /*- src_mat.cols / 2*/) * 45;
		At.y = (A.y - src_mat.rows / 2) * 45;

		Bt.x = (B.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (B.y - src_mat.rows / 2) * 45;

		Ct.x = (C.x /*- src_mat.cols / 2*/) * 45;
		Ct.y = (C.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x > src_mat.cols * 15
			|| At.y > src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x > src_mat.cols * 15
			|| Bt.y > src_mat.rows * 15)
		{
			continue;
		}

		if (Ct.x < 0
			|| Ct.y < 0
			|| Ct.x > src_mat.cols * 15
			|| Ct.y > src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		cv::Scalar color(0, 255, 0);
		int thickness = 1;

		line(showMat, startPoint1, endPoint1, color, thickness, cv::LINE_AA);

		cv::Point startPoint2(Bt.x, Bt.y);
		cv::Point endPoint2(Ct.x, Ct.y);

		line(showMat, startPoint2, endPoint2, color, thickness, cv::LINE_AA);

		cv::Point startPoint3(Ct.x, Ct.y);
		cv::Point endPoint3(At.x, At.y);

		line(showMat, startPoint3, endPoint3, color, thickness, cv::LINE_AA);

	}

	for (int i = 0; i != ts_1p1.size(); ++i)
	{
		Point2d A = ts_1p1[i][0];
		Point2d B = ts_1p1[i][1];
		Point2d C = ts_1p1[i][2];

		At.x = (A.x /*- src_mat.cols / 2*/) * 45;
		At.y = (A.y - src_mat.rows / 2) * 45;

		Bt.x = (B.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (B.y - src_mat.rows / 2) * 45;

		Ct.x = (C.x /*- src_mat.cols / 2*/) * 45;
		Ct.y = (C.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x>src_mat.cols * 15
			|| At.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x>src_mat.cols * 15
			|| Bt.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Ct.x < 0
			|| Ct.y < 0
			|| Ct.x>src_mat.cols * 15
			|| Ct.y>src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		cv::Scalar color(203, 192, 255);
		int thickness = 1;

		line(showMat, startPoint1, endPoint1, color, thickness, cv::LINE_AA);

		cv::Point startPoint2(Bt.x, Bt.y);
		cv::Point endPoint2(Ct.x, Ct.y);

		line(showMat, startPoint2, endPoint2, color, thickness, cv::LINE_AA);

		cv::Point startPoint3(Ct.x, Ct.y);
		cv::Point endPoint3(At.x, At.y);

		line(showMat, startPoint3, endPoint3, color, thickness, cv::LINE_AA);
	}

	for (int i = 0; i != ts_2p2.size(); ++i)
	{
		Point2d A = ts_2p2[i][0];
		Point2d B = ts_2p2[i][1];
		Point2d C = ts_2p2[i][2];

		At.x = (A.x /*- src_mat.cols / 2*/) * 45;
		At.y = (A.y - src_mat.rows / 2) * 45;

		Bt.x = (B.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (B.y - src_mat.rows / 2) * 45;

		Ct.x = (C.x /*- src_mat.cols / 2*/) * 45;
		Ct.y = (C.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x>src_mat.cols * 15
			|| At.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x>src_mat.cols * 15
			|| Bt.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Ct.x < 0
			|| Ct.y < 0
			|| Ct.x>src_mat.cols * 15
			|| Ct.y>src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		cv::Scalar color(128, 0, 128);
		int thickness = 1;

		line(showMat, startPoint1, endPoint1, color, thickness, cv::LINE_AA);

		cv::Point startPoint2(Bt.x, Bt.y);
		cv::Point endPoint2(Ct.x, Ct.y);

		line(showMat, startPoint2, endPoint2, color, thickness, cv::LINE_AA);

		cv::Point startPoint3(Ct.x, Ct.y);
		cv::Point endPoint3(At.x, At.y);

		line(showMat, startPoint3, endPoint3, color, thickness, cv::LINE_AA);
	}


	for (int i = 0; i != ts_1p2.size(); ++i)
	{
		Point2d A = ts_1p2[i][0];
		Point2d B = ts_1p2[i][1];
		Point2d C = ts_1p2[i][2];

		At.x = (A.x /*- src_mat.cols / 2*/) * 45;
		At.y = (A.y - src_mat.rows / 2) * 45;

		Bt.x = (B.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (B.y - src_mat.rows / 2) * 45;

		Ct.x = (C.x /*- src_mat.cols / 2*/) * 45;
		Ct.y = (C.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x>src_mat.cols * 15
			|| At.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x>src_mat.cols * 15
			|| Bt.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Ct.x < 0
			|| Ct.y < 0
			|| Ct.x>src_mat.cols * 15
			|| Ct.y>src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		cv::Scalar color(0, 255, 255);
		int thickness = 1;

		line(showMat, startPoint1, endPoint1, color, thickness, cv::LINE_AA);

		cv::Point startPoint2(Bt.x, Bt.y);
		cv::Point endPoint2(Ct.x, Ct.y);

		line(showMat, startPoint2, endPoint2, color, thickness, cv::LINE_AA);

		cv::Point startPoint3(Ct.x, Ct.y);
		cv::Point endPoint3(At.x, At.y);

		line(showMat, startPoint3, endPoint3, color, thickness, cv::LINE_AA);
	}

	for (int i = 0; i != ts_2p3.size(); ++i)
	{
		Point2d A = ts_2p3[i][0];
		Point2d B = ts_2p3[i][1];
		Point2d C = ts_2p3[i][2];

		At.x = (A.x /*- src_mat.cols / 2*/) * 45;
		At.y = (A.y - src_mat.rows / 2) * 45;

		Bt.x = (B.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (B.y - src_mat.rows / 2) * 45;

		Ct.x = (C.x /*- src_mat.cols / 2*/) * 45;
		Ct.y = (C.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x>src_mat.cols * 15
			|| At.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x>src_mat.cols * 15
			|| Bt.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Ct.x < 0
			|| Ct.y < 0
			|| Ct.x>src_mat.cols * 15
			|| Ct.y>src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		cv::Scalar color(130, 0, 75);
		int thickness = 1;

		line(showMat, startPoint1, endPoint1, color, thickness, cv::LINE_AA);

		cv::Point startPoint2(Bt.x, Bt.y);
		cv::Point endPoint2(Ct.x, Ct.y);

		line(showMat, startPoint2, endPoint2, color, thickness, cv::LINE_AA);

		cv::Point startPoint3(Ct.x, Ct.y);
		cv::Point endPoint3(At.x, At.y);

		line(showMat, startPoint3, endPoint3, color, thickness, cv::LINE_AA);
	}

	for (int i = 0; i != ts_1p3.size(); ++i)
	{
		Point2d A = ts_1p3[i][0];
		Point2d B = ts_1p3[i][1];
		Point2d C = ts_1p3[i][2];

		At.x = (A.x /*- src_mat.cols / 2*/) * 45;
		At.y = (A.y - src_mat.rows / 2) * 45;

		Bt.x = (B.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (B.y - src_mat.rows / 2) * 45;

		Ct.x = (C.x /*- src_mat.cols / 2*/) * 45;
		Ct.y = (C.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x>src_mat.cols * 15
			|| At.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x>src_mat.cols * 15
			|| Bt.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Ct.x < 0
			|| Ct.y < 0
			|| Ct.x>src_mat.cols * 15
			|| Ct.y>src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		cv::Scalar color(226, 43, 138);
		int thickness = 1;

		line(showMat, startPoint1, endPoint1, color, thickness, cv::LINE_AA);

		cv::Point startPoint2(Bt.x, Bt.y);
		cv::Point endPoint2(Ct.x, Ct.y);

		line(showMat, startPoint2, endPoint2, color, thickness, cv::LINE_AA);

		cv::Point startPoint3(Ct.x, Ct.y);
		cv::Point endPoint3(At.x, At.y);

		line(showMat, startPoint3, endPoint3, color, thickness, cv::LINE_AA);
	}

	int bdyNum = orderedBdy.size();
	for (int i = 0; i != bdyNum; ++i)
	{
		Point2d pre = orderedBdy[(i - 1 + bdyNum) % bdyNum];
		Point2d p = orderedBdy[i];

		At.x = (pre.x /*- src_mat.cols / 2*/) * 45;
		At.y = (pre.y - src_mat.rows / 2) * 45;

		Bt.x = (p.x /*- src_mat.cols / 2*/) * 45;
		Bt.y = (p.y - src_mat.rows / 2) * 45;

		if (At.x < 0
			|| At.y < 0
			|| At.x>src_mat.cols * 15
			|| At.y>src_mat.rows * 15)
		{
			continue;
		}

		if (Bt.x < 0
			|| Bt.y < 0
			|| Bt.x>src_mat.cols * 15
			|| Bt.y>src_mat.rows * 15)
		{
			continue;
		}

		cv::Point startPoint1(At.x, At.y);
		cv::Point endPoint1(Bt.x, Bt.y);

		showMat.at<Vec3b>(Bt.y, Bt.x) = Vec3b(0, 0, 255);
		line(showMat, startPoint1, endPoint1, color2, thickness2, cv::LINE_AA);
	}

	imwrite("F:\\showmat.bmp", showMat);
	TriClassify();
}

void ImgProc::OutputBdy(double2 *ps, int *idx, int &pN)
{
	int N = orderedBdy.size();
	pN = N;
	for (int i = 0; i != N; ++i)
	{
		ps[i][0] = orderedBdy[i].x;
		ps[i][1] = orderedBdy[i].y;
		idx[i] = i;
	}
}

void ImgProc::OutputPsBdy(double2 *innPs, int &pN, int *bdy, int &bN)
{
	bool *p_tag_s = new bool[ps.size()];
	memset(p_tag_s, 0, sizeof(bool)*ps.size());
	pN = 0;
	int idx;
	for (int i = 0; i != avail_tria.size(); ++i)
	{
		idx = avail_tria[i][0];
		if (false == p_tag_s[idx])
		{
			innPs[pN][0] = ps[idx].x;
			innPs[pN][1] = ps[idx].y;
			pN++;
			p_tag_s[idx] = true;
		}

		idx = avail_tria[i][1];
		if (false == p_tag_s[idx])
		{
			innPs[pN][0] = ps[idx].x;
			innPs[pN][1] = ps[idx].y;
			pN++;
			p_tag_s[idx] = true;
		}

		idx = avail_tria[i][2];
		if (false == p_tag_s[idx])
		{
			innPs[pN][0] = ps[idx].x;
			innPs[pN][1] = ps[idx].y;
			pN++;
			p_tag_s[idx] = true;
		}
	}
	// copy model points to M_data
	M_data4P.resize(boost::extents[pN][2]);
	for (int i = 0; i != pN; ++i)
	{
		M_data4P[i][0] = innPs[i][0];
		M_data4P[i][1] = innPs[i][1];
	}

	// build a kd tree from the model point cloud
	M_tree4P = new kdtree::KDTree(M_data4P);
	kdtree::KDTreeResultVector result;
	// kd tree query + result
	std::vector<float> query(2);

	bN = orderedBdy.size();
	for (int i = 0; i != bN; ++i)
	{
		query[0] = orderedBdy[i].x;
		query[1] = orderedBdy[i].y;

		// search nearest neighbor
		M_tree4P->n_nearest(query, 1, result);

		bdy[i] = result[0].idx + 1;
	}

	delete[]p_tag_s;
}

double ImgProc::EuclidDist(double2 &v1, double2 &v2)
{
	return sqrt(pow(v1[0] - v2[0], 2)
		+ pow(v1[1] - v2[1], 2));
}

void ImgProc::AvailModel4Output(double2 *ps_op, int &pN, uchar *gs_op, int &gN, int3 *tri_set, int &tN)
{
	assert(ps.size() == gs.size());
	pN = ps.size();
	gN = gs.size();
	for (int i = 0; i != pN; ++i)
	{
		ps_op[i][0] = this->ps[i].x;
		ps_op[i][1] = this->ps[i].y;
		gs_op[i] = this->gs[i];
	}

	int j = 0;
	for (int i = 0; i != avail_tria.size(); ++i)
	{
		Point2d ctr = findCentroid(this->ps[this->avail_tria[i][0]]
			, this->ps[this->avail_tria[i][1]]
			, this->ps[this->avail_tria[i][2]]);
		if (ctr.x > c / 2 - 1
			&& ctr.x < c / 2 + 1
			&& ctr.y > r / 2 - 1
			&& ctr.y < r / 2 + 1)
		{
			continue;
		}
		tri_set[j][0] = this->avail_tria[i][0];
		tri_set[j][1] = this->avail_tria[i][2];
		tri_set[j][2] = this->avail_tria[i][1];
		j++;
	}
	tN = j;
	printf("avail_tria's number is %d\n", j);


}

void ImgProc::OuterModel4Output(double2 *ps_op, int &pN, uchar *gs_op, int &gN, int3 *tri_set, int &tN)
{
	assert(ps.size() == gs.size());
	pN = ps.size();
	gN = gs.size();
	for (int i = 0; i != pN; ++i)
	{
		ps_op[i][0] = this->ps[i].x;
		ps_op[i][1] = this->ps[i].y;
		gs_op[i] = this->gs[i];
	}

	tN = outer_tria.size();
	printf("outer_tria's number is %d\n", tN);
	for (int i = 0; i != tN; ++i)
	{
		tri_set[i][0] = this->outer_tria[i][0];
		tri_set[i][1] = this->outer_tria[i][2];
		tri_set[i][2] = this->outer_tria[i][1];
	}
}

void ImgProc::Sobel4Edge(Mat &sMat, Mat &srcMat)
{
	//cv::Mat blurred_image;
	//cv::GaussianBlur(srcMat, blurred_image, cv::Size(5, 5), 1.5);
	//cv::GaussianBlur(blurred_image, blurred_image, cv::Size(5, 5), 1.5);
	cv::Mat sobel_x, sobel_y;
	cv::Sobel(srcMat, sobel_x, CV_64F, 1, 0, 3);
	cv::Sobel(srcMat, sobel_y, CV_64F, 0, 1, 3);

	cv::Mat sobel_magnitude;
	cv::magnitude(sobel_x, sobel_y, sobel_magnitude);
	cv::normalize(sobel_magnitude, sMat, 0, 255, cv::NORM_MINMAX);
	sMat.convertTo(sMat, CV_8U);
}


uchar ImgProc::GetinterpoGray(Mat &mat, Point2d offset, int i, int j)
{
	uchar g00, g01, g10, g11;
	double a, b;
	if (offset.x < 0 && offset.y < 0)
	{
		g00 = mat.at<uchar>(i - 1, j - 1);
		g01 = mat.at<uchar>(i - 1, j);
		g10 = mat.at<uchar>(i, j - 1);
		g11 = mat.at<uchar>(i, j);
		a = 1 + offset.y;
		b = 1 + offset.x;
		return BilineInterpola(g00, g01, g10, g11, a, b);
	}
	else if (offset.x > 0 && offset.y < 0)
	{
		g00 = mat.at<uchar>(i - 1, j);
		g01 = mat.at<uchar>(i - 1, j + 1);
		g10 = mat.at<uchar>(i, j);
		g11 = mat.at<uchar>(i, j + 1);
		a = 1 + offset.y;
		b = offset.x;
		return BilineInterpola(g00, g01, g10, g11, a, b);
	}
	else if (offset.x < 0 && offset.y > 0)
	{
		g00 = mat.at<uchar>(i, j - 1);
		g01 = mat.at<uchar>(i, j);
		g10 = mat.at<uchar>(i + 1, j - 1);
		g11 = mat.at<uchar>(i + 1, j);
		a = offset.y;
		b = 1 + offset.x;
		return BilineInterpola(g00, g01, g10, g11, a, b);
	}
	else if (offset.x > 0 && offset.y > 0)
	{
		g00 = mat.at<uchar>(i, j);
		g01 = mat.at<uchar>(i, j + 1);
		g10 = mat.at<uchar>(i + 1, j);
		g11 = mat.at<uchar>(i + 1, j + 1);
		a = offset.y;
		b = offset.x;
		return BilineInterpola(g00, g01, g10, g11, a, b);
	}
	else
		return 0;
}


void ImgProc::Convert2IdxS(vector<vector<int>> &idxSet
	, vector<Point2d> &addPs, vector<vector<Point2d>> &addTs)
{
	int M_num = addPs.size();
	uchar gc;
	for (int m = 0; m < M_num; m++)
	{
		ps.push_back(addPs[m]);
		gc = getGray(src_mat, addPs[m]);
		gs.push_back(gc);
	}
	M_num = ps.size();
	// copy model points to M_data
	M_data4Padd.resize(boost::extents[M_num][2]);
	for (int m = 0; m < M_num; m++)
	{
		M_data4Padd[m][0] = ps[m].x;
		M_data4Padd[m][1] = ps[m].y;
	}

	// build a kd tree from the model point cloud
	M_tree4Padd = new kdtree::KDTree(M_data4Padd);
	kdtree::KDTreeResultVector result;
	// kd tree query + result
	std::vector<float>         query(2);
	int t_n = addTs.size();

	for (int i = 0; i != t_n; ++i)
	{
		vector<int> idx;
		vector<Point2d> t = addTs[i];

		query[0] = t[0].x;
		query[1] = t[0].y;
		// search nearest neighbor
		M_tree4Padd->n_nearest(query, 1, result);
		idx.push_back(result[0].idx);

		query[0] = t[1].x;
		query[1] = t[1].y;
		// search nearest neighbor
		M_tree4Padd->n_nearest(query, 1, result);
		idx.push_back(result[0].idx);

		query[0] = t[2].x;
		query[1] = t[2].y;
		// search nearest neighbor
		M_tree4Padd->n_nearest(query, 1, result);
		idx.push_back(result[0].idx);

		idxSet.push_back(idx);
	}

	delete M_tree4Padd;
}

void ImgProc::replacePs(vector<Point2d> &ps, double thres)
{
	float dist;
	int M_num = ps.size();
	// copy model points to M_data
	M_data4P.resize(boost::extents[M_num][2]);
	for (int m = 0; m < M_num; m++)
	{
		M_data4P[m][0] = ps[m].x;
		M_data4P[m][1] = ps[m].y;
	}

	// build a kd tree from the model point cloud
	M_tree4P = new kdtree::KDTree(M_data4P);
	kdtree::KDTreeResultVector result;
	// kd tree query + result
	std::vector<float> query(2);
	int b_n = orderedBdy.size();
	for (int i = 0; i != b_n; ++i)
	{
		int pre_i = (b_n + i - 1) % b_n;
		query[0] = (orderedBdy[i].x + orderedBdy[pre_i].x) / 2;
		query[1] = (orderedBdy[i].y + orderedBdy[pre_i].y) / 2;

		// search nearest neighbor
		M_tree4P->n_nearest(query, 16, result);

		for (int j = 0; j != 16; ++j)
		{
			dist = pointToSegmentDistance(ps[result[j].idx]
				, orderedBdy[i], orderedBdy[pre_i]);

			if (dist > EPSILON && dist < thres)
			{
				FootResult fr = repelFromSegmentIfTooClose(ps[result[j].idx]
					, orderedBdy[i], orderedBdy[pre_i], thres);

				if (fr.adjusted)
				{
					ps[result[j].idx] = fr.adjustedPoint;
					gs[result[j].idx] = getGray(src_mat, fr.adjustedPoint);
				}

			}
		}

	}

	delete M_tree4P;
}

void ImgProc::replaceBdy(vector<Point2d> &ps, vector<int> &untreatedIds, double thres)
{
	float dist;
	int M_num = ps.size();
	// copy model points to M_data
	M_data4P.resize(boost::extents[M_num][2]);
	for (int m = 0; m < M_num; m++)
	{
		M_data4P[m][0] = ps[m].x;
		M_data4P[m][1] = ps[m].y;
	}

	// build a kd tree from the model point cloud
	M_tree4P = new kdtree::KDTree(M_data4P);
	kdtree::KDTreeResultVector result;
	// kd tree query + result
	std::vector<float>         query(2);
	int b_n = orderedBdy.size();
	for (int i = 0; i != b_n; ++i)
	{
		query[0] = orderedBdy[i].x;
		query[1] = orderedBdy[i].y;

		// search nearest neighbor
		M_tree4P->n_nearest(query, 1, result);
		dist = norm(ps[result[0].idx] - orderedBdy[i]);

		if (dist < thres)
		{
			ps[result[0].idx] = orderedBdy[i];
			uchar data = src_mat.at<uchar>(
				floor(orderedBdy[i].y), floor(orderedBdy[i].x));
			gs[result[0].idx] = data;
		}
		else
		{
			untreatedIds.push_back(i);
		}
	}
	delete M_tree4P;
}

void ImgProc::DrawBdy(const char *filePath)
{
	int bN = orderedBdy.size();
	for (int i = 0; i != bN; ++i)
	{
		int idx_pre = (bN + i - 1) % bN;
		Point2d p = orderedBdy[i];
		Point2d p_pre = orderedBdy[idx_pre];
		line(ori_mat, Point(p.x - 0.5, p.y - 0.5)
			, Point(p_pre.x - 0.5, p_pre.y - 0.5), Scalar(0, 0, 255), 1, LINE_AA);
	}

	imwrite(filePath, ori_mat);
}

void ImgProc::InitKdTree4Tri(vector<Point2d> &triCtrs)
{
	int M_num = triCtrs.size();
	// copy model points to M_data
	M_data4Tri.resize(boost::extents[M_num][2]);
	for (int m = 0; m < M_num; m++)
	{
		M_data4Tri[m][0] = triCtrs[m].x;
		M_data4Tri[m][1] = triCtrs[m].y;
	}

	// build a kd tree from the model point cloud
	M_tree4Tri = new kdtree::KDTree(M_data4Tri);
}

void ImgProc::HandleNearEdgePt(vector<vector<Point2d>> &addTriS, vector<Point2d> &addPs
	, vector<int> &delTriS, int idx, double thres)
{
	kdtree::KDTreeResultVector result;
	// kd tree query + result
	std::vector<float> query(2);
	Point2d p, A, B, C;
	int i;
	int wBnum = orderedBdy.size();
	double distAB, distBC, distCA;

	p = orderedBdy[idx];
	query[0] = p.x;
	query[1] = p.y;

	// search nearest neighbor
	M_tree4Tri->n_nearest(query, SEARCHITEMNUM, result);
	vector<int> triId;
	double minDist = MAXV;
	Point2d ectr;

	for (int j = 0; j != SEARCHITEMNUM; ++j)
	{
		triId = triangle_set[result[j].idx];
		A = ps[triId[0]];
		B = ps[triId[1]];
		C = ps[triId[2]];

		if (norm(A - p) < EPSILON || norm(B - p) < EPSILON || norm(C - p) < EPSILON)
		{
			break;
		}

		if (1 == pointInTriangle(A, B, C, p))
		{
			distAB = pointToSegmentDistance(p, A, B);
			distBC = pointToSegmentDistance(p, B, C);
			distCA = pointToSegmentDistance(p, C, A);

			minDist = min(min(distAB, distBC), distCA);

			if (abs(distAB - minDist) < EPSILON)
			{
				ectr.x = (A.x + B.x) / 2;
				ectr.y = (A.y + B.y) / 2;
			}
			else if (abs(distBC - minDist) < EPSILON)
			{
				ectr.x = (B.x + C.x) / 2;
				ectr.y = (B.y + C.y) / 2;
			}
			else if (abs(distCA - minDist) < EPSILON)
			{
				ectr.x = (C.x + A.x) / 2;
				ectr.y = (C.y + A.y) / 2;
			}

			break;
		}
	}

	if (minDist < thres)
	{
		int interpolNum = 0;
		for (int j = 0; j != SEARCHITEMNUM; ++j)
		{
			triId = triangle_set[result[j].idx];
			A = ps[triId[0]];
			B = ps[triId[1]];
			C = ps[triId[2]];

			Point2d ctrAB, ctrBC, ctrCA;
			ctrAB.x = (A.x + B.x) / 2;
			ctrAB.y = (A.y + B.y) / 2;

			ctrBC.x = (B.x + C.x) / 2;
			ctrBC.y = (B.y + C.y) / 2;

			ctrCA.x = (C.x + A.x) / 2;
			ctrCA.y = (C.y + A.y) / 2;

			if (norm(ctrAB - ectr) < EPSILON)
			{
				delTriS.push_back(result[j].idx);
				splitEdge1(addTriS, A, B, C, p);
				interpolNum++;
			}
			else if (norm(ctrBC - ectr) < EPSILON)
			{
				delTriS.push_back(result[j].idx);
				splitEdge1(addTriS, B, C, A, p);
				interpolNum++;
			}
			else if (norm(ctrCA - ectr) < EPSILON)
			{
				delTriS.push_back(result[j].idx);
				splitEdge1(addTriS, C, A, B, p);
				interpolNum++;
			}
			else
			{

			}

		}

		if (1 <= interpolNum)
		{
			printf("No %d Error:degree of interpolation is %d\n\n", idx, interpolNum);
		}
	}
}

void ImgProc::ClipTri(vector<int> &intersectTridxSet
	, vector<vector<LineSegIntersec>> &lSegIntersecSet
	, vector<Point2d> &addPs, vector<Point2d> &edgeCtrs
	, vector<Point2d> &triCtrs, vector<int> &untreatedIds)
{

	int M_num = triCtrs.size();
	// copy model points to M_data
	M_data4Tri.resize(boost::extents[M_num][2]);
	for (int m = 0; m < M_num; m++)
	{
		M_data4Tri[m][0] = triCtrs[m].x;
		M_data4Tri[m][1] = triCtrs[m].y;
	}

	// build a kd tree from the model point cloud
	M_tree4Tri = new kdtree::KDTree(M_data4Tri);
	kdtree::KDTreeResultVector result;
	// kd tree query + result
	std::vector<float>         query(2);

	int wBnum = orderedBdy.size();
	Point2d p1, p2, A, B, C, interior;
	int i_p, i, res1, res2, idx;

	for (i = 0; i != wBnum; ++i)
	{
		i_p = (wBnum + i - 1) % wBnum;
		p1 = orderedBdy[i_p];
		p2 = orderedBdy[i];

		query[0] = edgeCtrs[i].x;
		query[1] = edgeCtrs[i].y;

		// search nearest neighbor
		M_tree4Tri->n_nearest(query, SEARCHITEMNUM, result);
		vector<int> triId;
		for (int j = 0; j != SEARCHITEMNUM; ++j)
		{
			triId = triangle_set[result[j].idx];
			A = ps[triId[0]];
			B = ps[triId[1]];
			C = ps[triId[2]];
			bool p1AB = containsPoint(A, B, p1);
			bool p2AB = containsPoint(A, B, p2);
			bool p1BC = containsPoint(B, C, p1);
			bool p2BC = containsPoint(B, C, p2);
			bool p1CA = containsPoint(C, A, p1);
			bool p2CA = containsPoint(C, A, p2);
			res1 = pointInTriangle(A, B, C, p1);
			res2 = pointInTriangle(A, B, C, p2);
			if (1 == res1)
			{
				interior = p1;
			}
			else if (1 == res2)
			{
				interior = p2;
			}
			int mask;
			vector<Point2d> intersections = triaLineSegIntersec(A, B, C, p1, p2, mask);
			int intersectNum = intersections.size();

			if (1 == intersectNum)
			{
				if (1 == res1 || 1 == res2)
				{
					idx = findPointIndex(addPs, intersections[0]);
					if (-1 == idx)
					{
						addPs.push_back(intersections[0]);
					}

					if (mask & 2 || mask & 4 || mask & 8)
					{
						idx = isExi(intersectTridxSet, result[j].idx);
						if (-1 == idx)
						{
							intersectTridxSet.push_back(result[j].idx);
							vector<LineSegIntersec> lsis;
							LineSegIntersec lsi;
							lsi.ps = intersections;
							lsi.p1 = p1;
							lsi.p2 = p2;
							lsi.res1 = res1;
							lsi.res2 = res2;
							lsi.mask = mask;
							lsis.push_back(lsi);
							lSegIntersecSet.push_back(lsis);
						}
						else
						{
							LineSegIntersec lsi;
							lsi.ps = intersections;
							lsi.p1 = p1;
							lsi.p2 = p2;
							lsi.res1 = res1;
							lsi.res2 = res2;
							lsi.mask = mask;
							lSegIntersecSet[idx].push_back(lsi);
						}
					}
					else
					{
						printf("-----------res1=%d res2=%d intersectnum=%d------------\n"
							, res1, res2, intersectNum);
					}
				}
				else if ((1 != res1 && 1 != res2) && (0 == res1 || 0 == res2))
				{
					if (2 == res1 || 2 == res2
						|| 3 == res1 || 3 == res2
						|| 4 == res1 || 4 == res2)
					{
						idx = isExi(intersectTridxSet, result[j].idx);
						if (-1 == idx)
						{
							intersectTridxSet.push_back(result[j].idx);
							vector<LineSegIntersec> lsis;
							LineSegIntersec lsi;
							lsi.ps = intersections;
							lsi.p1 = p1;
							lsi.p2 = p2;
							lsi.res1 = res1;
							lsi.res2 = res2;
							lsi.mask = mask;
							lsis.push_back(lsi);
							lSegIntersecSet.push_back(lsis);
						}
						else
						{
							LineSegIntersec lsi;
							lsi.ps = intersections;
							lsi.p1 = p1;
							lsi.p2 = p2;
							lsi.res1 = res1;
							lsi.res2 = res2;
							lsi.mask = mask;
							lSegIntersecSet[idx].push_back(lsi);
						}
					}
					else
					{
						printf("-----------res1=%d res2=%d intersectnum=%d------------\n"
							, res1, res2, intersectNum);
					}
				}
				else
				{
					printf("-----------res1=%d res2=%d intersectnum=%d------------\n"
						, res1, res2, intersectNum);
				}
			}
			else if (2 == intersectNum)
			{
				double dist = norm(intersections[1] - intersections[0]);
				if (dist < EPSILON)
				{
					if (1 == res1 || 1 == res2)
					{
						if (norm(A - intersections[0]) < EPSILON
							|| norm(B - intersections[0]) < EPSILON
							|| norm(C - intersections[0]) < EPSILON)
						{
							idx = isExi(intersectTridxSet, result[j].idx);
							if (-1 == idx)
							{
								intersectTridxSet.push_back(result[j].idx);
								vector<LineSegIntersec> lsis;
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lsis.push_back(lsi);
								lSegIntersecSet.push_back(lsis);
							}
							else
							{
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lSegIntersecSet[idx].push_back(lsi);
							}
						}
						else
						{
							printf("res1=%d res2=%d intersectnum=%d\n", res1, res2, intersectNum);
						}
					}
					else
					{
						printf("res1=%d res2=%d intersectnum=%d\n", res1, res2, intersectNum);
					}
				}
				else //the two intersection points are far enough apart
				{

					if (0 == res1 && 0 == res2)
					{
						if (norm(A - intersections[0]) < EPSILON
							|| norm(B - intersections[0]) < EPSILON
							|| norm(C - intersections[0]) < EPSILON)
						{
							printf("res1=%d res2=%d intersectnum=%d\n", res1, res2, intersectNum);
						}
						else
						{
							idx = findPointIndex(addPs, intersections[0]);
							if (-1 == idx)
							{
								addPs.push_back(intersections[0]);
							}

							idx = findPointIndex(addPs, intersections[1]);
							if (-1 == idx)
							{
								addPs.push_back(intersections[1]);
							}

							if ((mask & 2 && mask & 4)
								|| (mask & 4 && mask & 8)
								|| (mask & 2 && mask & 8))
							{
								idx = isExi(intersectTridxSet, result[j].idx);
								if (-1 == idx)
								{
									intersectTridxSet.push_back(result[j].idx);
									vector<LineSegIntersec> lsis;
									LineSegIntersec lsi;
									lsi.ps = intersections;
									lsi.p1 = p1;
									lsi.p2 = p2;
									lsi.res1 = res1;
									lsi.res2 = res2;
									lsi.mask = mask;
									lsis.push_back(lsi);
									lSegIntersecSet.push_back(lsis);
								}
								else
								{
									LineSegIntersec lsi;
									lsi.ps = intersections;
									lsi.p1 = p1;
									lsi.p2 = p2;
									lsi.res1 = res1;
									lsi.res2 = res2;
									lsi.mask = mask;
									lSegIntersecSet[idx].push_back(lsi);
								}
							}
							else
							{
								printf("res1=%d res2=%d intersectnum=%d\n", res1, res2, intersectNum);
							}
						}
					}
					else if (areCollinear(A, B, p1, p2)
						|| areCollinear(B, C, p1, p2)
						|| areCollinear(C, A, p1, p2))
					{
						if (p1AB || p2AB || p1BC || p2BC || p1CA || p2CA)
						{
							idx = isExi(intersectTridxSet, result[j].idx);
							if (-1 == idx)
							{
								intersectTridxSet.push_back(result[j].idx);
								vector<LineSegIntersec> lsis;
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lsis.push_back(lsi);
								lSegIntersecSet.push_back(lsis);
							}
							else
							{
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lSegIntersecSet[idx].push_back(lsi);
							}
						}
						else
						{
							printf("res1=%d res2=%d intersectnum=%d\n"
								, res1, res2, intersectNum);
						}
					}
					else if (p1AB || p1BC || p1CA || p2AB || p2BC || p2CA)
					{
						if (norm(p1 - intersections[0]) < EPSILON
							|| norm(p2 - intersections[0]) < EPSILON)
						{
							idx = findPointIndex(addPs, intersections[1]);
							if (-1 == idx)
							{
								addPs.push_back(intersections[1]);
							}
						}
						else if (norm(p1 - intersections[1]) < EPSILON
							|| norm(p2 - intersections[1]) < EPSILON)
						{
							idx = findPointIndex(addPs, intersections[0]);
							if (-1 == idx)
							{
								addPs.push_back(intersections[0]);
							}
						}

						if ((mask & 2 && mask & 4)
							|| (mask & 4 && mask & 8)
							|| (mask & 2 && mask & 8))
						{
							idx = isExi(intersectTridxSet, result[j].idx);
							if (-1 == idx)
							{
								intersectTridxSet.push_back(result[j].idx);
								vector<LineSegIntersec> lsis;
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lsis.push_back(lsi);
								lSegIntersecSet.push_back(lsis);
							}
							else
							{
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lSegIntersecSet[idx].push_back(lsi);
							}
						}
						else
						{
							printf("res1=%d res2=%d intersectnum=%d\n"
								, res1, res2, intersectNum);
						}
					}
					else
					{
						printf("res1=%d res2=%d intersectnum=%d\n"
							, res1, res2, intersectNum);
					}
				}
			}
			else if (3 == intersectNum)
			{
				vector<Point2d> colliPs;
				vector<Point2d> unColliPs;
				if (containsPoint(p1, p2, A))
				{
					colliPs.push_back(A);
				}
				else
				{
					unColliPs.push_back(A);
				}

				if (containsPoint(p1, p2, B))
				{
					colliPs.push_back(B);
				}
				else
				{
					unColliPs.push_back(B);
				}

				if (containsPoint(p1, p2, C))
				{
					colliPs.push_back(C);
				}
				else
				{
					unColliPs.push_back(C);
				}

				int colliNum = colliPs.size();
				if (1 == colliNum)
				{
					if (norm(p1 - colliPs[0]) < EPSILON
						|| norm(p2 - colliPs[0]) < EPSILON)
					{
						if (mask & 2 && mask & 4 && mask & 8)
						{
							if (norm(intersections[1] - intersections[2]) < EPSILON)
							{
								idx = findPointIndex(addPs, intersections[0]);
								if (-1 == idx)
								{
									addPs.push_back(intersections[0]);
								}
							}
							else if (norm(intersections[0] - intersections[2]) < EPSILON)
							{
								idx = findPointIndex(addPs, intersections[1]);
								if (-1 == idx)
								{
									addPs.push_back(intersections[1]);
								}
							}
							else if (norm(intersections[0] - intersections[1]) < EPSILON)
							{
								idx = findPointIndex(addPs, intersections[2]);
								if (-1 == idx)
								{
									addPs.push_back(intersections[2]);
								}
							}

							idx = isExi(intersectTridxSet, result[j].idx);
							if (-1 == idx)
							{
								intersectTridxSet.push_back(result[j].idx);
								vector<LineSegIntersec> lsis;
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lsis.push_back(lsi);
								lSegIntersecSet.push_back(lsis);
							}
							else
							{
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lSegIntersecSet[idx].push_back(lsi);
							}
						}
						else
						{
							printf("res1=%d res2=%d intersectnum=%d\n"
								, res1, res2, intersectNum);
						}
					}
					else if (areCollinear(p1, colliPs[0], colliPs[0], p2))
					{
						if (p1AB
							|| p1BC
							|| p1CA
							|| p2AB
							|| p2BC
							|| p2CA)
						{
							idx = isExi(intersectTridxSet, result[j].idx);
							if (-1 == idx)
							{
								intersectTridxSet.push_back(result[j].idx);
								vector<LineSegIntersec> lsis;
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lsis.push_back(lsi);
								lSegIntersecSet.push_back(lsis);
							}
							else
							{
								LineSegIntersec lsi;
								lsi.ps = intersections;
								lsi.p1 = p1;
								lsi.p2 = p2;
								lsi.res1 = res1;
								lsi.res2 = res2;
								lsi.mask = mask;
								lSegIntersecSet[idx].push_back(lsi);
							}
						}
						else
						{
							printf("res1=%d res2=%d intersectnum=%d\n"
								, res1, res2, intersectNum);
						}
					}
					else
					{
						printf("res1=%d res2=%d intersectnum=%d\n"
							, res1, res2, intersectNum);
					}
				}
				else
				{
					printf("res1=%d res2=%d intersectnum=%d\n"
						, res1, res2, intersectNum);
				}
			}
			else
			{
				printf("res1=%d res2=%d intersectnum=%d\n"
					, res1, res2, intersectNum);
			}
		}
	}

	delete M_tree4Tri;
}

void ImgProc::TriClassify()
{
	int tN = triangle_set.size();
	vector<int> idx;
	Point2d p1, p2, p3;

	for (int i = 0; i != tN; ++i)
	{
		idx = triangle_set[i];
		p1 = ps[idx[0]];
		p2 = ps[idx[1]];
		p3 = ps[idx[2]];
		Point2d ctr = findCentroid(p1, p2, p3);

		if (isPointInsidePolygon(ctr, orderedBdy))
		{
			avail_tria.push_back(idx);
		}
		else
		{
			outer_tria.push_back(idx);
		}
	}
}

void ImgProc::Triangula4Edge(vector<LineSegIntersec> &lSegIntersec, int &tridx)
{
	Point2d A = ps[triangle_set[tridx][0]];
	Point2d B = ps[triangle_set[tridx][1]];
	Point2d C = ps[triangle_set[tridx][2]];

	int lSegNum = lSegIntersec.size();
	if (2 == lSegNum)
	{
		Triangula4Edge2(lSegIntersec, A, B, C);
	}
	else if (1 == lSegNum)
	{
		Triangula4Edge1(lSegIntersec, A, B, C);
	}
	else
	{
		printf("Error: The number of adjacent line ");
		printf("segments intersecting with the triangle is %d.\n", lSegNum);
	}
}

void ImgProc::Triangula4Edge1(vector<LineSegIntersec> &lSegIntersec
	, Point2d A, Point2d B, Point2d C)
{
	int num = lSegIntersec[0].ps.size();

	if (2 == num)
	{
		if (norm(lSegIntersec[0].ps[0] - lSegIntersec[0].ps[1]) < EPSILON)
		{
			printf("The case where a line segment intersects at a vertex.\n");
		}
		else
		{
			if ((lSegIntersec[0].mask & 2) && (lSegIntersec[0].mask & 4))
			{
				split2(ts_1, A, B, C, lSegIntersec[0].ps[0], lSegIntersec[0].ps[1]);
			}
			else if ((lSegIntersec[0].mask & 4) && (lSegIntersec[0].mask & 8))
			{
				split2(ts_1, B, C, A, lSegIntersec[0].ps[0], lSegIntersec[0].ps[1]);
			}
			else if ((lSegIntersec[0].mask & 2) && (lSegIntersec[0].mask & 8))
			{
				split2(ts_1, C, A, B, lSegIntersec[0].ps[1], lSegIntersec[0].ps[0]);
			}
			else
			{
				printf("Error:single segment's two endpoints.\n");
			}
		}

	}
	else if (3 == num)
	{

		if (areCollinear(lSegIntersec[0].p1, lSegIntersec[0].p2, A, B)
			|| areCollinear(lSegIntersec[0].p1, lSegIntersec[0].p2, B, C)
			|| areCollinear(lSegIntersec[0].p1, lSegIntersec[0].p2, C, A))
		{
			printf("The edge coincides with the triangle edge.\n");
		}
		else
		{
			Point2d bc;
			int bcId;
			if (norm(lSegIntersec[0].ps[0] - lSegIntersec[0].ps[1]) < EPSILON)
			{
				bc = lSegIntersec[0].ps[2];
				bcId = 2;
			}
			else if (norm(lSegIntersec[0].ps[0] - lSegIntersec[0].ps[2]) < EPSILON)
			{
				bc = lSegIntersec[0].ps[1];
				bcId = 1;
			}
			else if (norm(lSegIntersec[0].ps[1] - lSegIntersec[0].ps[2]) < EPSILON)
			{
				bc = lSegIntersec[0].ps[0];
				bcId = 0;
			}
			else
			{
				printf("Error:single segment's intersection with edge.\n");
			}

			if (0 == bcId)
			{
				splitEdge1(ts_1, A, B, C, bc);
			}
			else if (1 == bcId)
			{
				splitEdge1(ts_1, B, C, A, bc);
			}
			else if (2 == bcId)
			{
				splitEdge1(ts_1, C, A, B, bc);
			}
			else
			{
				printf("Error:single segment's location of intersection.\n");
			}
		}
	}
	else
	{
		printf("Error:single segment's intersecting number-- %d.\n", num);
	}
}

void ImgProc::Triangula4Edge2(vector<LineSegIntersec> &lSegIntersec
	, Point2d A, Point2d B, Point2d C)
{
	int lnum = lSegIntersec[0].ps.size();
	int rnum = lSegIntersec[1].ps.size();

	if (1 == lnum && 1 == rnum)
	{
		if (Triangula4Edge2_1p1(ts_1p1, lSegIntersec, A, B, C))
		{
			//printf("1,1 intersecting is successful.\n");
		}
		else
		{
			printf("1,1 intersecting is failing.\n");
		}
	}
	else if (2 == max(lnum, rnum) && 1 == min(lnum, rnum))
	{
		if (Triangula4Edge2_1p2(ts_1p2, lSegIntersec, A, B, C))
		{
			//printf("1,2 intersecting is successful.\n");
		}
		else
		{
			printf("1,2 intersecting is failing.\n");
		}
	}
	else if (2 == lnum && 2 == rnum)
	{
		if (Triangula4Edge2_2p2(ts_2p2, lSegIntersec, A, B, C))
		{
			//printf("2,2 intersecting is successful.\n");
		}
		else
		{
			printf("2,2 intersecting is failing.\n");
		}
	}
	else if (3 == max(lnum, rnum) && 1 == min(lnum, rnum))
	{
		if (Triangula4Edge2_1p3(ts_1p3, lSegIntersec, A, B, C))
		{
			//printf("1,3 intersecting is successful.\n");
		}
		else
		{
			printf("1,3 intersecting is failing.\n");
		}
	}
	else if (3 == max(lnum, rnum) && 2 == min(lnum, rnum))
	{
		if (Triangula4Edge2_2p3(ts_2p3, lSegIntersec, A, B, C))
		{
			//printf("2,3 intersecting is successful.\n");
		}
		else
		{
			printf("2,3 intersecting is failing.\n");
		}
	}
	else
	{
		printf("Error: lnum=%d rnum=%d.\n", lnum, rnum);
	}


}

int ImgProc::Triangula4Edge2_1p1(vector<vector<Point2d>> &ts_1p1
	, vector<LineSegIntersec> &lSegIntersec
	, Point2d A, Point2d B, Point2d C)
{
	Point2d intersec;
	int res1, res2;

	if (!GetIntersecParams(intersec, res1, res2, lSegIntersec))
	{
		printf("Error: intersecting params");
		return 0;
	}

	if (1 == res1 && 1 == res2)
	{
		if (lSegIntersec[0].mask == lSegIntersec[1].mask)
		{
			Point2d pa, pb;
			if (2 == lSegIntersec[0].mask)
			{
				if (norm(A - lSegIntersec[0].ps[0])
			> norm(A - lSegIntersec[1].ps[0]))
				{
					pa = lSegIntersec[1].ps[0];
					pb = lSegIntersec[0].ps[0];
				}
				else
				{
					pa = lSegIntersec[0].ps[0];
					pb = lSegIntersec[1].ps[0];
				}

				split1p1interior1(ts_1p1, A, B, C, intersec, pa, pb);
				return 1;
			}
			else if (4 == lSegIntersec[0].mask)
			{
				if (norm(B - lSegIntersec[0].ps[0])
			> norm(B - lSegIntersec[1].ps[0]))
				{
					pa = lSegIntersec[1].ps[0];
					pb = lSegIntersec[0].ps[0];
				}
				else
				{
					pa = lSegIntersec[0].ps[0];
					pb = lSegIntersec[1].ps[0];
				}

				split1p1interior1(ts_1p1, B, C, A, intersec, pa, pb);
				return 1;
			}
			else if (8 == lSegIntersec[0].mask)
			{
				if (norm(C - lSegIntersec[0].ps[0])
			> norm(C - lSegIntersec[1].ps[0]))
				{
					pa = lSegIntersec[1].ps[0];
					pb = lSegIntersec[0].ps[0];
				}
				else
				{
					pa = lSegIntersec[0].ps[0];
					pb = lSegIntersec[1].ps[0];
				}

				split1p1interior1(ts_1p1, C, A, B, intersec, pa, pb);
				return 1;
			}
			else
			{
				printf("An error occurred for two equal masks.\n");
				return 0;
			}
		}
		else
		{
			Point2d ab, bc;
			if (2 == lSegIntersec[0].mask && 4 == lSegIntersec[1].mask)
			{
				ab = lSegIntersec[0].ps[0];
				bc = lSegIntersec[1].ps[0];
				split1p1interior2(ts_1p1, A, B, C, intersec, ab, bc);
				return 1;
			}
			else if (2 == lSegIntersec[1].mask && 4 == lSegIntersec[0].mask)
			{
				ab = lSegIntersec[1].ps[0];
				bc = lSegIntersec[0].ps[0];
				split1p1interior2(ts_1p1, A, B, C, intersec, ab, bc);
				return 1;
			}
			else if (4 == lSegIntersec[0].mask && 8 == lSegIntersec[1].mask)
			{
				ab = lSegIntersec[0].ps[0];
				bc = lSegIntersec[1].ps[0];
				split1p1interior2(ts_1p1, B, C, A, intersec, ab, bc);
				return 1;
			}
			else if (4 == lSegIntersec[1].mask && 8 == lSegIntersec[0].mask)
			{
				ab = lSegIntersec[1].ps[0];
				bc = lSegIntersec[0].ps[0];
				split1p1interior2(ts_1p1, B, C, A, intersec, ab, bc);
				return 1;
			}
			else if (8 == lSegIntersec[0].mask && 2 == lSegIntersec[1].mask)
			{
				ab = lSegIntersec[0].ps[0];
				bc = lSegIntersec[1].ps[0];
				split1p1interior2(ts_1p1, C, A, B, intersec, ab, bc);
				return 1;
			}
			else if (8 == lSegIntersec[1].mask && 2 == lSegIntersec[0].mask)
			{
				ab = lSegIntersec[1].ps[0];
				bc = lSegIntersec[0].ps[0];
				split1p1interior2(ts_1p1, C, A, B, intersec, ab, bc);
				return 1;
			}
			else
			{
				printf("Two different masks are incorrect.\n");
				return 0;
			}
		}
	}
	else if (1 < res1 && res1 == res2)
	{
		if (2 == res1)
		{
			split1p1edge(ts_1p1, A, B, C, intersec);
			return 1;
		}
		else if (3 == res1)
		{
			split1p1edge(ts_1p1, B, C, A, intersec);
			return 1;
		}
		else if (4 == res1)
		{
			split1p1edge(ts_1p1, C, A, B, intersec);
			return 1;
		}
		else
		{
			printf("Two masks that intersect on one edge have an error.\n");
			return 0;
		}

	}
	else
	{
		printf("1,1---Endpoint position abnormal.");
		return 0;
	}
}

void ImgProc::GetAdjacFarPs(vector<Point2d> &adjacentPoints, vector<Point2d> &farSet
	, vector<LineSegIntersec> &lSegIntersec, Point2d intersec)
{
	adjacentPoints.clear();
	farSet.clear();
	if (norm(intersec - lSegIntersec[0].ps[0])
		< norm(intersec - lSegIntersec[0].ps[1]))
	{
		adjacentPoints.push_back(lSegIntersec[0].ps[0]);
		farSet.push_back(lSegIntersec[0].ps[1]);
	}
	else
	{
		adjacentPoints.push_back(lSegIntersec[0].ps[1]);
		farSet.push_back(lSegIntersec[0].ps[0]);
	}

	if (norm(intersec - lSegIntersec[1].ps[0])
		< norm(intersec - lSegIntersec[1].ps[1]))
	{
		adjacentPoints.push_back(lSegIntersec[1].ps[0]);
		farSet.push_back(lSegIntersec[1].ps[1]);
	}
	else
	{
		adjacentPoints.push_back(lSegIntersec[1].ps[1]);
		farSet.push_back(lSegIntersec[1].ps[0]);
	}
}

int ImgProc::Triangula4Edge2_2p3(vector<vector<Point2d>> &addTs
	, vector<LineSegIntersec> &lSegIntersec
	, Point2d A, Point2d B, Point2d C)
{
	Point2d intersec;
	int res1, res2;

	if (!GetIntersecParams(intersec, res1, res2, lSegIntersec))
	{
		printf("Error(2,3): intersecting params");
		return 0;
	}

	if (IsVertexCoincident(intersec, A, B, C))
	{
		if (3 == lSegIntersec[0].ps.size())
		{
			if (areCollinear(lSegIntersec[0].p1, lSegIntersec[0].p2, A, B))
			{
				printf("(2,3): Edge coincidence-A,B");
				return 1;
			}
			else if (areCollinear(lSegIntersec[0].p1, lSegIntersec[0].p2, B, C))
			{
				printf("(2,3): Edge coincidence-B,C");
				return 1;
			}
			else if (areCollinear(lSegIntersec[0].p1, lSegIntersec[0].p2, C, A))
			{
				printf("(2,3): Edge coincidence-C,A");
				return 1;
			}
		}
		else
		{
			if (areCollinear(lSegIntersec[1].p1, lSegIntersec[1].p2, A, B))
			{
				printf("(2,3): Edge coincidence-A,B");
				return 1;
			}
			else if (areCollinear(lSegIntersec[1].p1, lSegIntersec[1].p2, B, C))
			{
				printf("(2,3): Edge coincidence-B,C");
				return 1;
			}
			else if (areCollinear(lSegIntersec[1].p1, lSegIntersec[1].p2, C, A))
			{
				printf("(2,3): Edge coincidence-C,A");
				return 1;
			}
		}

		Point2d ab;
		double maxDist = 0;

		if (3 == lSegIntersec[0].ps.size())
		{
			for (int i = 0; i < lSegIntersec[0].ps.size(); ++i)
			{
				double d = norm(lSegIntersec[0].ps[i] - intersec);
				if (maxDist < d)
				{
					maxDist = d;
					ab = lSegIntersec[0].ps[i];
				}
			}
		}
		else
		{
			for (int i = 0; i < lSegIntersec[1].ps.size(); ++i)
			{
				double d = norm(lSegIntersec[1].ps[i] - intersec);
				if (maxDist < d)
				{
					maxDist = d;
					ab = lSegIntersec[1].ps[i];
				}
			}
		}

		if (norm(C - intersec) < EPSILON)
		{
			splitEdge1(addTs, A, B, C, ab);
			return 1;
		}
		else if (norm(A - intersec) < EPSILON)
		{
			splitEdge1(addTs, B, C, A, ab);
			return 1;
		}
		else if (norm(B - intersec) < EPSILON)
		{
			splitEdge1(addTs, C, A, B, ab);
			return 1;
		}
		else
		{
			printf("Error(2,3): the intersection in edge");
			return 0;
		}
	}
	else
	{
		if (1 < res1 && res1 == res2)
		{
			int idx = -1;
			if (2 == lSegIntersec[0].ps.size())
			{
				idx = 0;
			}
			else if (2 == lSegIntersec[1].ps.size())
			{
				idx = 1;
			}
			if (areCollinear(lSegIntersec[idx].p1, lSegIntersec[idx].p2, A, B))
			{
				splitEdge1(addTs, A, B, C, intersec);
				return 1;
			}
			else if (areCollinear(lSegIntersec[idx].p1, lSegIntersec[idx].p2, B, C))
			{
				splitEdge1(addTs, B, C, A, intersec);
				return 1;
			}
			else if (areCollinear(lSegIntersec[idx].p1, lSegIntersec[idx].p2, C, A))
			{
				splitEdge1(addTs, C, A, B, intersec);
				return 1;
			}
			else
			{
				if (isPointOnSegment(A, lSegIntersec[1 - idx].p1, lSegIntersec[1 - idx].p2))
				{
					if (isPointOnSegment(lSegIntersec[idx].ps[0], B, C)
						&& isPointOnSegment(lSegIntersec[idx].ps[1], C, A))
					{
						split2p3_bc(addTs
							, B, C, A, intersec
							, lSegIntersec[idx].ps[1]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[1], B, C)
						&& isPointOnSegment(lSegIntersec[idx].ps[0], C, A))
					{
						split2p3_bc(addTs
							, B, C, A, intersec
							, lSegIntersec[idx].ps[0]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[0], B, C)
						&& isPointOnSegment(lSegIntersec[idx].ps[1], A, B))
					{
						split2p3_ca(addTs
							, B, C, A, intersec
							, lSegIntersec[idx].ps[1]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[1], B, C)
						&& isPointOnSegment(lSegIntersec[idx].ps[0], A, B))
					{
						split2p3_ca(addTs
							, B, C, A, intersec
							, lSegIntersec[idx].ps[0]);
						return 1;
					}
					else
					{
						printf("Error(2,3): intersec is not on BC edge.");
						return 0;
					}
				}
				else if (isPointOnSegment(B, lSegIntersec[1 - idx].p1, lSegIntersec[1 - idx].p2))
				{
					if (isPointOnSegment(lSegIntersec[idx].ps[0], C, A)
						&& isPointOnSegment(lSegIntersec[idx].ps[1], A, B))
					{
						split2p3_bc(addTs
							, C, A, B, intersec
							, lSegIntersec[idx].ps[1]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[1], C, A)
						&& isPointOnSegment(lSegIntersec[idx].ps[0], A, B))
					{
						split2p3_bc(addTs
							, C, A, B, intersec
							, lSegIntersec[idx].ps[0]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[0], C, A)
						&& isPointOnSegment(lSegIntersec[idx].ps[1], B, C))
					{
						split2p3_ca(addTs
							, C, A, B, intersec
							, lSegIntersec[idx].ps[1]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[1], C, A)
						&& isPointOnSegment(lSegIntersec[idx].ps[0], B, C))
					{
						split2p3_ca(addTs
							, C, A, B, intersec
							, lSegIntersec[idx].ps[0]);
						return 1;
					}
					else
					{
						printf("Error(2,3): intersec is not on CA edge.");
						return 0;
					}
				}
				else if (isPointOnSegment(C, lSegIntersec[1 - idx].p1, lSegIntersec[1 - idx].p2))
				{
					if (isPointOnSegment(lSegIntersec[idx].ps[0], A, B)
						&& isPointOnSegment(lSegIntersec[idx].ps[1], B, C))
					{
						split2p3_bc(addTs
							, A, B, C, intersec
							, lSegIntersec[idx].ps[1]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[1], A, B)
						&& isPointOnSegment(lSegIntersec[idx].ps[0], B, C))
					{
						split2p3_bc(addTs
							, A, B, C, intersec
							, lSegIntersec[idx].ps[0]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[0], A, B)
						&& isPointOnSegment(lSegIntersec[idx].ps[1], C, A))
					{
						split2p3_ca(addTs
							, A, B, C, intersec
							, lSegIntersec[idx].ps[1]);
						return 1;
					}
					else if (isPointOnSegment(lSegIntersec[idx].ps[1], A, B)
						&& isPointOnSegment(lSegIntersec[idx].ps[0], C, A))
					{
						split2p3_ca(addTs
							, A, B, C, intersec
							, lSegIntersec[idx].ps[0]);
						return 1;
					}
					else
					{
						printf("Error(2,3): intersec is not on AB edge.");
						return 0;
					}
				}
				else
				{
					printf("Error(2,3): vertex is not on edge which has three intersects.");
					return 0;
				}

			}
		}
		else
		{
			printf("Error(2,3): two intersections");
			return 0;
		}
	}

}

int ImgProc::Triangula4Edge2_1p3(vector<vector<Point2d>> &addTs
	, vector<LineSegIntersec> &lSegIntersec
	, Point2d A, Point2d B, Point2d C)
{
	Point2d intersec;
	int res1, res2;

	if (!GetIntersecParams(intersec, res1, res2, lSegIntersec))
	{
		printf("Error(1,3): intersecting params");
		return 0;
	}

	if (1 < res1 && res1 == res2)
	{
		if ((lSegIntersec[0].mask & 2) && (lSegIntersec[1].mask & 2))
		{
			splitEdge1(addTs, A, B, C, intersec);
			return 1;
		}
		else if ((lSegIntersec[0].mask & 4) && (lSegIntersec[1].mask & 4))
		{
			splitEdge1(addTs, B, C, A, intersec);
			return 1;
		}
		else if ((lSegIntersec[0].mask & 8) && (lSegIntersec[1].mask & 8))
		{
			splitEdge1(addTs, C, A, B, intersec);
			return 1;
		}
		else
		{
			printf("Error(1,3): intersec is not on these edges.");
			return 0;
		}
	}
	else
	{
		printf("Error(1,3): two intersections");
		return 0;
	}
}

int ImgProc::Triangula4Edge2_2p2(vector<vector<Point2d>> &addTs
	, vector<LineSegIntersec> &lSegIntersec
	, Point2d A, Point2d B, Point2d C)
{
	Point2d intersec;
	int res1, res2;

	if (!GetIntersecParams(intersec, res1, res2, lSegIntersec))
	{
		printf("Error(2,2): intersecting params");
		return 0;
	}

	if (1 == res1 && 1 == res2)
	{
		split2p2interior(addTs, A, B, C, intersec);
		return 1;
	}
	else if (0 == res1 && 0 == res2)
	{
		if (lSegIntersec[0].mask == lSegIntersec[1].mask)
		{
			Point2d intersecNearAB1;
			Point2d intersecNearAB2;
			Point2d intersecNearBC1;
			Point2d intersecNearBC2;
			Point2d intersecNearCA1;
			Point2d intersecNearCA2;
			if (6 == lSegIntersec[0].mask)
			{
				if (norm(A - lSegIntersec[0].ps[0])
					< norm(A - lSegIntersec[1].ps[0]))
				{
					intersecNearAB1 = lSegIntersec[0].ps[0];
					intersecNearBC1 = lSegIntersec[0].ps[1];
					intersecNearAB2 = lSegIntersec[1].ps[0];
					intersecNearBC2 = lSegIntersec[1].ps[1];
				}
				else
				{
					intersecNearAB1 = lSegIntersec[1].ps[0];
					intersecNearBC1 = lSegIntersec[1].ps[1];
					intersecNearAB2 = lSegIntersec[0].ps[0];
					intersecNearBC2 = lSegIntersec[0].ps[1];
				}
				if (norm(lSegIntersec[0].ps[0] - lSegIntersec[1].ps[0])
					< norm(lSegIntersec[0].ps[1] - lSegIntersec[1].ps[1]))
				{
					split2p2outer1_l(addTs
						, A, B, C, intersecNearAB1
						, intersecNearAB2, intersecNearBC1
						, intersecNearBC2);
					return 1;
				}
				else if (norm(lSegIntersec[0].ps[0] - lSegIntersec[1].ps[0])
					> norm(lSegIntersec[0].ps[1] - lSegIntersec[1].ps[1]))
				{
					split2p2outer1_r(addTs
						, A, B, C, intersecNearAB1
						, intersecNearAB2, intersecNearBC1
						, intersecNearBC2);
					return 1;
				}
				else
				{
					printf("Error(2,2): The projections of two sides are equal.");
					return 0;
				}

			}
			else if (12 == lSegIntersec[0].mask)
			{
				if (norm(B - lSegIntersec[0].ps[0])
					< norm(B - lSegIntersec[1].ps[0]))
				{
					intersecNearBC1 = lSegIntersec[0].ps[0];
					intersecNearCA1 = lSegIntersec[0].ps[1];
					intersecNearBC2 = lSegIntersec[1].ps[0];
					intersecNearCA2 = lSegIntersec[1].ps[1];
				}
				else
				{
					intersecNearBC1 = lSegIntersec[1].ps[0];
					intersecNearCA1 = lSegIntersec[1].ps[1];
					intersecNearBC2 = lSegIntersec[0].ps[0];
					intersecNearCA2 = lSegIntersec[0].ps[1];
				}

				if (norm(lSegIntersec[0].ps[0] - lSegIntersec[1].ps[0])
					< norm(lSegIntersec[0].ps[1] - lSegIntersec[1].ps[1]))
				{
					split2p2outer1_l(addTs
						, B, C, A, intersecNearBC1
						, intersecNearBC2, intersecNearCA1
						, intersecNearCA2);
					return 1;
				}
				else  if (norm(lSegIntersec[0].ps[0] - lSegIntersec[1].ps[0])
					> norm(lSegIntersec[0].ps[1] - lSegIntersec[1].ps[1]))
				{
					split2p2outer1_r(addTs
						, B, C, A, intersecNearBC1
						, intersecNearBC2, intersecNearCA1
						, intersecNearCA2);
					return 1;
				}
				else
				{
					printf("Error(2,2): The projections of two sides are equal.");
					return 0;
				}
			}
			else if (10 == lSegIntersec[0].mask)
			{
				if (norm(C - lSegIntersec[0].ps[0])
					< norm(C - lSegIntersec[1].ps[0]))
				{
					intersecNearCA1 = lSegIntersec[0].ps[1];
					intersecNearAB1 = lSegIntersec[0].ps[0];
					intersecNearCA2 = lSegIntersec[1].ps[1];
					intersecNearAB2 = lSegIntersec[1].ps[0];
				}
				else
				{
					intersecNearCA1 = lSegIntersec[1].ps[1];
					intersecNearAB1 = lSegIntersec[1].ps[0];
					intersecNearCA2 = lSegIntersec[0].ps[1];
					intersecNearAB2 = lSegIntersec[0].ps[0];
				}

				if (norm(lSegIntersec[0].ps[1] - lSegIntersec[1].ps[1])
					< norm(lSegIntersec[0].ps[0] - lSegIntersec[1].ps[0]))
				{
					split2p2outer1_l(addTs
						, C, A, B, intersecNearCA1
						, intersecNearCA2, intersecNearAB1
						, intersecNearAB2);
					return 1;
				}
				else  if (norm(lSegIntersec[0].ps[1] - lSegIntersec[1].ps[1])
					> norm(lSegIntersec[0].ps[0] - lSegIntersec[1].ps[0]))
				{
					split2p2outer1_l(addTs
						, C, A, B, intersecNearCA1
						, intersecNearCA2, intersecNearAB1
						, intersecNearAB2);
					return 1;
				}
				else
				{
					printf("Error(2,2): The projections of two sides are equal.");
					return 0;
				}
			}
			else
			{
				printf("Error(2,2): intersecting params");
				return 0;
			}

		}
		else
		{
			if (norm(lSegIntersec[0].ps[0] - lSegIntersec[0].ps[1]) < EPSILON
				&& norm(lSegIntersec[1].ps[0] - lSegIntersec[1].ps[1]) < EPSILON)
			{
				printf("Intersects only with two vertices.");
				return 1;
			}
			else if (norm(lSegIntersec[0].ps[0] - lSegIntersec[0].ps[1]) < EPSILON)
			{
				if ((lSegIntersec[1].mask & 2) && (lSegIntersec[1].mask & 4))
				{
					split2(addTs, A, B, C, lSegIntersec[1].ps[0], lSegIntersec[1].ps[1]);
					return 1;
				}
				else if ((lSegIntersec[1].mask & 4) && (lSegIntersec[1].mask & 8))
				{
					split2(addTs, B, C, A, lSegIntersec[1].ps[0], lSegIntersec[1].ps[1]);
					return 1;
				}
				else if ((lSegIntersec[1].mask & 2) && (lSegIntersec[1].mask & 8))
				{
					split2(addTs, C, A, B, lSegIntersec[1].ps[1], lSegIntersec[1].ps[0]);
					return 1;
				}
				else
				{
					printf("Error:single segment's two endpoints.\n");
					return 0;
				}
			}
			else if (norm(lSegIntersec[1].ps[0] - lSegIntersec[1].ps[1]) < EPSILON)
			{
				if ((lSegIntersec[0].mask & 2) && (lSegIntersec[0].mask & 4))
				{
					split2(addTs, A, B, C, lSegIntersec[0].ps[0], lSegIntersec[0].ps[1]);
					return 1;
				}
				else if ((lSegIntersec[0].mask & 4) && (lSegIntersec[0].mask & 8))
				{
					split2(addTs, B, C, A, lSegIntersec[0].ps[0], lSegIntersec[0].ps[1]);
					return 1;
				}
				else if ((lSegIntersec[0].mask & 2) && (lSegIntersec[0].mask & 8))
				{
					split2(addTs, C, A, B, lSegIntersec[0].ps[1], lSegIntersec[0].ps[0]);
					return 1;
				}
				else
				{
					printf("Error:single segment's two endpoints.\n");
					return 0;
				}
			}
			else
			{
				if ((lSegIntersec[0].mask & 2) && (lSegIntersec[1].mask & 2))
				{
					Point2d nearA, nearB, ca, bc;

					if ((lSegIntersec[0].mask & 8) && (lSegIntersec[1].mask & 4))
					{
						nearA = lSegIntersec[0].ps[0];
						nearB = lSegIntersec[1].ps[0];
						ca = lSegIntersec[0].ps[1];
						bc = lSegIntersec[1].ps[1];
					}
					else if ((lSegIntersec[1].mask & 8) && (lSegIntersec[0].mask & 4))
					{
						nearA = lSegIntersec[1].ps[0];
						nearB = lSegIntersec[0].ps[0];
						ca = lSegIntersec[1].ps[1];
						bc = lSegIntersec[0].ps[1];
					}
					else
					{
						printf("Error:two edge segments intersect at the same triangle edge.\n");
						return 0;
					}

					split2p2outer2(addTs, A, B, C, nearA, nearB, ca, bc);
					return 1;
				}
				else if ((lSegIntersec[0].mask & 4) && (lSegIntersec[1].mask & 4))
				{
					Point2d nearB, nearC, ab, ca;
					if ((lSegIntersec[0].mask & 2) && (lSegIntersec[1].mask & 8))
					{
						nearB = lSegIntersec[0].ps[1];
						nearC = lSegIntersec[1].ps[0];
						ab = lSegIntersec[0].ps[0];
						ca = lSegIntersec[1].ps[1];
					}
					else if ((lSegIntersec[1].mask & 2) && (lSegIntersec[0].mask & 8))
					{
						nearB = lSegIntersec[1].ps[1];
						nearC = lSegIntersec[0].ps[0];
						ab = lSegIntersec[1].ps[0];
						ca = lSegIntersec[0].ps[1];
					}
					else
					{
						printf("Error:two edge segments intersect at the same triangle edge.\n");
						return 0;
					}

					split2p2outer2(addTs, B, C, A, nearB, nearC, ab, ca);
					return 1;
				}
				else if ((lSegIntersec[0].mask & 8) && (lSegIntersec[1].mask & 8))
				{
					Point2d nearC, nearA, bc, ab;
					if ((lSegIntersec[0].mask & 4) && (lSegIntersec[1].mask & 2))
					{
						nearC = lSegIntersec[0].ps[1];
						nearA = lSegIntersec[1].ps[1];
						bc = lSegIntersec[0].ps[0];
						ab = lSegIntersec[1].ps[0];
					}
					else if ((lSegIntersec[1].mask & 4) && (lSegIntersec[0].mask & 2))
					{
						nearC = lSegIntersec[1].ps[1];
						nearA = lSegIntersec[0].ps[1];
						bc = lSegIntersec[1].ps[0];
						ab = lSegIntersec[0].ps[0];
					}
					else
					{
						printf("Error:two edge segments intersect at the same triangle edge.\n");
						return 0;
					}

					split2p2outer2(addTs, C, A, B, nearC, nearA, bc, ab);
					return 1;
				}
				else
				{
					printf("Error(2,2 outer):\n");
					printf("The case of intersecting with two different edges.\n");
					return 0;
				}
			}
		}
	}
	else if (1 < res1 && res1 == res2)
	{
		bool resAB = (lSegIntersec[0].mask & 2) && (lSegIntersec[1].mask & 2);
		bool resBC = (lSegIntersec[0].mask & 4) && (lSegIntersec[1].mask & 4);
		bool resCA = (lSegIntersec[0].mask & 8) && (lSegIntersec[1].mask & 8);
		Point2d fastIntersec1, fastIntersec2;
		if (norm(intersec - lSegIntersec[0].ps[0]) < TOLERANCE)
		{
			fastIntersec1 = lSegIntersec[0].ps[1];
		}
		else
		{
			fastIntersec1 = lSegIntersec[0].ps[0];
		}

		if (norm(intersec - lSegIntersec[1].ps[0]) < TOLERANCE)
		{
			fastIntersec2 = lSegIntersec[1].ps[1];
		}
		else
		{
			fastIntersec2 = lSegIntersec[1].ps[0];
		}

		if (IsVertexCoincident(fastIntersec1, A, B, C)
			&& IsVertexCoincident(fastIntersec2, A, B, C))
		{
			if (resAB)
			{
				splitEdge1(addTs, A, B, C, intersec);
				return 1;
			}
			else if (resBC)
			{
				splitEdge1(addTs, B, C, A, intersec);
				return 1;
			}
			else if (resCA)
			{
				splitEdge1(addTs, C, A, B, intersec);
				return 1;
			}
			else
			{
				printf("Error(2,2 edge two vertex)\n");
				return 0;
			}

		}
		else if (IsVertexCoincident(fastIntersec1, A, B, C))
		{
			if (resAB)
			{
				if (lSegIntersec[1].mask & 4)
				{
					split2p2edge_Vt1_bc(addTs, A, B, C, intersec, fastIntersec2);
					return 1;
				}
				else if (lSegIntersec[1].mask & 8)
				{
					split2p2edge_Vt1_ca(addTs, A, B, C, intersec, fastIntersec2);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge and one vertex in AB)\n");
					return 0;
				}
			}
			else if (resBC)
			{
				if (lSegIntersec[1].mask & 8)
				{
					split2p2edge_Vt1_bc(addTs, B, C, A, intersec, fastIntersec2);
					return 1;
				}
				else if (lSegIntersec[1].mask & 2)
				{
					split2p2edge_Vt1_ca(addTs, B, C, A, intersec, fastIntersec2);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge and one vertex in BC)\n");
					return 0;
				}
			}
			else if (resCA)
			{
				if (lSegIntersec[1].mask & 2)
				{
					split2p2edge_Vt1_bc(addTs, C, A, B, intersec, fastIntersec2);
					return 1;
				}
				else if (lSegIntersec[1].mask & 4)
				{
					split2p2edge_Vt1_ca(addTs, C, A, B, intersec, fastIntersec2);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge and one vertex in CA)\n");
					return 0;
				}
			}
			else
			{
				printf("Error(2,2 edge and one vertex in fastIntersec2)\n");
				return 0;
			}
		}
		else if (IsVertexCoincident(fastIntersec2, A, B, C))
		{
			if (resAB)
			{
				if (lSegIntersec[0].mask & 4)
				{
					split2p2edge_Vt1_bc(addTs, A, B, C, intersec, fastIntersec1);
					return 1;
				}
				else if (lSegIntersec[0].mask & 8)
				{
					split2p2edge_Vt1_ca(addTs, A, B, C, intersec, fastIntersec1);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge and one vertex in AB)\n");
					return 0;
				}
			}
			else if (resBC)
			{
				if (lSegIntersec[0].mask & 8)
				{
					split2p2edge_Vt1_bc(addTs, B, C, A, intersec, fastIntersec1);
					return 1;
				}
				else if (lSegIntersec[0].mask & 2)
				{
					split2p2edge_Vt1_ca(addTs, B, C, A, intersec, fastIntersec1);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge and one vertex in BC)\n");
					return 0;
				}
			}
			else if (resCA)
			{
				if (lSegIntersec[0].mask & 2)
				{
					split2p2edge_Vt1_bc(addTs, C, A, B, intersec, fastIntersec1);
					return 1;
				}
				else if (lSegIntersec[0].mask & 4)
				{
					split2p2edge_Vt1_ca(addTs, C, A, B, intersec, fastIntersec1);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge and one vertex in CA)\n");
					return 0;
				}
			}
			else
			{
				printf("Error(2,2 edge and one vertex in fastIntersec1)\n");
				return 0;
			}
		}
		else if (lSegIntersec[0].mask == lSegIntersec[1].mask)
		{
			Point2d intersecNearA;
			Point2d intersecNearB;
			Point2d intersecNearC;
			if (resAB)
			{
				if (lSegIntersec[0].mask & 4)
				{
					if (norm(B - lSegIntersec[0].ps[1])
						< norm(B - lSegIntersec[1].ps[1]))
					{
						intersecNearB = lSegIntersec[0].ps[1];
						intersecNearC = lSegIntersec[1].ps[1];
					}
					else if (norm(B - lSegIntersec[0].ps[1])
						> norm(B - lSegIntersec[1].ps[1]))
					{
						intersecNearB = lSegIntersec[1].ps[1];
						intersecNearC = lSegIntersec[0].ps[1];
					}
					else
					{
						printf("Error(2,2):distance in BC\n");
						return 0;
					}
					split2p2edge2_bc(addTs, A, B, C, intersec
						, intersecNearB, intersecNearC);
					return 1;
				}
				else if (lSegIntersec[0].mask & 8)
				{
					if (norm(C - lSegIntersec[0].ps[1])
						< norm(C - lSegIntersec[1].ps[1]))
					{
						intersecNearC = lSegIntersec[0].ps[1];
						intersecNearA = lSegIntersec[1].ps[1];
					}
					else if (norm(C - lSegIntersec[0].ps[1])
						> norm(C - lSegIntersec[1].ps[1]))
					{
						intersecNearC = lSegIntersec[1].ps[1];
						intersecNearA = lSegIntersec[0].ps[1];
					}
					else
					{
						printf("Error(2,2):distance in CA\n");
						return 0;
					}
					split2p2edge2_ca(addTs, A, B, C, intersec
						, intersecNearC, intersecNearA);
					return 1;
				}
				else
				{
					printf("Error(2,2):intersections(one point in edge and one edge)\n");
					return 0;
				}
			}
			else if (resBC)
			{
				if (lSegIntersec[0].mask & 8)
				{
					if (norm(C - lSegIntersec[0].ps[1])
						< norm(C - lSegIntersec[1].ps[1]))
					{
						intersecNearC = lSegIntersec[0].ps[1];
						intersecNearA = lSegIntersec[1].ps[1];
					}
					else if (norm(C - lSegIntersec[0].ps[1])
						> norm(C - lSegIntersec[1].ps[1]))
					{
						intersecNearC = lSegIntersec[1].ps[1];
						intersecNearA = lSegIntersec[0].ps[1];
					}
					else
					{
						printf("Error(2,2):distance in CA\n");
						return 0;
					}
					split2p2edge2_bc(addTs, B, C, A, intersec
						, intersecNearC, intersecNearA);
					return 1;
				}
				else if (lSegIntersec[0].mask & 2)
				{
					if (norm(A - lSegIntersec[0].ps[0])
						< norm(A - lSegIntersec[1].ps[0]))
					{
						intersecNearA = lSegIntersec[0].ps[0];
						intersecNearB = lSegIntersec[1].ps[0];
					}
					else if (norm(A - lSegIntersec[0].ps[0])
						> norm(A - lSegIntersec[1].ps[0]))
					{
						intersecNearA = lSegIntersec[1].ps[0];
						intersecNearB = lSegIntersec[0].ps[0];
					}
					else
					{
						printf("Error(2,2):distance in CA\n");
						return 0;
					}
					split2p2edge2_ca(addTs, B, C, A, intersec
						, intersecNearA, intersecNearB);
					return 1;
				}
				else
				{
					printf("Error(2,2):intersections(one point in edge and one edge)\n");
					return 0;
				}
			}
			else if (resCA)
			{
				if (lSegIntersec[0].mask & 2)
				{
					if (norm(A - lSegIntersec[0].ps[0])
						< norm(A - lSegIntersec[1].ps[0]))
					{
						intersecNearA = lSegIntersec[0].ps[0];
						intersecNearB = lSegIntersec[1].ps[0];
					}
					else if (norm(A - lSegIntersec[0].ps[0])
						> norm(A - lSegIntersec[1].ps[0]))
					{
						intersecNearA = lSegIntersec[1].ps[0];
						intersecNearB = lSegIntersec[0].ps[0];
					}
					else
					{
						printf("Error(2,2):distance in CA\n");
						return 0;
					}
					split2p2edge2_bc(addTs, C, A, B, intersec
						, intersecNearA, intersecNearB);
					return 1;
				}
				else if (lSegIntersec[0].mask & 4)
				{
					if (norm(B - lSegIntersec[0].ps[0])
						< norm(B - lSegIntersec[1].ps[0]))
					{
						intersecNearB = lSegIntersec[0].ps[0];
						intersecNearC = lSegIntersec[1].ps[0];
					}
					else if (norm(B - lSegIntersec[0].ps[0])
						> norm(B - lSegIntersec[1].ps[0]))
					{
						intersecNearB = lSegIntersec[1].ps[0];
						intersecNearC = lSegIntersec[0].ps[0];
					}
					else
					{
						printf("Error(2,2):distance in CA\n");
						return 0;
					}
					split2p2edge2_ca(addTs, C, A, B, intersec
						, intersecNearB, intersecNearC);
					return 1;
				}
				else
				{
					printf("Error(2,2):intersections(one point in edge and one edge)\n");
					return 0;
				}
			}
			else
			{
				printf("Error(2,2):The intersection point is not on the edges\n");
				return 0;
			}
		}
		else
		{
			if (resAB)
			{
				if ((lSegIntersec[0].mask & 8)
					&& (lSegIntersec[1].mask & 4))
				{
					split2p2edge_noVt(addTs
						, A, B, C, intersec
						, fastIntersec1, fastIntersec2);
					return 1;
				}
				else if ((lSegIntersec[0].mask & 4)
					&& (lSegIntersec[1].mask & 8))
				{
					split2p2edge_noVt(addTs
						, A, B, C, intersec
						, fastIntersec2, fastIntersec1);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge no vertex in AB)\n");
					return 0;
				}
			}
			else if (resBC)
			{
				if ((lSegIntersec[0].mask & 2)
					&& (lSegIntersec[1].mask & 8))
				{
					split2p2edge_noVt(addTs
						, B, C, A, intersec
						, fastIntersec1, fastIntersec2);
					return 1;
				}
				else if ((lSegIntersec[0].mask & 8)
					&& (lSegIntersec[1].mask & 2))
				{
					split2p2edge_noVt(addTs
						, B, C, A, intersec
						, fastIntersec2, fastIntersec1);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge no vertex in BC)\n");
					return 0;
				}
			}
			else if (resCA)
			{
				if ((lSegIntersec[0].mask & 4)
					&& (lSegIntersec[1].mask & 2))
				{
					split2p2edge_noVt(addTs
						, C, A, B, intersec
						, fastIntersec1, fastIntersec2);
					return 1;
				}
				else if ((lSegIntersec[0].mask & 2)
					&& (lSegIntersec[1].mask & 4))
				{
					split2p2edge_noVt(addTs
						, C, A, B, intersec
						, fastIntersec2, fastIntersec1);
					return 1;
				}
				else
				{
					printf("Error(2,2 edge no vertex in CA)\n");
					return 0;
				}
			}
			else
			{
				printf("Error(2,2 edge no vertex)\n");
				return 0;
			}
		}
	}
	else
	{
		printf("Error(2,2) res1 %d is different from res2 %d\n", res1, res2);
		return 0;
	}
}

bool ImgProc::IsVertexCoincident(Point2d p, Point2d A, Point2d B, Point2d C)
{
	if (norm(p - A) < TOLERANCE || norm(p - B) < TOLERANCE || norm(p - C) < TOLERANCE)
	{
		return true;
	}
	else
	{
		return false;
	}
}

int ImgProc::Triangula4Edge2_1p2(vector<vector<Point2d>> &addTs
	, vector<LineSegIntersec> &lSegIntersec
	, Point2d A, Point2d B, Point2d C)
{
	Point2d intersec;
	int res1, res2;

	if (!GetIntersecParams(intersec, res1, res2, lSegIntersec))
	{
		printf("Error: intersecting params");
		return 0;
	}

	if (1 == res1 && 1 == res2)
	{
		Point2d intersecNearEdge;

		if (2 == lSegIntersec[0].ps.size()
			&& 1 == lSegIntersec[1].ps.size())
		{
			if (arePointsEqual(lSegIntersec[0].ps[0], lSegIntersec[0].ps[1]))
			{
				if (arePointsEqual(lSegIntersec[0].ps[0], A))
				{
					if (4 == lSegIntersec[1].mask)
					{
						split1p2interior2(addTs
							, A, B, C, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else if (2 == lSegIntersec[1].mask)
					{
						split1p2interior1_ab(addTs
							, A, B, C, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else if (8 == lSegIntersec[1].mask)
					{
						split1p2interior1_ca(addTs
							, A, B, C, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else
					{
						printf("The case of the interior point \
							   						at (1,2) intersecting with A has an error.\n");
						return 0;
					}
				}
				else if (arePointsEqual(lSegIntersec[0].ps[0], B))
				{
					if (8 == lSegIntersec[1].mask)
					{
						split1p2interior2(addTs
							, B, C, A, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else if (4 == lSegIntersec[1].mask)
					{
						split1p2interior1_ab(addTs
							, B, C, A, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else if (2 == lSegIntersec[1].mask)
					{
						split1p2interior1_ca(addTs
							, B, C, A, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else
					{
						printf("The case of the interior point \
							   						at (1,2) intersecting with B has an error.\n");
						return 0;
					}
				}
				else if (arePointsEqual(lSegIntersec[0].ps[0], C))
				{
					if (2 == lSegIntersec[1].mask)
					{
						split1p2interior2(addTs
							, C, A, B, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else if (8 == lSegIntersec[1].mask)
					{
						split1p2interior1_ab(addTs
							, C, A, B, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else if (4 == lSegIntersec[1].mask)
					{
						split1p2interior1_ca(addTs
							, C, A, B, intersec
							, lSegIntersec[1].ps[0]);
						return 1;
					}
					else
					{
						printf("The case of the interior point ");
						printf("at(1, 2) intersecting with C has an error.\n");
						return 0;
					}
				}
				else
				{
					printf("Error(1,2 inner):does not pass through the vertex.\n");
					return 0;
				}

			}
			else
			{
				printf("Error: The two intersection points ");
				printf("of intersection line 1 do not coincide.\n");
				return 0;
			}

		}
		else if (2 == lSegIntersec[1].ps.size()
			&& 1 == lSegIntersec[0].ps.size())
		{
			if (arePointsEqual(lSegIntersec[1].ps[0], lSegIntersec[1].ps[1]))
			{
				if (arePointsEqual(lSegIntersec[1].ps[0], A))
				{
					if (4 == lSegIntersec[0].mask)
					{
						split1p2interior2(addTs
							, A, B, C, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else if (2 == lSegIntersec[0].mask)
					{
						split1p2interior1_ab(addTs
							, A, B, C, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else if (8 == lSegIntersec[0].mask)
					{
						split1p2interior1_ca(addTs
							, A, B, C, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else
					{
						printf("The case of the interior point \
							   							  at (1,2) intersecting with A has an error.\n");
						return 0;
					}
				}
				else if (arePointsEqual(lSegIntersec[1].ps[0], B))
				{
					if (8 == lSegIntersec[0].mask)
					{
						split1p2interior2(addTs
							, B, C, A, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else if (4 == lSegIntersec[0].mask)
					{
						split1p2interior1_ab(addTs
							, B, C, A, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else if (2 == lSegIntersec[0].mask)
					{
						split1p2interior1_ca(addTs
							, B, C, A, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else
					{
						printf("The case of the interior point \
							   							   at (1,2) intersecting with B has an error.\n");
						return 0;
					}
				}
				else if (arePointsEqual(lSegIntersec[1].ps[0], C))
				{
					if (2 == lSegIntersec[0].mask)
					{
						split1p2interior2(addTs
							, C, A, B, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else if (8 == lSegIntersec[0].mask)
					{
						split1p2interior1_ab(addTs
							, C, A, B, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else if (4 == lSegIntersec[0].mask)
					{
						split1p2interior1_ca(addTs
							, C, A, B, intersec
							, lSegIntersec[0].ps[0]);
						return 1;
					}
					else
					{
						printf("The case of the interior point \
							   							  at (1,2) intersecting with C has an error.\n");
						return 0;
					}
				}
				else
				{
					printf("Error(1,2 inner):does not pass through the vertex.\n");
					return 0;
				}

			}
			else
			{
				printf("Error: The two intersection points\
					   					   of intersection line 1 do not coincide.\n");
				return 0;
			}
		}
		else
		{
			printf("Error: Does not conform to the (1,2) rule.\n");
			return 0;
		}

	}
	else if (1 < res1 && res1 == res2)
	{
		Point2d intersec2;
		if (1 == lSegIntersec[0].ps.size()
			&& 2 == lSegIntersec[1].ps.size())
		{
			double norm0 = norm(lSegIntersec[1].ps[0] - intersec);
			double norm1 = norm(lSegIntersec[1].ps[1] - intersec);
			if (norm0 < EPSILON)
			{
				intersec2 = lSegIntersec[1].ps[1];
			}
			else if (norm1 < EPSILON)
			{
				intersec2 = lSegIntersec[1].ps[0];
			}
			else
			{
				printf("Error: The second intersection point is abnormal.\n");
				return 0;
			}

			if (2 == lSegIntersec[0].mask && 2 == res1)
			{
				if (norm(intersec2 - A) < EPSILON
					|| norm(intersec2 - B) < EPSILON)
				{
					splitEdge1(addTs, A, B, C, intersec);
					return 1;
				}
				else
				{
					if (lSegIntersec[1].mask & 4)
					{
						split1p2Edge2_BC(addTs, A, B, C, intersec, intersec2);
						return 1;
					}
					else if (lSegIntersec[1].mask & 8)
					{
						split1p2Edge2_CA(addTs, A, B, C, intersec, intersec2);
						return 1;
					}
					else
					{
						printf("Error: 1,2 edge intersection case 2.\n");
						return 0;
					}
				}
			}
			else if (4 == lSegIntersec[0].mask && 3 == res1)
			{
				if (norm(intersec2 - B) < EPSILON
					|| norm(intersec2 - C) < EPSILON)
				{
					splitEdge1(addTs, B, C, A, intersec);
					return 1;
				}
				else
				{
					if (lSegIntersec[1].mask & 8)
					{
						split1p2Edge2_BC(addTs, B, C, A, intersec, intersec2);
						return 1;
					}
					else if (lSegIntersec[1].mask & 2)
					{
						split1p2Edge2_CA(addTs, B, C, A, intersec, intersec2);
						return 1;
					}
					else
					{
						printf("Error: 1,2 edge intersection case 2.\n");
						return 0;
					}
				}
			}
			else if (8 == lSegIntersec[0].mask && 4 == res1)
			{
				if (norm(intersec2 - C) < EPSILON
					|| norm(intersec2 - A) < EPSILON)
				{
					splitEdge1(addTs, C, A, B, intersec);
					return 1;
				}
				else
				{
					if (lSegIntersec[1].mask & 2)
					{
						split1p2Edge2_BC(addTs, C, A, B, intersec, intersec2);
						return 1;
					}
					else if (lSegIntersec[1].mask & 4)
					{
						split1p2Edge2_CA(addTs, C, A, B, intersec, intersec2);
						return 1;
					}
					else
					{
						printf("Error: 1,2 edge intersection case 2.\n");
						return 0;
					}
				}
			}
			else
			{
				printf("Error: Incomplete enumeration of 1,2 edge intersection cases.\n");
				return 0;
			}

		}
		else if (1 == lSegIntersec[1].ps.size()
			&& 2 == lSegIntersec[0].ps.size())
		{
			if (norm(lSegIntersec[0].ps[0] - intersec) < EPSILON)
			{
				intersec2 = lSegIntersec[0].ps[1];
			}
			else if (norm(lSegIntersec[0].ps[1] - intersec) < EPSILON)
			{
				intersec2 = lSegIntersec[0].ps[0];
			}
			else
			{
				printf("Error(1,2): The second intersection point is abnormal.\n");
				return 0;
			}

			if (2 == lSegIntersec[1].mask && 2 == res1)
			{
				if (norm(intersec2 - A) < EPSILON
					|| norm(intersec2 - B) < EPSILON)
				{
					splitEdge1(addTs, A, B, C, intersec);
					return 1;
				}
				else
				{
					if (lSegIntersec[0].mask & 4)
					{
						split1p2Edge2_BC(addTs, A, B, C, intersec, intersec2);
						return 1;
					}
					else if (lSegIntersec[0].mask & 8)
					{
						split1p2Edge2_CA(addTs, A, B, C, intersec, intersec2);
						return 1;
					}
					else
					{
						printf("Error: 1,2 edge intersection case 2.\n");
						return 0;
					}
				}
			}
			else if (4 == lSegIntersec[1].mask && 3 == res1)
			{
				if (norm(intersec2 - B) < EPSILON
					|| norm(intersec2 - C) < EPSILON)
				{
					splitEdge1(addTs, B, C, A, intersec);
					return 1;
				}
				else
				{
					if (lSegIntersec[0].mask & 8)
					{
						split1p2Edge2_BC(addTs, B, C, A, intersec, intersec2);
						return 1;
					}
					else if (lSegIntersec[0].mask & 2)
					{
						split1p2Edge2_CA(addTs, B, C, A, intersec, intersec2);
						return 1;
					}
					else
					{
						printf("Error: 1,2 edge intersection case 2.\n");
						return 0;
					}
				}
			}
			else if (8 == lSegIntersec[1].mask && 4 == res1)
			{
				if (norm(intersec2 - C) < EPSILON
					|| norm(intersec2 - A) < EPSILON)
				{
					splitEdge1(addTs, C, A, B, intersec);
					return 1;
				}
				else
				{
					if (lSegIntersec[0].mask & 2)
					{
						split1p2Edge2_BC(addTs, C, A, B, intersec, intersec2);
						return 1;
					}
					else if (lSegIntersec[0].mask & 4)
					{
						split1p2Edge2_CA(addTs, C, A, B, intersec, intersec2);
						return 1;
					}
					else
					{
						printf("Error: 1,2 edge intersection case 2.\n");
						return 0;
					}
				}
			}
			else
			{
				printf("Error: Incomplete enumeration of 1,2 edge intersection cases.\n");
				return 0;
			}

		}
		else
		{
			printf("Error(1,2): Number of intersection points\n");
			return 0;
		}
	}
	else
	{
		printf("Error(1,2): res1 and res2\n");
		return 0;
	}
}

int ImgProc::GetIntersecParams(Point2d &intersec, int &res1, int &res2
	, vector<LineSegIntersec> &lSegIntersec)
{
	Point2d ps[4];
	ps[0] = lSegIntersec[0].p1;
	ps[1] = lSegIntersec[0].p2;
	ps[2] = lSegIntersec[1].p1;
	ps[3] = lSegIntersec[1].p2;
	int idx1, idx2;

	if (false == findDuplicateIndices(ps, idx1, idx2))
	{
		printf("1,2---No intersection points.\n");
		return 0;
	}

	if (idx1 + idx2 < 2 || idx1 + idx2>4)
	{
		printf("1,2---Intersection point index error.\n");
		return 0;
	}

	if (0 == min(idx1, idx2))
	{
		res1 = lSegIntersec[0].res1;
		intersec = lSegIntersec[0].p1;
	}
	else
	{
		res1 = lSegIntersec[0].res2;
		intersec = lSegIntersec[0].p2;
	}

	if (2 == max(idx1, idx2))
	{
		res2 = lSegIntersec[1].res1;
	}
	else
	{
		res2 = lSegIntersec[1].res2;
	}
	return 1;
}

void ImgProc::DrawOrigMesh(Mat &mat, double zoom)
{
	int r = mat.rows;
	int c = mat.cols;

	int rC = r / zoom;
	int cC = c / zoom;
	Point pt1, pt2;
	int thickness = 1;
	cv::Scalar color(0, 255, 0);
	int lineType = cv::LINE_AA;
	for (int i = 1; i != rC - 1; ++i)
	{
		for (int j = 1; j != cC - 1; ++j)
		{
			pt1.x = (j - 1)*zoom;
			pt1.y = i*zoom;
			pt2.x = j*zoom;
			pt2.y = i*zoom;

			line(mat, pt1, pt2, color, thickness, lineType);

			pt1.x = (j + 1)*zoom;
			pt1.y = i*zoom;
			pt2.x = j*zoom;
			pt2.y = i*zoom;

			line(mat, pt1, pt2, color, thickness, lineType);

			pt1.x = j*zoom;
			pt1.y = (i - 1)*zoom;
			pt2.x = j*zoom;
			pt2.y = i*zoom;

			line(mat, pt1, pt2, color, thickness, lineType);

			pt1.x = j*zoom;
			pt1.y = (i + 1)*zoom;
			pt2.x = j*zoom;
			pt2.y = i*zoom;

			line(mat, pt1, pt2, color, thickness, lineType);

		}
	}
}


void ImgProc::DrawTriMesh(Mat &mat, double zoom)
{
	Point pt1, pt2;
	int thickness = 1;
	cv::Scalar color(0, 0, 255);
	int lineType = cv::LINE_AA;

	int tSize = triangle_set.size();
	int A, B, C;

	for (int i = 0; i != tSize; ++i)
	{
		A = triangle_set[i][0];
		B = triangle_set[i][1];
		C = triangle_set[i][2];

		pt1.x = ps[A].x*zoom;
		pt1.y = ps[A].y*zoom;
		pt2.x = ps[B].x*zoom;
		pt2.y = ps[B].y*zoom;

		line(mat, pt1, pt2, color, thickness, lineType);


		pt1.x = ps[B].x*zoom;
		pt1.y = ps[B].y*zoom;
		pt2.x = ps[C].x*zoom;
		pt2.y = ps[C].y*zoom;

		line(mat, pt1, pt2, color, thickness, lineType);

		pt1.x = ps[C].x*zoom;
		pt1.y = ps[C].y*zoom;
		pt2.x = ps[A].x*zoom;
		pt2.y = ps[A].y*zoom;

		line(mat, pt1, pt2, color, thickness, lineType);
	}
}

void ImgProc::GeoQualEval4Avail(StatisPara &statisPara)
{
	statisPara.minAngle = MAXV;
	statisPara.minArea = MAXV;
	statisPara.sliverNum = 0;
	statisPara.areaVari = 0;
	statisPara.equilateralCount = 0;

	int n = avail_tria.size();
	Point2d A, B, C;
	Point2d marA, marB, marC;
	Point2d manA, manB, manC;
	double averArea = 0, averAngle = 0;
	for (int i = 0; i != n; ++i)
	{
		A = ps[avail_tria[i][0]];
		B = ps[avail_tria[i][1]];
		C = ps[avail_tria[i][2]];

		if (isEquilateralTriangle(A, B, C))
		{
			statisPara.equilateralCount++;
		}

		double area = triangleArea(A, B, C);
		averArea += area / n;
		if (area < statisPara.minArea)
		{
			statisPara.minArea = area;
			marA = A;
			marB = B;
			marC = C;
		}

		double angle = triangleMinAngle(A, B, C);
		averAngle += angle / n;
		if (angle < statisPara.minAngle)
		{
			statisPara.minAngle = angle;
			manA = A;
			manB = B;
			manC = C;
		}
		if (angle < 0.08727)
		{
			statisPara.sliverNum++;
		}
	}

	for (int i = 0; i != n; ++i)
	{
		A = ps[avail_tria[i][0]];
		B = ps[avail_tria[i][1]];
		C = ps[avail_tria[i][2]];

		double area = triangleArea(A, B, C);
		statisPara.areaVari += pow(averArea - area, 2) / n;
	}

	printf("Point A of minArea: (%.6f, %.6f)\n", marA.x, marA.y);
	printf("Point B of minArea: (%.6f, %.6f)\n", marB.x, marB.y);
	printf("Point C of minArea: (%.6f, %.6f)\n", marC.x, marC.y);

	printf("Point A of minAngle: (%.6f, %.6f)\n", manA.x, manA.y);
	printf("Point B of minAngle: (%.6f, %.6f)\n", manB.x, manB.y);
	printf("Point C of minAngle: (%.6f, %.6f)\n", manC.x, manC.y);
}