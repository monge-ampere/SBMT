// SBMT.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "ImgProc.h"
#include "FileProc.h"

typedef double double3[3];

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	ImgProc *proc = new ImgProc("data//cross-section2.bmp"
		, "data//cross-section-bdy.bmp", "E:\\testPic\\cross-section_ori.bmp");
	proc->Analysis();

	StatisPara statisPara;
	proc->GeoQualEval4Avail(statisPara);
	printf("minAngle=%.8f, minArea=%.8f------------\n", statisPara.minAngle, statisPara.minArea);
	printf("sliver number=%d------------\n", statisPara.sliverNum);
	printf("equilateral number=%d------------\n", statisPara.equilateralCount);
	printf("area variance=%f------------\n", statisPara.areaVari);
	double2 *ps = new double2[MAX_VERTEX_NUM];
	int pN = 0;
	uchar *gs = new uchar[MAX_VERTEX_NUM];
	int gN = 0;
	int3 *tri_set = new int3[MAX_TRIANGLE_NUM];
	int tN = 0;
	proc->AvailModel4Output(ps, pN, gs, gN, tri_set, tN);
	Ps2Txt(ps, pN, "E:\\testPic\\points.txt");
	Tris2Txt(tri_set, tN, "E:\\testPic\\triangles.txt");
	write2m(ps, pN, gs, gN, tri_set, tN, "drop.m");

	int *bdy = new int[4096];
	int bN = 0;
	proc->OutputBdy(ps, bdy, bN);

	//tN = 0;
	//proc->OutputTris(tri_set, tN);

	write2f(ps, bN, "E:\\testPic\\points.csv");
	writeBdy2f(bdy, bN, "E:\\testPic\\edges.csv");

	proc->OutputBdyIds(bdy, bN);
	Bdy2Txt(bdy, bN, "E:\\testPic\\boundary.txt");
	

	printf("finishing----------------------------\n\n");
	system("pause");
	delete proc;
	return 0;
}