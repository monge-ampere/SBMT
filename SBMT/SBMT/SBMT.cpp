// SBMT.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "ImgProc.h"
#include "FileProc.h"

typedef double double3[3];

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	ImgProc *proc = new ImgProc("data//drop.bmp"
		, "data//drop_bdy.bmp", "data//drop.bmp");
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
	write2m(ps, pN, gs, gN, tri_set, tN, "drop.m");
	write2m(ps, pN, gs, gN, tri_set, tN, "drop.open.m");

	int *bdy = new int[4096];
	int bN = 0;
	proc->OutputBdy(ps, bdy, bN);

	write2f(ps, bN, "F:\\points.csv");
	writeBdy2f(bdy, bN, "F:\\edges.csv");

	printf("finishing----------------------------\n\n");
	system("pause");
	delete proc;
	return 0;
}

