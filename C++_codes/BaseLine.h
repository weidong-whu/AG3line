#pragma once
#include<math.h>
#include<iostream>
#include<opencv.hpp>

using namespace::cv;
using namespace::std;
class BaseLine
{
	float c1,c2,c3,c4,t;
	Mat pt,v1,v2;//for line_point distance calculation
	/*
	[c1 c2
	 c3 c4]
	*/
	float a,a_1,x1_mx,y1_my,x1_mx_y1_my,len;
	float l_min,l_max,w_min,w_max,l,w;
	float sumdx,sumdy,lambda;//line parameters
	double get_theta( Point3i *reg, int ed1,int ed2,int ed3,int ed4, double x, double y,
                         double reg_angle, double prec );
	public:
	int counter,xsize,ysize,minlength;
	float xxs[3];
	float yys[3];
	float pdense[3];
	float mx,my,dx,dy,segang,theta,sx,sy,x1,y1,x2,y2,length,midcounter,sxx,syy,ui,uj,uk,sqrtab;
	BaseLine(int m,int n,int minl,float p1,float p2,float p3);
	~BaseLine(void);
	//initial the line parameter

	void initial(int cx,int cy,float pang);
	void updateAnglePrams(float angle);
	void updateLinePrams(int x1,int y1);
	void getLine(Point3i*reg);
	float pt_line_dis(float xx,float yy);
	void BaseLine::getlength();
	void region2line( Point3i *reg, int ed1,int ed2,int ed3,int ed4,
                         double prec);
	void newValidationPixel(int i);

	bool lengthSatisfy();

	float getAnchorThreshold();
	void BaseLine::newValidationPixel2(int i);
	float withinLength(float x, float y);
	void reverseDir(int cx,int cy);
};

