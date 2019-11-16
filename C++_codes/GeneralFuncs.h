/** Label for pixels with undefined gradient. */
#include<math.h>
#include<iostream>
#include <io.h>
#include <string>
#include <vector>
#include <fstream>
#include <opencv2/opencv.hpp>
#include"BaseLine.h"
#include"fastatan2.h"
#define NOTDEF -1024.0
#define PI_8 0.3927
#define PI_16 0.1963495408
#define PI_32 0.0981747704246
#define M_PI   3.1415926535
#define PI_2 1.570796
#define M_2__PI  6.28318530718
/** 3/2 pi */
#define M_3_2_PI 4.71238898038
#define NOTUSED 0
#define USED 1
#define GRADE_2 27.0400
#define SQRT_MAGIC_F 0x5f3759df
double angle_diff(double a, double b);
int round(float number);
float mean(float*value,int ed);
void sdeviation(float*value,int ed,float *std,float *m);
void getAllFiles(std::string path, std::vector<std::string>&files, std::string fileType);
void getPixels4Dir(float theta,float* xxs,float* yys);
int direction(double regang);
void Bresenham(int x1,
	int y1,
	int const x2,
	int const y2,cv::Point3i*linereg,int *size);

bool goVertical(float regang);
void EnumerateRectPoints1(double sx, double sy, double ex, double ey, cv::Point3i*linereg, int * pNoPoints,double width,int xsize);
void getRoundPixels(BaseLine*bl,Point3i*reg,int *size);

 float  sqrt2(const float x);