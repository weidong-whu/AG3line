#include"GeneralFuncs.h"
#include"BaseLine.h"

 float  sqrt2(const float x)
{
  const float xhalf = 0.5f*x;
 
  union // get bits for floating value
  {
    float x;
    int i;
  } u;
  u.x = x;
  u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
  return x*u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy

}
double angle_diff(double a, double b)
{
    a -= b;
    while ( a <= -M_PI ) a += M_2__PI;
    while ( a >   M_PI ) a -= M_2__PI;
    if ( a < 0.0 ) a = -a;
    return a;
}

int round(float number)
{
    return (number > 0.0) ? (number + 0.5) : (number - 0.5); 
}


 
void getAllFiles(std::string path, std::vector<std::string>&files, std::string fileType){
 
    //ÎÄ¼þ¾ä±ú
    long hFile = 0;
    struct _finddata_t  fileInfo;
    std::string p;
 
    if ((hFile = _findfirst(p.assign(path).append("\\*" + fileType).c_str(), &fileInfo)) != -1){
        do{
            files.push_back(p.assign(path).append("\\").append(fileInfo.name));
        } while (_findnext(hFile, &fileInfo) == 0);
 
        _findclose(hFile);//¹Ø±Õ¾ä±ú
 
    }
 
}
void getPixels4Dir(float theta,float* xxs,float* yys)
{
		switch (direction(theta))
			{
			case 1:
				xxs[1]=xxs[0]+1;
				yys[1]=yys[0];
				xxs[2]=xxs[0]-1;
				yys[2]=yys[0];
				break;
			case 4:
				xxs[1]=xxs[0]+1;
				yys[1]=yys[0]+1;
				xxs[2]=xxs[0]-1;
				yys[2]=yys[0]-1;
			
				break;
			case 3:
				xxs[1]=xxs[0];
				yys[1]=yys[0]+1;
				xxs[2]=xxs[0];
				yys[2]=yys[0]-1;
				break;
			case 2:
				xxs[1]=xxs[0]-1;
				yys[1]=yys[0]+1;
				xxs[2]=xxs[0]+1;
				yys[2]=yys[0]-1;
				break;
			}


}

int direction(double regang)
{	/*
	1:
	0 0 0
	1 o 1
	0 0 0

	2:
	1 0 0
	0 o 0
	0 0 1

	3:
	0 1 0
	0 o 0
	0 1 0

	4:
	0 0 1
	0 o 0
	1 0 0
	*/
	if ((regang>-0.3927&&regang<0.3927)||
		(regang>2.7489||regang<-2.7489))
		return 1;

	if ((regang>=0.3927&&regang<=1.1781)||
		(regang<=-1.9635&&regang>= -2.7489))
		return 2;

	if ((regang>=1.1781&&regang<= 1.9635)||
		(regang>=-1.9635&&regang<= -1.1781))
		return 3;


	if ((regang>=1.9635&&regang<= 2.7489)||
		(regang<= -0.3927&&regang>=-1.1781))
		return 4;

	return 0;
}


void Bresenham(int x1,
	int y1,
	int const x2,
	int const y2,Point3i*linereg,int *size)
{
	*size=0;
	int delta_x(x2 - x1);
	// if x1 == x2, then it does not matter what we set here
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	// if y1 == y2, then it does not matter what we set here
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;
	linereg[*size].x=x1;
	linereg[*size].y=y1;
	//plot(x1, y1);

	if (delta_x >= delta_y)
	{
		// error may go below zero
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2)
		{
			// reduce error, while taking into account the corner case of error == 0
			if ((error > 0) || (!error && (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
			}
			// else do nothing

			error += delta_y;
			x1 += ix;

			//plot(x1, y1);
			linereg[*size].x=x1;
			linereg[*size].y=y1;
			*size=*size+1;
		}
	}
	else
	{
		// error may go below zero
		int error(delta_x - (delta_y >> 1));

		while (y1 != y2)
		{
			// reduce error, while taking into account the corner case of error == 0
			if ((error > 0) || (!error && (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}
			// else do nothing

			error += delta_x;
			y1 += iy;

			// plot(x1, y1);
			linereg[*size].x=x1;
			linereg[*size].y=y1;
			*size=*size+1;
		}
	}
}


bool goVertical(float regang)
{	//check if the angle goes horizon

	bool ish=(regang>-0.7854&&regang<0.7854)||
		(regang>2.3563&&regang<3.15)||
		(regang<-2.3563&&regang>-3.15);
	return ish;
}


void EnumerateRectPoints1(double sx, double sy, double ex, double ey, Point3i*linereg,int*pNoPoints , double width,int xsize)
{
	double vxTmp[4], vyTmp[4];
	double vx[4], vy[4];
	int n, offset;

	double x1 = sx;
	double y1 = sy;
	double x2 = ex;
	double y2 = ey;
	

	double dx = x2 - x1;
	double dy = y2 - y1;
	double vLen = sqrt2(dx*dx + dy*dy);

	// make unit vector
	dx = dx / vLen;
	dy = dy / vLen;

	/* build list of rectangle corners ordered
	in a circular way around the rectangle */
	vxTmp[0] = x1 - dy * width / 2.0;
	vyTmp[0] = y1 + dx * width / 2.0;
	vxTmp[1] = x2 - dy * width / 2.0;
	vyTmp[1] = y2 + dx * width / 2.0;
	vxTmp[2] = x2 + dy * width / 2.0;
	vyTmp[2] = y2 - dx * width / 2.0;
	vxTmp[3] = x1 + dy * width / 2.0;
	vyTmp[3] = y1 - dx * width / 2.0;

	/* compute rotation of index of corners needed so that the first
	Point3i has the smaller x.

	if one side is vertical, thus two corners have the same smaller x
	value, the one with the largest y value is selected as the first.
	*/
	if (x1 < x2 && y1 <= y2) offset = 0;
	else if (x1 >= x2 && y1 < y2) offset = 1;
	else if (x1 > x2 && y1 >= y2) offset = 2;
	else                          offset = 3;

	/* apply rotation of index. */
	for (n = 0; n<4; n++) {
		vx[n] = vxTmp[(offset + n) % 4];
		vy[n] = vyTmp[(offset + n) % 4];
	} //end-for

	  /* Set a initial condition.

	  The values are set to values that will cause 'ri_inc' (that will
	  be called immediately) to initialize correctly the first 'column'
	  and compute the limits 'ys' and 'ye'.

	  'y' is set to the integer value of vy[0], the starting corner.

	  'ys' and 'ye' are set to very small values, so 'ri_inc' will
	  notice that it needs to start a new 'column'.

	  The smaller integer coordinate inside of the rectangle is
	  'ceil(vx[0])'. The current 'x' value is set to that value minus
	  one, so 'ri_inc' (that will increase x by one) will advance to
	  the first 'column'.
	  */
	int x = (int)ceil(vx[0]) - 1;
	int y = (int)ceil(vy[0]);
	double ys = -DBL_MAX, ye = -DBL_MAX;

	int noPoints = 0;
	while (1) {
		/* if not at end of exploration,
		increase y value for next pixel in the 'column' */
		y++;

		/* if the end of the current 'column' is reached,
		and it is not the end of exploration,
		advance to the next 'column' */
		while (y > ye && x <= vx[2]) {
			/* increase x, next 'column' */
			x++;

			/* if end of exploration, return */
			if (x > vx[2]) break;

			/* update lower y limit (start) for the new 'column'.

			We need to interpolate the y value that corresponds to the
			lower side of the rectangle. The first thing is to decide if
			the corresponding side is

			vx[0],vy[0] to vx[3],vy[3] or
			vx[3],vy[3] to vx[2],vy[2]

			Then, the side is interpolated for the x value of the
			'column'. But, if the side is vertical (as it could happen if
			the rectangle is vertical and we are dealing with the first
			or last 'columns') then we pick the lower value of the side
			by using 'inter_low'.
			*/
			if ((double)x < vx[3]) {
				/* interpolation */
				if (fabs(vx[0] - vx[3]) <= 0.01) {
					if (vy[0]<vy[3]) ys = vy[0];
					else if (vy[0]>vy[3]) ys = vy[3];
					else     ys = vy[0] + (x - vx[0]) * (vy[3] - vy[0]) / (vx[3] - vx[0]);
				}
				else
					ys = vy[0] + (x - vx[0]) * (vy[3] - vy[0]) / (vx[3] - vx[0]);

			}
			else {
				/* interpolation */
				if (fabs(vx[3] - vx[2]) <= 0.01) {
					if (vy[3]<vy[2]) ys = vy[3];
					else if (vy[3]>vy[2]) ys = vy[2];
					else     ys = vy[3] + (x - vx[3]) * (y2 - vy[3]) / (vx[2] - vx[3]);
				}
				else
					ys = vy[3] + (x - vx[3]) * (vy[2] - vy[3]) / (vx[2] - vx[3]);
			} //end-else

			  /* update upper y limit (end) for the new 'column'.

			  We need to interpolate the y value that corresponds to the
			  upper side of the rectangle. The first thing is to decide if
			  the corresponding side is

			  vx[0],vy[0] to vx[1],vy[1] or
			  vx[1],vy[1] to vx[2],vy[2]

			  Then, the side is interpolated for the x value of the
			  'column'. But, if the side is vertical (as it could happen if
			  the rectangle is vertical and we are dealing with the first
			  or last 'columns') then we pick the lower value of the side
			  by using 'inter_low'.
			  */
			if ((double)x < vx[1]) {
				/* interpolation */
				if (fabs(vx[0] - vx[1]) <= 0.01) {
					if (vy[0]<vy[1]) ye = vy[1];
					else if (vy[0]>vy[1]) ye = vy[0];
					else     ye = vy[0] + (x - vx[0]) * (vy[1] - vy[0]) / (vx[1] - vx[0]);
				}
				else
					ye = vy[0] + (x - vx[0]) * (vy[1] - vy[0]) / (vx[1] - vx[0]);

			}
			else {
				/* interpolation */
				if (fabs(vx[1] - vx[2]) <= 0.01) {
					if (vy[1]<vy[2]) ye = vy[2];
					else if (vy[1]>vy[2]) ye = vy[1];
					else     ye = vy[1] + (x - vx[1]) * (vy[2] - vy[1]) / (vx[2] - vx[1]);
				}
				else
					ye = vy[1] + (x - vx[1]) * (vy[2] - vy[1]) / (vx[2] - vx[1]);
			} //end-else

			  /* new y */
			y = (int)ceil(ys);
		} //end-while

		  // Are we done?
		if (x > vx[2]) break;

		linereg[noPoints].x = x;
		linereg[noPoints].y = y;
		linereg[noPoints].z=y*xsize+x;
		noPoints++;
	} //end-while

	*pNoPoints = noPoints;
}

void getRoundPixels(BaseLine*bl,Point3i*reg,int *size)
{
	float dx=(bl->x2-bl->x1)/bl->length;
	float dy=(bl->y2-bl->y1)/bl->length;

	*size=round(bl->length);
	int locat;
	for(int i=0;i<*size;i++)
	{
		
		reg[i].x=round(bl->x1+i*dx);
		reg[i].y=round(bl->y1+i*dy);
		locat=reg[i].y*bl->xsize+reg[i].x;
		reg[i].z=locat;
	}
}

int partition(float *array, int left, int right)
{
	if (array == NULL)
		return -1;

	int pos = right;
	right--;
	while (left <= right)
	{
		while (left < pos && array[left] <= array[pos])
			left++;

		while (right >= 0 && array[right] > array[pos])
			right--;

		if (left >= right)
			break;

		swap(array[left], array[right]);
	}
	swap(array[left], array[pos]);

	return left;
}

float getMidIndex(float *array1, int size)
{
	if (array1 == NULL || size <= 0)
		return -1;

	int left = 0;
	int right = size - 1;
	int midPos = right >> 1;
	int index = -1;

	while (index != midPos)
	{
		index = partition(array1, left, right);

		if (index < midPos)
		{
			left = index + 1;
		}
		else if (index > midPos)
		{
			right = index - 1;
		} 
		else
		{
			break;
		}
	}

	assert(index == midPos);
	return array1[index];
}

//mean value of the arry
float mean(float*value,int ed)
{
	float sum=0;
	for(int i=0;i<ed;i++)
	{
		sum+=value[i];
	}
	return sum/(float)ed;
}

void sdeviation(float*value,int ed,float*m,float*std)
{
	//*getMidIndex(value,ed);
	
	float meanv=mean(value,ed);
	*m=meanv;
	float cell,sum;
	sum=0;
	for(int i=0;i<ed;i++)
	{
		cell=value[i]-(meanv);
		cell=cell*cell;
		sum+=cell;
	}
	ed=ed-1;
	*std=sqrt2(sum/(float)ed);
}


