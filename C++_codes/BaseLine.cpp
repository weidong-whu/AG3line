#include "BaseLine.h"
#include "GeneralFuncs.h"

BaseLine::BaseLine(int m,int n,int minl,float p1,float p2,float p3)
{
	minlength=minl;
	xsize=n;
	ysize=m;
	pdense[0]=p1;
	pdense[1]=p2;
	pdense[2]=p3;
}

BaseLine::~BaseLine(void)
{
}

void BaseLine::initial(int cx,int cy,float pang)
{

	counter=1;
	mx=cx;
	my=cy;
	segang=pang;
	sumdx=std::cos(pang);
	sumdy=std::sin(pang);
	c1=0;
	c2=0;
	c3=0;
	c4=0;

	dx=-sin(pang);
	dy=cos(pang);

	theta=atan(dy/dx);

	sx=cx;
	sy=cy;

	ui = (my+10*dy)  - my ;
    uj = mx  - (mx+10*dx) ;
    uk = (mx+10*dx) * my - mx * (my+10*dy);
	sqrtab=sqrt2( ui*ui+uj*uj);
	//cout<<sqrt( ui*ui+uj*uj)-sqrt2( ui*ui+uj*uj)<<endl;
}
void BaseLine::updateAnglePrams(float angle)
{
	this->sumdx+=std::cos(angle);
	this->sumdy+=std::sin(angle);
	this->segang=atan2approx(sumdy,sumdx);
}
void BaseLine::updateLinePrams(int x1,int y1)
{
	counter=counter+1;
	a=1/(float)counter;
	x1_mx=(x1-mx);
	y1_my=(y1-my);
	a_1=(1-a);

	mx=a_1*mx+a*x1;
	my=a_1*my+a*y1;

	c1=(x1_mx*x1_mx*a+c1)*a_1;
	c2=(x1_mx*y1_my*a+c2)*a_1;
	c3=(x1_mx*y1_my*a+c3)*a_1;
	c4=(y1_my*y1_my*a+c4)*a_1;
	
    lambda=0.5*(c1+c4-sqrt2((c1-c4)*(c1-c4)+4*c2*c2));
	if (c1>c4)
		theta=atan2approx(lambda-c1,-c2);
	else
		theta=atan2approx(-c2,lambda-c4);

	//keep the direction
    //if ( angle_diff(theta,segang) > PI_8 ) theta += M_PI;
	dx=std::sin(theta);
    dy=std::cos(theta);

	if(abs(dx)>abs(dy))
	{
		if((x1_mx)*(dx)<0)
		{
			dx=-dx;
			dy=-dy;
		}
		t=(x1-mx)/dx;
		sx=x1;
		sy=my+t*dy;

	}
	else
	{
		if((y1_my)*(dy)<0)
		{
			dx=-dx;
			dy=-dy;
		}
		t=(y1-my)/dy;
		sy=y1;
		sx=mx+t*dx;
	}
	ui = sy  - my ;
    uj = mx  - sx ;
    uk = sx * my - mx * sy;
	sqrtab=sqrt2( ui*ui+uj*uj);
	sx=round(sx);
	sy=round(sy);
	// cout<<sqrt( ui*ui+uj*uj)-sqrt2( ui*ui+uj*uj)<<endl;
}

void BaseLine::getLine(Point3i*reg)
{
	

    l_min = l_max  = 0.0;
    for (int i = 0; i < counter; i++)
    {
        l =  ( (double) reg[i].x - mx) * dx + ( (double) reg[i].y - my) * dy;
    
        if ( l > l_max ) l_max = l;
        if ( l < l_min ) l_min = l;
       
    }

	  /* store values */
    x1 = mx + l_min * dx;
    y1 = my + l_min * dy;
    x2 = mx + l_max * dx;
    y2 = my + l_max * dy;
	length=sqrt2((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

}
float BaseLine::pt_line_dis(float xx,float yy)
	{
		pt.at<float>(0,0)=xx;
		pt.at<float>(1,0)=yy;

		v1.at<float>(0,0)=mx;
		v1.at<float>(1,0)=my;

	    v2.at<float>(0,0)=sx;
		v2.at<float>(1,0)=sy;
		/* a = v1 - v2;
		   b = pt - v2;
           d = norm(cross(a,b)) / norm(a);
	    */
		return cv::norm((v1-v2).cross(pt-v2))/cv::norm(v1-v2);
	}
void BaseLine::getlength()
{
	length=sqrt2((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	
}

double BaseLine::get_theta( Point3i *reg, int ed1,int ed2,int ed3,int ed4, double x, double y,
                         double reg_angle, double prec )
{
    double lambda, theta, weight;
    double Ixx = 0.0;
    double Iyy = 0.0;
    double Ixy = 0.0;
    int i;

   double x_x,y_y;
    /* compute inertia matrix */
    for (i = ed1; i < ed2; i++)
    {
        x_x=reg[i].x - x;
		y_y=reg[i].y - y;
        Ixx += y_y * y_y;
        Iyy +=x_x * x_x;
        Ixy -= x_x * y_y;
    }
	    /* compute inertia matrix */
    for (i = ed3; i < ed4; i++)
    {
        
        x_x=reg[i].x - x;
		y_y=reg[i].y - y;
        Ixx += y_y * y_y;
        Iyy +=x_x * x_x;
        Ixy -= x_x * y_y;
    }

    /* compute smallest eigenvalue */
    lambda = 0.5 * ( Ixx + Iyy - sqrt2( (Ixx - Iyy) * (Ixx - Iyy) + 4.0 * Ixy * Ixy ) );

    /* compute angle */
    theta = fabs(Ixx) > fabs(Iyy) ? atan2approx(lambda - Ixx, Ixy) :atan2approx(Ixy, lambda - Iyy);

    /* The previous procedure don't cares about orientation,
       so it could be wrong by 180 degrees. Here is corrected if necessary. */
    if ( angle_diff(theta, reg_angle) > prec ) theta += M_PI;

    return theta;
}
void BaseLine::region2line( Point3i *reg, int ed1,int ed2,int ed3,int ed4,
                         double prec)
{
	double x, y, dx, dy, l, w, theta, weight, sum, l_min, l_max, w_min, w_max;
	int i,reg_size;


	/* center of the region:

	It is computed as the weighted sum of the coordinates
	of all the pixels in the region. The norm of the gradient
	is used as the weight of a pixel. The sum is as follows:
	cx = \sum_i G(i).x_i
	cy = \sum_i G(i).y_i
	where G(i) is the norm of the gradient of pixel i
	and x_i,y_i are its coordinates.
	*/
	x = y = sum =reg_size= 0.0;

	for (i = ed1; i <ed2; i++)
	{

		x += (double) reg[i].x ;
		y += (double) reg[i].y ;
		reg_size++;

	}
	for (i = ed3; i <ed4; i++)
	{

		x += (double) reg[i].x ;
		y += (double) reg[i].y ;
		reg_size++;
	}
	x /= reg_size;
	y /= reg_size;

	/* theta */
	theta = get_theta(reg, ed1,ed2,ed3,ed4, x, y, segang, prec);

	/* length and width:

	'l' and 'w' are computed as the distance from the center of the
	region to pixel i, projected along the rectangle axis (dx,dy) and
	to the orthogonal axis (-dy,dx), respectively.

	The length of the rectangle goes from l_min to l_max, where l_min
	and l_max are the minimum and maximum values of l in the region.
	Analogously, the width is selected from w_min to w_max, where
	w_min and w_max are the minimum and maximum of w for the pixels
	in the region.
	*/
	dx = cos(theta);
	dy = sin(theta);
	l_min = l_max = w_min = w_max = 0.0;
	double x_x,y_y;
	for (i = ed1; i <ed2; i++)
	{
		x_x=reg[i].x - x;
		y_y=reg[i].y - y;
		l =  x_x * dx + y_y * dy;
		w = - x_x * dy +  y_y * dx;

		if ( l > l_max ) l_max = l;
		if ( l < l_min ) l_min = l;
		if ( w > w_max ) w_max = w;
		if ( w < w_min ) w_min = w;
	}


	for (i = ed3; i <ed4; i++)
	{
		x_x=reg[i].x - x;
		y_y=reg[i].y - y;
		l =  x_x * dx + y_y * dy;
		w = - x_x * dy +  y_y * dx;

		if ( l > l_max ) l_max = l;
		if ( l < l_min ) l_min = l;
		if ( w > w_max ) w_max = w;
		if ( w < w_min ) w_min = w;
	}

	/* store values */
	x1 = x + l_min * dx;
	y1 = y + l_min * dy;
	x2 = x + l_max * dx;
	y2 = y + l_max * dy;
	theta = theta;
	dx = dx;
	dy = dy;
	mx=x;
	my=y;


	ui = y1 * 1 - my * 1;
    uj = mx * 1 - x1 * 1;
    uk = x1 * my - mx * y1;
	sqrtab=sqrt2( ui*ui+uj*uj);
}

void BaseLine::newValidationPixel(int i)
{
	xxs[0]=round((int)sx+(float)i*dx);
	yys[0]=round((int)sy+(float)i*dy);
	//cout<<direction(bl->theta)<<endl;
	if (goVertical(theta))
	{
		xxs[1]=xxs[0]+1;
		yys[1]=yys[0];
		xxs[2]=xxs[0]-1;
		yys[2]=yys[0];

	}
	else
	{

		xxs[1]=xxs[0];
		yys[1]=yys[0]+1;
		xxs[2]=xxs[0];
		yys[2]=yys[0]-1;

	}
}

void BaseLine::newValidationPixel2(int i)
{
	xxs[0]=round((int)sx+(float)i*dx);
	yys[0]=round((int)sy+(float)i*dy);
	getPixels4Dir(theta,xxs,yys);
	
}

bool BaseLine::lengthSatisfy()
{
		return length>minlength;
}
float  BaseLine::getAnchorThreshold()
{
	return (float)(pdense[0] *std::pow(length,pdense[1] )+pdense[2]) ;
}

float BaseLine::withinLength(float x, float y)
{
	return std::abs(ui*x+uj*y+uk)/sqrtab;
}
void  BaseLine::reverseDir(int cx,int cy)
{
	dx=-dx;dy=-dy;
	sx=cx,sy=cy;
}