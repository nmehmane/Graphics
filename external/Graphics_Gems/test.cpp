#include <cstdlib>
#include <cstdio>
#include "NearestPoint.h"
#include "GraphicsGems.h"
/*
 *  main :
 *	Given a cubic Bezier curve (i.e., its control points), and some
 *	arbitrary point in the plane, find the point on the curve
 *	closest to that arbitrary point.
 */
int main()
{
   
 static Point2 bezCurve[4] = {	/*  A cubic Bezier curve	*/
	{ 0.0, 0.0 },
	{ 1.0, 2.0 },
	{ 3.0, 3.0 },
	{ 4.0, 2.0 },
    };
    static Point2 arbPoint = { 3.5, 2.0 }; /*Some arbitrary point*/
    Point2	pointOnCurve;		 /*  Nearest point on the curve */
    double t_param = 0; 
    /*  Find the closest point */
    pointOnCurve = NearestPointOnCurve(arbPoint, bezCurve, &t_param);
    printf("pointOnCurve : (%4.4f, %4.4f)\n", pointOnCurve.x,
		pointOnCurve.y);
    printf("parameter of the point : %lf\n",t_param);
}
