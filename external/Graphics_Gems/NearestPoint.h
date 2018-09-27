#ifndef N_POINT
#define N_POINT

#include "GraphicsGems.h"

Point2  NearestPointOnCurve(Point2 P, Point2 *V, double *tt);
static	int	FindRoots(Point2 *w, int degree, double *t, int depth);
static	Point2	*ConvertToBezierForm(Point2 P,Point2 *V);
static double ComputeXIntercept(Point2 *V, int degree);
static int ControlPolygonFlatEnough( const Point2* V, int degree );
static int CrossingCount(Point2 *V, int degree);
static Point2 Bezier(Point2 *V, int degree, double t, Point2 *Left, Point2 *Right);
static Vector2 V2ScaleII(Vector2 *v, double s);

#endif
