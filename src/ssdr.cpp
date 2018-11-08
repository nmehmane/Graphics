#include "ssdr.h"
#include <assert.h>
#include <cassert>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <gsl/gsl_linalg.h> // for svd
#include <Eigen/Dense>
#include <iostream>
#include "gco/GCoptimization.h"
#include "Graphics_Gems/GraphicsGems.h"
#include "Graphics_Gems/NearestPoint.h"

namespace ssdr_ns
{
ssdr::ssdr()
{
}

ssdr::~ssdr()
{
}


/*******************************************************************************************************/
void ssdr::perform_ssdr()
{
    rp_SamplePoints = sample_curves( rp_CurveEndPoints, 
                                     rp_CurveMiddlePoints,
                                     rp_Curves, 
                                     1, 
                                     &rp_SamplePoint_tg,
                                     &rp_SamplePoint_Curves );

    std::cout << " Points \n" << std::endl;
    for( auto& p : rp_SamplePoints )
    {
        std::cout << p << std::endl;
    }

    std::cout << " Tangents \n" << std::endl;
    for( auto& q : rp_SamplePoint_tg )
    {
        std::cout << q << std::endl;
    }
    printf("size = %lu\n",rp_SamplePoint_tg.size());
    rp_SamplePoint_tg =  normalize_tangents( rp_SamplePoint_tg );

    for( auto & ntg : rp_SamplePoint_tg )
    {
        std::cout << ntg << std::endl;
    }

    dfp_SamplePoints = sample_curves( dfp_CurveEndPoints, 
                                      dfp_CurveMiddlePoints,
                                      dfp_Curves, 
                                      1, 
                                      &dfp_SamplePoint_tg,
                                      &dfp_SamplePoint_Curves );

    std::cout << " Points \n" << std::endl;
    for( auto& p : dfp_SamplePoints )
    {
        std::cout << p << std::endl;
    }

    std::cout << " Tangents \n" << std::endl;
    for( auto& q : dfp_SamplePoint_tg )
    {
        std::cout << q << std::endl;
    }
    printf(" size = %lu \n", dfp_SamplePoint_tg.size() );
    dfp_SamplePoint_tg = normalize_tangents( dfp_SamplePoint_tg );

    for( auto & ntg : dfp_SamplePoint_tg )
    {
        std::cout << ntg << std::endl;
    }
    // find out what curves are neighbors in the rest pose
    for( auto & pp : rp_Curves )
    {
        std::cout << " curve \n" << pp(0) << "  " << pp(1) << "  " << pp(2) << "  " << pp(3) << std::endl;
    }
    Eigen::MatrixXd neighbor_curves = set_curve_neighborhood( rp_Curves );
    
    // will hold the n x n Transformations
    std::vector<Eigen::Matrix2d> Rs;
    std::vector<Eigen::Vector2d> Ts;
    
    // we compute n x n candidate R,T by pairing each point of the rest pose with 
    // each point in the deformed pose
    initialize_transformations( rp_SamplePoints, dfp_SamplePoints, 
                                rp_SamplePoint_tg, dfp_SamplePoint_tg,
                                &Rs, &Ts );
    /*std::cout << " Rs \n" << std::endl;

    for( auto & r : Rs )
    {
        std::cout <<  r << std::endl;
    }

    std::cout << " Ts \n" << std::endl;
    
    for( auto & t : Ts )
    {
        std::cout <<  t << std::endl;
    }
    */
    int iter = 0;
    // now we have n x n Labels we copute the optimal label assignment.
    do{
        // to compute the the optimal label assignment we will compute the data term,
        // the smoothness term 
        int* data = set_data_term( rp_SamplePoints, dfp_SamplePoints, Rs, Ts );
        //for( int i = 0 ; i < rp_SamplePoints.size() * Rs.size() ; i++ )
        //{
        //    printf(" data term %d = %d \n",i ,data[i]);
        //}
        int *smooth = set_smoothness_term( Rs.size() );
       // for( int i = 0 ; i < Rs.size()*Rs.size() ; i++ )
       // {
       //     printf(" smooth term %d = %d \n",i ,smooth[i]);
       // }
        // use the returned matrix to calculate the neighborhood information needed
        // for optainnig the optimal label assignment
        int *numNeighbors;
        int **neighborsIndexes;
        int** neighborsWeights;
   /*     printf(" curve numbers \n");
        for( int ii : rp_SamplePoint_Curves)
        {
            std::cout << ii << std::endl ;
        }
        printf(" hi\n");
     */   compute_neighborhood_info( rp_SamplePoint_Curves, neighbor_curves, 
                                   &numNeighbors, 
                                   &neighborsIndexes, 
                                   &neighborsWeights );
/*        for( int i = 0; i < rp_SamplePoints.size() ; i++ )
        {
            printf(" \nnum neighbors %d = %d \n",i, numNeighbors[i] );
            for( int j = 0 ; j < numNeighbors[i] ; j++ )
            {
                printf(" %d  ",neighborsIndexes[i][j]);
            }
        }
        for( int i = 0; i < rp_SamplePoints.size() ; i++ )
        {
            printf(" \nnum neighbors %d = %d \n",i, numNeighbors[i] );
            for( int j = 0 ; j < numNeighbors[i] ; j++ )
            {
                printf(" %d  ",neighborsWeights[i][j]);
            }
        }
  */    
        int num_points = rp_SamplePoints.size();
        int num_labels = Rs.size();
        std::vector<int> label_assignment = assign_labels( num_points, num_labels, data, smooth,
                                                 numNeighbors, neighborsIndexes, neighborsWeights );
        printf("label assignment\n");
        for( int l : label_assignment )
        {
            printf(" %d   ",l);
        }
        // now we have a label assignment we separate the rest pose clusters based on this
        std::vector<std::vector<int>> pose_1_clusters;
        std::vector<int> label_numbers;

        set_restpose_clusters( label_assignment, num_points, 
                               &pose_1_clusters, &label_numbers); 
        int uu = 0;
        for( auto& rr : pose_1_clusters )
        {
            printf("\ncluster label %d\n",label_numbers.at(uu));
            
            for( auto& oo : rr )
            {
                printf(" %d  ", oo );
            }
            uu++;
        }


        std::vector<std::vector<int>> pose_2_clusters;
        // given the restpose clusters we compute the corresponding deformed pose clusters
        compute_corresponding_clusters( rp_SamplePoints, dfp_SamplePoints, label_numbers, Rs, Ts,
                                        pose_1_clusters, &pose_2_clusters );
        for( auto& rr : pose_2_clusters )
        {
            printf("\ncluster 2\n");
            for( auto& oo : rr )
            {
                printf(" %d  ", oo );
            }
        }
        // now having clusters in the rest pose and the corresponding clusters in the deformed pose
        // we recompute Transforms
        Rs.clear();
        Ts.clear();

        recompute_RT(  rp_SamplePoints, dfp_SamplePoints, &Rs, &Ts, pose_1_clusters, pose_2_clusters );
        std::cout << " Rs \n" << std::endl;

        for( auto & r : Rs )
        {
            std::cout <<  r << std::endl;
        }

        std::cout << " Ts \n" << std::endl;
        
        for( auto & t : Ts )
        {
            std::cout <<  t << std::endl;
        }
    
        
    }while( (++iter) < 5 );

                                             
}
/*******************************************************************************************************/
// used to sample the inputed vector of curves at samples_per_curve rate 
// Returns the vector of sample points and calculates the tangents at the sampled points and 
// returns them in the 5th argument
// ++ modify to output the curve number of each point
std::vector<Point> ssdr::sample_curves( const std::vector<Point>& end_points, 
                                        const std::vector<Point>& middle_points, 
                                        const std::vector<Eigen::Vector4i>& cubic_curves,
                                        const int samples_per_curve,
                                        std::vector<Eigen::Vector2d> *p_samples_tg,
                                        std::vector<int>* p_curves)
{
    assert( samples_per_curve > 0 );

    //will hold the pose's sampled points
    std::vector<Point> p_samples;

    
    const double delta_t = 1.0 / samples_per_curve;
    
    double t = 0;
    Point p0, p1, p2, p3;
    
    int j = 0;

    // each curve is sampled samples_per_curve times
    for( auto & curve : cubic_curves )
    {
        t = 0;
        //end points
        p0 = end_points.at(curve[0]);
        p3 = end_points.at(curve[3]);
        
        //tangent control points
        p1 = middle_points.at(curve[1]);
        p2 = middle_points.at(curve[2]);
        
        for( int i = 0 ; i < samples_per_curve ; i++ )
        {
            (*p_curves).push_back(j);
            p_samples.push_back( evaluate_curve( p0, p1, p2, p3, t ) );
            (*p_samples_tg).push_back( evaluate_curve_tg( p0, p1, p2, p3, t ) );
            
            t += delta_t;
        }
        j++;
        
    }
    return p_samples;
}

/*******************************************************************************************************/
// given the 4 control points evaluates the curve at parameter t
Point ssdr::evaluate_curve( const Point& p0, const Point& p1,
                            const Point& p2, const Point& p3,
                            const double t )
{
    /*
              [ 1  t  t^2  t^3 ] * [  1  0  0  0 ] * [  P0 ]
                                   [ -3  3  0  0 ]   [  P1 ]
                                   [  3 -6  3  0 ]   [  P2 ]
                                   [ -1  3 -3  1 ]   [  P3 ]
    */
    Eigen::Matrix4d B;
    B <<  1, 0, 0, 0,
         -3, 3, 0, 0,
          3,-6, 3, 0,
         -1, 3,-3, 1;

    Eigen::MatrixXd U(1,4);
    U << 1, t, t * t, pow(t,3);

    Eigen::MatrixXd T = U * B;

    return ( T(0,0) * p0 ) + ( T(0,1) * p1 ) + ( T(0,2) * p2 ) + ( T(0,3) * p3 ) ;
}   

/*******************************************************************************************************/
// given the 4 control points evaluates the tangent curve at parameter t
Point ssdr::evaluate_curve_tg( const Point& p0, const Point& p1, 
                               const Point& p2, const Point& p3, 
                               const double t )
{
    Point tg =  3 * ( 1 - t ) * ( 1 - t ) * ( p1 - p0 );
    tg += 6 * ( 1 - t ) * t * ( p2 - p1 );
    tg += 3 * t * t * ( p3 - p2 );
    
    return tg;
}   

/*******************************************************************************************************/
// given a vector of tangent points it returns a vector of normalized tangent points
std::vector<Point>  ssdr::normalize_tangents( const std::vector<Point>& sp_tg )
{
    printf("size = %lu\n",sp_tg.size());
    std::vector<Point> sp_normalized_tg;
    for( const auto& td_point : sp_tg )
    {
        std::cout << td_point << std::endl;

        sp_normalized_tg.push_back( td_point.normalized() );
    }
    
    return sp_normalized_tg;
}

/*******************************************************************************************************/
// given a set of points it returns a vector of the points t parameters and e
// a vector of curve numbers to that sppecifies which curve the point belongs to
std::vector<double> ssdr::calculate_param(const std::vector<Point>& points,
                                          const std::vector<Point>& tg_points,
                                          const std::vector<Eigen::Vector4i>& cubic_curves,
                                          std::vector<int>* curve_numbers)
{   
    
    // for each sample point in the rest pose, find the nearest point on each curve and see if the
    // distance is smaller than some epsilon. if if is consider that the point is on the curve and so 
    // record its curve number and its t . knowing t now we can find out what the tangent is
    
    std::vector<double> points_param_t;
    
    for(int i = 0; i < points.size(); i++)
    {
        // get each point 
        Point2 arbPoint = { (points.at(i))(0), (points.at(i))(1) };
        
        int p0_num, p1_num, p2_num, p3_num;

        for(int j = 0; j < cubic_curves.size(); j++)
        {
            p0_num = (cubic_curves.at(j))(0);
            p3_num = (cubic_curves.at(j))(3);
            p1_num = (cubic_curves.at(j))(1);
            p2_num = (cubic_curves.at(j))(2);

            //  get each cubic Bezier curve to test if the point belongs to this curve
            Point2 bezCurve[4] = {	
                { (points.at(p0_num))(0), (points.at(p0_num))(1) },
                { (tg_points.at(p1_num))(0), (tg_points.at(p1_num))(1) },
                { (tg_points.at(p2_num))(0), (tg_points.at(p2_num))(1) },
                { (points.at(p3_num))(0), (points.at(p3_num))(1) },
                };
            
            // will hold the nearest point on the curve
            Point2	pointOnCurve;
            // will hold the parameter of the closest curve point to this point 
            double t_param = 0; 

            //Find the closest point 
            pointOnCurve = NearestPointOnCurve(arbPoint, bezCurve, &t_param);

            // find the distance of the point to this curve point if the distance
            // is small enough consider the point to be on the curve 
            double dista = ( pointOnCurve.x - arbPoint.x ) * ( pointOnCurve.x - arbPoint.x );
            dista += ( ( pointOnCurve.y - arbPoint.y ) * ( pointOnCurve.y - arbPoint.y ) );
            if( dista <= 1 )
            {
                printf("pointOnCurve : (%4.4f, %4.4f)\n", pointOnCurve.x,
                                                      pointOnCurve.y);
                printf("parameter of the point : %lf\n",t_param);
                // record its curve number and its parameter t
                points_param_t.push_back(t_param);
                curve_numbers->push_back(j);
                break;
            }
            
        }
    }
    return points_param_t;
}

/*******************************************************************************************************/
// given a vector of parameters and their respective curve numbers it returns a vector of tangents at t
 std::vector<Point>  ssdr::compute_tangents( const std::vector<double>& points_param_t,
                                             const std::vector<int>& curve_numbers,
                                             const std::vector<Eigen::Vector4i>& cubic_curves,
                                             const std::vector<Point>& end_points,
                                             const std::vector<Point>& tangent_points)
 {
     Point p0, p1, p2, p3;
     int index = 0;
     std::vector<Point> p_tgs;

     for( auto & curve : cubic_curves )
     {
        //end points
        p0 = end_points.at(curve[0]);
        p3 = end_points.at(curve[3]);
        
        //tangent control points
        p1 = tangent_points.at(curve[1]);
        p2 = tangent_points.at(curve[2]);

        p_tgs.push_back( evaluate_curve_tg( p0, p1, p2, p3, points_param_t.at(index) ) );
     }
     return p_tgs;
 }

/*******************************************************************************************************/

// given a cluster pose_1 and a corresponding cluster pose_2 calculate the Rotation and the
// Translation between the two clusters and return them in the 3rd and 4th argument

void ssdr::compute_cluster_RT( std::vector<Point> pose_1, 
                               std::vector<Point> pose_2, 
                               Eigen::Matrix2d *R, 
                               Eigen::Vector2d *T )
{
    assert( pose_1.size() == pose_2.size() );
    
    int n = pose_1.size();

    // hold the centers of each cluster
    Point p_bar(0, 0);
    Point q_bar(0, 0);
    
    // calculate centers
    for( int i = 0; i < n ; i++)
    {
        p_bar += pose_1.at(i);
        q_bar += pose_2.at(i);
    } 
    
    p_bar /= n;
    q_bar /= n;
    
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2,n);
    Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(2,n);

    //set the columns of the Matrices
    for( int i = 0; i < n ; i++)
    {
        X.col(i) = pose_1.at(i) - p_bar;
        Y.col(i) = pose_2.at(i) - q_bar;
    }
    
  //  std::cout << " X" << std::endl << X << std::endl;
  //  std::cout << " Y" << std::endl << Y << std::endl;

    Eigen::MatrixXd S_tmp = X * Y.transpose();
    
 //   std::cout << " covariance matrix" << std::endl << S_tmp << std::endl;
    
    // Compute SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.compute( S_tmp, Eigen::ComputeThinU | Eigen::ComputeThinV);
 //   std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
   // std::cout << "Its left singular vectors are the columns of the thin U matrix:"
     //         << std::endl << svd.matrixU() << std::endl;
   // std::cout << "Its right singular vectors are the columns of the thin V matrix:"
    //          << std::endl << svd.matrixV() << std::endl;

    // The optimal rotation is V diag( 1,1,det(V*U.traspose) ) U.transpose

    Eigen::MatrixXd Ut = svd.matrixU().transpose();
    Eigen::MatrixXd diag_mat = Eigen::MatrixXd::Identity(2,2);
    diag_mat(1,1) = (svd.matrixV() * Ut).determinant();
    *R = (svd.matrixV() * diag_mat ) * Ut;
    *T = q_bar - ((*R)*p_bar);

   // std::cout << " R" << std::endl << R << std::endl;
   // std::cout << " T" << std::endl << T << std::endl;

}

/*******************************************************************************************************/

void ssdr::initialize_transformations(const std::vector<Point>& pose1_samples,
                                      const std::vector<Point>& pose2_samples,
                                      const std::vector<Point>& pose1_samples_tg,
                                      const std::vector<Point>& pose2_samples_tg,
                                      std::vector<Eigen::Matrix2d>* Rs,
                                      std::vector<Eigen::Vector2d>* Ts )
{
    int i = 0;
    int j = 0;

    // for each sample in the rest-pose we pair it up with each sample from the
    // deformed pose
    for( auto& p : pose1_samples )
    {
        //holds the restpose cluster
        std::vector<Point> cluster1;
        cluster1.push_back( p );
        cluster1.push_back( p + pose1_samples_tg.at(i) );
        i++; 
        j = 0;

        for( auto& q : pose2_samples )
        {
            //holds the deformed pose cluster
            std::vector<Point> cluster2;
            cluster2.push_back( q );
            cluster2.push_back( q + pose2_samples_tg.at(j) );
            j++; 
            Eigen::Matrix2d R;
            Eigen::Vector2d T;
            compute_cluster_RT( cluster1, cluster2, &R, &T );
            (*Rs).push_back( R );
            (*Ts).push_back( T );
        }
    }
}

/*******************************************************************************************************/

// given a vector of points and and point p find the closest point to p among all the 
// points in the vector and return the point's index. if the third argument is not null the 
// minimum distance is also returned
int ssdr::find_closest_point( const Point& p , const std::vector<Point>& points , double* dist )
{
    double min_dist = (points.at(0) - p ).norm( );
    int min_index = 0;

    double curr_dist;
    int curr_index = 0;

    for( auto& curr_p : points )
    {
        curr_dist = ( curr_p - p ).norm( );
       // std::cout << " curr_p \n" << curr_p << "\np\n" << p << " \ncurr dist\n " << curr_dist << "\nmin dist \n" << min_dist << std::endl;
  
        if( curr_dist <= min_dist )
        {
            min_dist = curr_dist;
            min_index = curr_index;
        }

        curr_index++;
    }
    if( dist != nullptr )
    {
        *dist = min_dist;
    }
    return min_index;
}

/*******************************************************************************************************/

int* ssdr::set_data_term( const std::vector<Point>& pose_1,
                          const std::vector<Point>& pose_2,
                          const std::vector<Eigen::Matrix2d>& Rs,
                          const std::vector<Eigen::Vector2d>& Ts )
{
    //stores the number of labels
    int num_labels = Rs.size();

    //number of sample points
    int num_samples = pose_1.size();

    // setting up the array for data costs i.e. the cost of using Label l_i ( Ri,Ti ) for point p
    int *data = new int[num_samples * num_labels];
    
    // For each vertex A in the rest position 
    for( int i = 0; i < num_samples; i++)
    {
        // For every R&T
        for( int label_index = 0; label_index < num_labels; label_index++ )
        {
            // Apply The current R and T to vertex A and get a position p.
            Point p = Rs.at(label_index) * pose_1.at(i) + Ts.at(label_index);
            // Find the closest vertex in the deformed pose to this point
            double temp;
            find_closest_point( p , pose_2 , &temp);
            data[ i * num_labels + label_index ] = static_cast<int>(temp);

        }
    }
    
    return data;
}

/*******************************************************************************************************/

int* ssdr::set_smoothness_term( int num_labels )
{
    // setting up the array for smooth costs
    int *smooth = new int[ num_labels * num_labels ];
    for ( int l1 = 0; l1 < num_labels; l1++ )
    {
        //printf("smoothness term for %d\n", l1);
        for (int l2 = 0; l2 < num_labels; l2++ )
        {
            smooth[l1+l2*num_labels] = (l1-l2)*(l1-l2) < 1  ? 0 : 10000;
            //printf(" %d ",smooth[l1+l2*num_labels]);
        }
        //printf("\n");
    }
    return smooth;
}

/*******************************************************************************************************/

// currently two curves are considered neighbors if they share an endpoint
// ++ needs to be extended
// sets up an n x n matrix that has M( i ,j ) = 1 if the curve i and j are neighbors

Eigen::MatrixXd  ssdr::set_curve_neighborhood( const std::vector<Eigen::Vector4i>& pose_curves )
{
    int num_curves = pose_curves.size();

    Eigen::MatrixXd neighbor_curves( num_curves, num_curves );

    int i1 = 0;
    int i2 = 0;

    for( auto& c1 : pose_curves )
    {
        i2 = 0;
        for( auto& c2 : pose_curves )
        {
            if( &c1 == &c2 || c1(0) == c2(3) || c1(3) == c2(0) )
            {
                neighbor_curves(i1,i2) = 1;
            }
            else
            {
                neighbor_curves(i1,i2) = 0;
            }
            i2++;
        }

        i1++;
    }
    std::cout << " neighbor_curves\n" << neighbor_curves << std::endl;
    return neighbor_curves;
}

/*******************************************************************************************************/

// num_neighbors = for each point it will have the neighbor count
// neighbors_Indexes = for each point it will hold the index of the point's neighbors
// neighbors_weights = ??

void ssdr::compute_neighborhood_info( const std::vector<int>& curve_numbers,
                                      const Eigen::MatrixXd& curve_nbh,
                                      int **numNeighbors, 
                                      int ***neighborsIndexes, 
                                      int*** neighborsWeights )
{
    int num_points = curve_numbers.size();
    std::vector<std::vector<int>> nbh( num_points );
    *numNeighbors = new int[ num_points ];
    *neighborsIndexes = new int*[ num_points ];
    *neighborsWeights = new int*[ num_points ];
    
    for( int i = 0; i < num_points; i++ )
    {
        (*numNeighbors)[i] = 0;
    }
    
    for( int i = 0; i < num_points; i++ )
    {
        int cn1 = curve_numbers.at(i);

        for( int j = i+1; j < curve_numbers.size(); j++ )
        {
            int cn2 = curve_numbers.at(j);

            if(  curve_nbh( cn1, cn2 ) )
            {
                (nbh.at(i)).push_back(j);
                (nbh.at(j)).push_back(i);

                (*numNeighbors)[i]++;
                (*numNeighbors)[j]++;
            }
        }
    }

    for( int i = 0 ; i < num_points; i++ )
    {
        (*neighborsIndexes)[i] = new int[ (nbh.at(i)).size() ];
        (*neighborsWeights)[i] = new int[ (nbh.at(i)).size() ];

        for( int j = 0 ; j < (nbh.at(i)).size() ; j++ )
        {
            (*neighborsIndexes)[i][j] = (nbh.at(i)).at(j);
            (*neighborsWeights)[i][j] = 1;
        }
    }
}

/*******************************************************************************************************/

std::vector<int> ssdr::assign_labels( int num_points,
                                      int num_labels,
                                      int *data_terms,
                                      int *smoothness_term,
                                      int *numNeighbors, 
                                      int **neighborsIndexes, 
                                      int** neighborsWeights )
{
    std::vector<int> result;
    
    try{

		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph( num_points , num_labels );
		gc->setDataCost( data_terms );
		gc->setSmoothCost( smoothness_term );

        //pass in all neighbor information at once
        gc->setAllNeighbors( numNeighbors , neighborsIndexes , neighborsWeights );

        gc->compute_energy();
		gc->expansion(2);
        gc->swap(2);
        gc->compute_energy();

		for ( int  i = 0; i < num_points; i++ )
        {
			result.push_back( gc->whatLabel(i) );
        }

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}
    return result;

}

/*******************************************************************************************************/

void ssdr::set_restpose_clusters( const std::vector<int>& label_assignment,
                                  const int num_points,
                                  std::vector<std::vector<int>>* pose_1_clusters,
                                  std::vector<int>* label_numbers)
{
    // eliminate duplicate labels to count the number of clusters
    std::vector<int> cp_la = label_assignment;
    std::sort( cp_la.begin(), cp_la.end() );
    cp_la.erase( unique( cp_la.begin(), cp_la.end() ), cp_la.end() );
    int num_clusters = cp_la.size();
    
    int j = 0;
    for( int ln : cp_la )
    {
        std::vector<int> temp;
        for( int i = 0; i < num_points; i++ )
        {
            if( label_assignment.at(i) == ln )
            {
                temp.push_back( i );
            }

        }
        (*pose_1_clusters).push_back( temp );
        (*label_numbers).push_back( ln );
        j++;
    }

}


/*******************************************************************************************************/

// uses the inputed rest pose clusters, applies the trasformation to each point in the
// cluster and finds the closest point in the deformed pose to it. doing this
// for every point gives us the deformed pose clusters
void ssdr::compute_corresponding_clusters( const std::vector<Point>& pose_1,
                                           const std::vector<Point>& pose_2, 
                                           const std::vector<int>& label_assignment,
                                           const std::vector<Eigen::Matrix2d>& Rs,
                                           const std::vector<Eigen::Vector2d>& Ts,
                                           std::vector<std::vector<int>>& pose_1_clusters,
                                           std::vector<std::vector<int>>* pose_2_clusters )
{

    int i = 0;
    // for each label
    for( int label : label_assignment )
    {
        //get the rotation and tanslation related to that label
        Eigen::Matrix2d R = Rs.at( label );
        Eigen::Vector2d T = Ts.at( label );
        
        std::vector<int> temp;
        
        for( int point_index : pose_1_clusters.at( i ) )
        {
            Point p = R * pose_1.at( point_index ) + T;
            // Find the closest vertex in the deformed pose to this point
            int j = find_closest_point( p , pose_2 , nullptr );
            temp.push_back( j );

        }
        (*pose_2_clusters).push_back( temp );
        i++;
    }
    
    //now we remove duplicate points in the corresponding cluster for instances were 2 
    // points mapped to the same point in the deformed pose
    for(auto &u : *pose_2_clusters )
    {
        std::sort( u.begin(), u.end() );
        u.erase( unique( u.begin(), u.end() ), u.end() );
    }


}

/*******************************************************************************************************/

// recompute cluster Transforms
void ssdr::recompute_RT( const std::vector<Point>& pose_1,
                         const std::vector<Point>& pose_2, 
                         std::vector<Eigen::Matrix2d>* Rs,
                         std::vector<Eigen::Vector2d>* Ts,
                         const std::vector<std::vector<int>>& pose_1_clusters,
                         const std::vector<std::vector<int>>& pose_2_clusters )
{

    for( int j = 0; j < pose_1_clusters.size() ; j++ )
    {

        std::vector<Point> cluster1;
        for( int i : pose_1_clusters.at(j) )
        {
            cluster1.push_back( pose_1.at(i) );
        }

        std::vector<Point> cluster2;
        for( int i : pose_2_clusters.at(j) )
        {
            cluster2.push_back( pose_2.at(i) );
        }
        Eigen::Matrix2d R;
        Eigen::Vector2d T;

        compute_cluster_RT( cluster1, cluster2, &R, &T );

        (*Rs).push_back( R );
        (*Ts).push_back( T );
    }
}


/*******************************************************************************************************/

/*
std::vector<Eigen::Vector2d> ssdr::transform_sample_points( )
{
    // get the rest pose coordinates
    Eigen::Vector2d point(0,0);

    std::vector<Eigen::Vector2d> transformed_points;

    int label_num = 0;

    for( int i = 0 ; i < rest_pose_samples ; i++ )
    {
        point(0) = rest_pose_samples(i,0);
        point(1) = rest_pose_samples(i,1);
        
        label_num = vertex_label_index.at(i);

        // calculate the transformed rest pose
        point = pair_R * point + pair_T;
        
        transformed_points.push_back( point );

    }
    return transformed_points;
}

void ssdr::draw_curves( )
{

}

*/
}
