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
#include <cairo.h>
#include <cairo-svg.h> 
#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "solve_it.hpp"

/******************************************************************************************/

int smoothFn_anydir( int p1, int p2, int l1, int l2, void* matrix_term )
{
    if( l1 == l2 )
    {
        return 0;
    }
    std::vector<std::vector<Point>>* matrix_term_2 = static_cast<std::vector<std::vector<Point>>*> (matrix_term);
    Point p1_prime = (matrix_term_2->at(p1)).at(l1); 
    Point p1_double_prime = (matrix_term_2->at(p1)).at(l2);
    
    Point p2_prime =  (matrix_term_2->at(p2)).at(l2);
    Point p2_double_prime =(matrix_term_2->at(p2)).at(l1);

    double d1 = 10 *(p1_prime - p1_double_prime).norm( );
    double d2 = 10 *(p2_prime - p2_double_prime).norm( );
//    printf(" &&&&&&&&&&&&&&&&&&&&&&& d1 = %lf , d2 = %lf\n",d1,d2);
    //return static_cast<int>((d1 + d2));
    return ( (d1 + d2 > 100) ?  100 : static_cast<int>((d1 + d2)) );
}
/******************************************************************************************/
/******************************************************************************************/

int smoothFn_anydir_2( int p1, int p2, int l1, int l2, void* matrix_term )
{
    if( l1 == l2 )
    {
        return 0;
    }
    std::vector<std::vector<Point>>* matrix_term_2 = static_cast<std::vector<std::vector<Point>>*> (matrix_term);
    Point p1_prime = (matrix_term_2->at(p1)).at(l1); 
    Point p1_double_prime = (matrix_term_2->at(p1)).at(l2);
    
    Point p2_prime =  (matrix_term_2->at(p2)).at(l2);
    Point p2_double_prime =(matrix_term_2->at(p2)).at(l1);

    double d1 = (p1_prime - p1_double_prime).norm( );
    double d2 = (p2_prime - p2_double_prime).norm( );

    double x = (d1) + (d2);
    double G = exp ( (-1 * x * x)/1000 ); 
    int cost = static_cast<int> (100 * ( 1 - G ));
    return cost;
}
/******************************************************************************************/
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
    srand (time(NULL));
    // sample the curves from the rest pose to use for solving 
    // the correspondance between the rest pose curves and the
    // deformed pose curves
    // this also computes the tangents at the sampled points
    // in addition to what curve it belongs to
    rp_SamplePoints = sample_curves( rp_CurveEndPoints, 
                                     rp_CurveMiddlePoints,
                                     rp_Curves, 
                                     4, 
                                     &rp_SamplePoint_tg,
                                     &rp_SamplePoint_Curves );


    
    // sample the doformed pose curves and find their tangets and
    // normalize them like above
    dfp_SamplePoints = sample_curves( dfp_CurveEndPoints, 
                                      dfp_CurveMiddlePoints,
                                      dfp_Curves, 
                                      4, 
                                      &dfp_SamplePoint_tg,
                                      &dfp_SamplePoint_Curves );


    Point temp_p1 =  calculate_scaling_factor( rp_SamplePoints );
    Point temp_p2 =  calculate_scaling_factor( dfp_SamplePoints );
    // this scaling factor specifies what the max x and y can be
    int sf_max = 200;
    Point scaling_p(0, 0);
    scaling_p(0) = ( (temp_p1(0) > temp_p2(0)) ? temp_p1(0) : temp_p2(0) ) / sf_max;
    scaling_p(1) = ( (temp_p1(1) > temp_p2(1)) ? temp_p1(1) : temp_p2(1) ) / sf_max;

    double max_diagnal = sf_max * sqrt(2);

    Point dividing_p(0, 0);
    dividing_p(0) = 1.0 / scaling_p(0);
    dividing_p(1) = 1.0 / scaling_p(1);

    global_num_points = rp_SamplePoints.size();

    rp_SamplePoints = scale_samples( rp_SamplePoints, dividing_p );
    rp_SamplePoint_tg = scale_samples( rp_SamplePoint_tg, dividing_p );
    dfp_SamplePoints = scale_samples( dfp_SamplePoints, dividing_p );
    dfp_SamplePoint_tg = scale_samples( dfp_SamplePoint_tg, dividing_p );

    std::cout << "\n dividing factor\n" << dividing_p << std::endl; 
    std::cout << "\n scaling factor\n" << scaling_p << std::endl; 
    std::cout << " scaled coord ************************************* " << std::endl;
    for( auto& boo : rp_SamplePoints )
    {
        std::cout << boo << std::endl;
    }
    std::cout << " scaled tgs ************************************* " << std::endl;
    for( auto& boo : rp_SamplePoint_tg )
    {
        std::cout << boo << std::endl;
    }
    std::cout << " scaled coord ************************************* " << std::endl;
    for( auto& boo : dfp_SamplePoints )
    {
        std::cout << boo << std::endl;
    }
    std::cout << " scaled tgs ************************************* " << std::endl;
    for( auto& boo : dfp_SamplePoint_tg )
    {
        std::cout << boo << std::endl;
    }
    std::cout << " scaled coord ************************************* " << std::endl;
    
    // makes the tangents of rhe rest samples  unit length
    rp_SamplePoint_tg =  normalize_tangents( rp_SamplePoint_tg );
    dfp_SamplePoint_tg = normalize_tangents( dfp_SamplePoint_tg );

    // since we know what curve each sample point belongs to we will
    // compute a matrix that specified which 2 points are considered
    // neighbors i.e. belong to curves that intersect
    Eigen::MatrixXd rp_neighbor_curves = set_curve_neighborhood( rp_Curves );
    Eigen::MatrixXd dfp_neighbor_curves = set_curve_neighborhood( dfp_Curves );
    
    std::vector<std::vector<int>> rp_nb_local(global_num_points);
	// holds distances
    std::vector<std::vector<double>> rp_nb_local_d(global_num_points);
    
    std::vector<std::vector<int>> dfp_nb_local(global_num_points);
	// holds distances
    std::vector<std::vector<double>> dfp_nb_local_d(global_num_points);

    set_nb_local(rp_nb_local , rp_nb_local_d, rp_SamplePoints, rp_SamplePoint_Curves, rp_neighbor_curves, max_diagnal );
    set_nb_local(dfp_nb_local , dfp_nb_local_d, dfp_SamplePoints, dfp_SamplePoint_Curves, dfp_neighbor_curves, max_diagnal );
    
/*    for( auto& vv : rp_nb_local )
    {
        for( int ii : vv )
        {
            std::cout << ii << "  ";
        }
        std::cout << "\n" ;
    } */
    // Compute n x n initial transformations from rest pose to 
    // deformed pose . then use the invers and negated versions
    // of these for transforming deformed pose sample points to rest pose


    // will hold the n x n Transformations
    std::vector<Eigen::Matrix2d> Rs;
    std::vector<Eigen::Vector2d> Ts;
    
    // will hold the n x n Transformations
    std::vector<Eigen::Matrix2d> all_possible_Rs;
    std::vector<Eigen::Vector2d> all_possible_Ts;

    // will hold the n x n Transformations ( inv(R) _T )
    std::vector<Eigen::Matrix2d> Rs_inv;
    std::vector<Eigen::Vector2d> Ts_neg;
    
    // we compute n x n candidate R,T by pairing each point of the 
    // rest pose with each point in the deformed pose
    initialize_transformations( rp_SamplePoints, dfp_SamplePoints, 
                                rp_SamplePoint_tg, dfp_SamplePoint_tg,
                                &Rs, &Ts );

    // based on the transformations from rest to deformed pose calculate
    // the trasformation the other way -->  R.inv * ( p_dfp - T ) = p_rp
    for( int l = 0 ; l < Rs.size() ;l++ )
    {
        all_possible_Rs.push_back( Rs.at(l) );
        all_possible_Ts.push_back( Ts.at(l) );
        Rs_inv.push_back( (Rs.at(l)).inverse() );
        Ts_neg.push_back( (Rs.at(l)).inverse() * (-1 * Ts.at(l)) );
    }


    std::vector<std::vector<int>> pose_1_clusters;
    std::vector<std::vector<int>> pose_2_clusters;
    
    std::vector<std::vector<int>> reverse_pose_1_clusters;
    std::vector<std::vector<int>> reverse_pose_2_clusters;

    int iter = 0;
    std::vector<int> label_assignment;
    std::vector<int> label_assignment_2;

    // now having n x n labels and having an initial label assignment 
    // we perform a bidirectional multilevel optimization
    do{
        
        // calculate a label assignment seperately in each direction

        int* data_forward = set_data_term( rp_SamplePoints, 
                                           dfp_SamplePoints, 
                                           Rs, Ts );
/*
            printf(" data term forward \n ");
            for( int j = 0 ; j <  rp_SamplePoints.size() ; j++ )
            {
                for( int i = 0 ; i < Rs.size() ; i++ )
                {
                    printf(" %d  ", data_forward[ j * Rs.size() + i ] );
                }
                printf("\n");
            } 
*/
        int* data_reverse = set_data_term( dfp_SamplePoints,
                                            rp_SamplePoints,
                                            Rs_inv, Ts_neg );
  /*          printf(" data term inverse \n ");
            for( int j = 0 ; j <  dfp_SamplePoints.size() ; j++ )
            {
                for( int i = 0 ; i < Rs_inv.size() ; i++ )
                {
                    printf(" %d  ", data_reverse[ j * Rs_inv.size() + i ] );
                }
                printf("\n");
            } 
*/
        int *smooth = set_smoothness_term( Rs.size() );

        int *numNeighbors_rp;
        int **neighborsIndexes_rp;
        int** neighborsWeights_rp;
        
        int *numNeighbors_dfp;
        int **neighborsIndexes_dfp;
        int** neighborsWeights_dfp;

        compute_neighborhood_info( rp_SamplePoint_Curves, 
                                   rp_neighbor_curves, 
                                   &numNeighbors_rp, 
                                   &neighborsIndexes_rp, 
                                   &neighborsWeights_rp );
           
        compute_neighborhood_info( dfp_SamplePoint_Curves, 
                                   dfp_neighbor_curves, 
                                   &numNeighbors_dfp, 
                                   &neighborsIndexes_dfp, 
                                   &neighborsWeights_dfp );
        
/*        compute_neighborhood_info_2( rp_nb_local,
									 rp_nb_local_d,
                                     &numNeighbors_rp, 
                                     &neighborsIndexes_rp, 
                                     &neighborsWeights_rp,
									 max_diagnal );

        compute_neighborhood_info_2( dfp_nb_local,
									 dfp_nb_local_d,
                                     &numNeighbors_dfp, 
                                     &neighborsIndexes_dfp, 
                                     &neighborsWeights_dfp,
									 max_diagnal );

            printf("\n\n deformed neighbor  \n\n");
            for( int i = 0; i < dfp_SamplePoints.size() ; i++ )
            {
                printf(" \nnum neighbors %d = %d \n",i, numNeighbors_dfp[i] );
                for( int j = 0 ; j < numNeighbors_dfp[i] ; j++ )
                {
                    printf(" %d  ",neighborsIndexes_dfp[i][j]);
                }
            }
            for( int i = 0; i < dfp_SamplePoints.size() ; i++ )
            {
                printf(" \nnum neighbors %d = %d \n",i, numNeighbors_dfp[i] );
                for( int j = 0 ; j < numNeighbors_dfp[i] ; j++ )
                {
                    printf(" %d  ",neighborsWeights_dfp[i][j]);
                }
            }

            printf("\n\n rest neighbor  \n\n");
            for( int i = 0; i < rp_SamplePoints.size() ; i++ )
            {
                printf(" \nnum neighbors %d = %d \n",i, numNeighbors_rp[i] );
                for( int j = 0 ; j < numNeighbors_rp[i] ; j++ )
                {
                    printf(" %d  ",neighborsIndexes_rp[i][j]);
                }
            }
            for( int i = 0; i < rp_SamplePoints.size() ; i++ )
            {
                printf(" \nnum neighbors %d = %d \n",i, numNeighbors_rp[i] );
                for( int j = 0 ; j < numNeighbors_rp[i] ; j++ )
                {
                    printf(" %d  ",neighborsWeights_rp[i][j]);
                }
            }
*/
        int num_points = rp_SamplePoints.size();
        int num_labels = Rs.size();

        std::vector<int> label_assign_forward;
        std::vector<int> label_assign_reverse;

        label_assign_forward = assign_labels( num_points, num_labels,
                                              data_forward, smooth,
                                              numNeighbors_rp,
                                              neighborsIndexes_rp, 
                                              neighborsWeights_rp );

        label_assign_reverse = assign_labels( num_points, num_labels, 
                                              data_reverse, smooth,
                                              numNeighbors_dfp, 
                                              neighborsIndexes_dfp, 
                                              neighborsWeights_dfp );
       
/*
        std::vector<std::vector<Point>> smooth_matrix_forward;
        std::vector<std::vector<Point>> smooth_matrix_reverse;
        std::vector<std::vector<Point>> smooth_matrix_bidir;
       // for every pair of label and point calculate the transformed point
       // and store it in the smoothness matrix M(Point # , label #) 
       for( int p1_ind = 0 ; p1_ind < global_num_points ; p1_ind++ )
       {
           // get points from rest and defrmed pose
           Point p1_coord = rp_SamplePoints.at(p1_ind);
           Point p2_coord = dfp_SamplePoints.at(p1_ind);
           std::vector<Point> temp_v_1;
           std::vector<Point> temp_v_2;

           for( int l1=0 ; l1 < Rs.size(); l1++ )
           {
               temp_v_1.push_back( Rs.at(l1) * p1_coord + Ts.at(l1) );
               temp_v_2.push_back( Rs_inv.at(l1) * p2_coord + Ts_neg.at(l1) );

           }
           smooth_matrix_forward.push_back(temp_v_1);
           smooth_matrix_reverse.push_back(temp_v_2);
           smooth_matrix_bidir.push_back(temp_v_1);
       }
       // setting the rest of the bidirectional matrix
       for( auto q_coord : dfp_SamplePoints )
       {
           std::vector<Point> temp_v_2;
           for( int l1=0 ; l1 < Rs.size(); l1++ )
           {
               temp_v_2.push_back( Rs_inv.at(l1) * q_coord + Ts_neg.at(l1) );

           }
           smooth_matrix_bidir.push_back(temp_v_2);

       }
  */    
      /*  global_Rs = Rs;
        global_Ts = Ts;
        global_SamplePoints = rp_SamplePoints; 
*/
       // int (ssdr::*fn)(int, int, int, int) = &ssdr::smoothFn_unidirectional;
/*        label_assign_forward = assign_labels_2( num_points, num_labels,
                                                data_forward, &::smoothFn_anydir_2,
                                                static_cast<void *>(&smooth_matrix_forward),
                                                numNeighbors_rp,
                                                neighborsIndexes_rp, 
                                                neighborsWeights_rp );

        global_Rs = Rs_inv;
        global_Ts = Ts_neg;
        global_SamplePoints = dfp_SamplePoints; 

        label_assign_reverse = assign_labels_2( num_points, num_labels, 
                                                data_reverse, &::smoothFn_anydir_2,
                                                static_cast<void *>(&smooth_matrix_reverse),
                                                numNeighbors_dfp, 
                                                neighborsIndexes_dfp, 
                                                neighborsWeights_dfp );
        global_Rs = Rs;
        global_Ts = Ts;
        global_Rs_2 = Rs_inv;
        global_Ts_2 = Ts_neg; */
 /*       printf("\nlabel assignment\n");
        for( int l : label_assign_forward )
        {
            printf(" %d   ",l);
        }
        
        printf("\nlabel assignment reverse\n");
        for( int l : label_assign_reverse )
        {
            printf(" %d   ",l);
        }
   */     
        // having a label assignment separately for each direction
        // the bidirectional label assignment is computed 

        // combine data terms by considering the first n points the 
        // points of the rest pose and the second n points as the points
        // of the deformed pose
        int* data_terms = combine_data_terms( data_forward, data_reverse,
                                         num_points, num_labels );
        
/*        printf(" data term\n ");
        for( int j = 0 ; j <  2*dfp_SamplePoints.size() ; j++ )
        {
            for( int i = 0 ; i < Rs.size() ; i++ )
            {
                printf(" %d  ", data_terms[ j * Rs.size() + i ] );
            }
            printf("\n");
        } 
*/
        // Having the numNeighbors, neighborsIndexes, neighborsWeights
        // separately for each direction, the bidirectional 
        // numNeighbors, neighborsIndexes, neighborsWeights is computed

        
        // having a label assignment in each direction separately 
        // the transformations are applied and the corresponding point 
        // to each sample point is found. rest pose sample points are
        // considered the first n spots and the second n spots are
        // the deformed pose sample points . each spot holds a vector
        // that has the indeces of the points which this sample point
        // corresponds to
        std::vector<std::vector<int>> nbh_temp;
        nbh_temp  = bidirectional_correspondance( rp_SamplePoints, 
                                                  dfp_SamplePoints, 
                                                  Rs, Ts,
                                                  label_assign_forward, 
                                                  label_assign_reverse );
        combine_neighborhood_info( &numNeighbors_rp,
                                   &neighborsIndexes_rp, 
                                   &neighborsWeights_rp,
                                   numNeighbors_dfp,
                                   neighborsIndexes_dfp, 
                                   neighborsWeights_dfp,
                                   nbh_temp, 
                                   num_points );

  /*      combine_neighborhood_info_2( &numNeighbors_rp,
                                   &neighborsIndexes_rp, 
                                   &neighborsWeights_rp,
                                   numNeighbors_dfp,
                                   neighborsIndexes_dfp, 
                                   neighborsWeights_dfp,
                                   nbh_temp, 
                                   num_points );   */
     /*   for( int i = 0; i < 2*rp_SamplePoints.size() ; i++ )
        {
            printf(" \nnum neighbors %d = %d \n",i, numNeighbors_rp[i] );
            for( int j = 0 ; j < numNeighbors_rp[i] ; j++ )
            {
                printf(" %d  ",neighborsIndexes_rp[i][j]);
            }
        }
        for( int i = 0; i < 2*rp_SamplePoints.size() ; i++ )
        {
            printf(" \nnum neighbors %d = %d \n",i, numNeighbors_rp[i] );
            for( int j = 0 ; j < numNeighbors_rp[i] ; j++ )
            {
                printf(" %d  ",neighborsWeights_rp[i][j]);
            }
        }
      
*/
        label_assignment = assign_labels( 2*num_points, num_labels, 
                                          data_terms, smooth,
                                          numNeighbors_rp, 
                                          neighborsIndexes_rp, 
                                          neighborsWeights_rp );
        
       /* label_assignment = assign_labels_2( 2*num_points, num_labels, 
                                            data_terms, &::smoothFn_anydir_2,
                                            static_cast<void *>(&smooth_matrix_bidir),
                                            numNeighbors_rp, 
                                            neighborsIndexes_rp, 
                                            neighborsWeights_rp );*/
        // having computed the bidirectional label assignment
        // the clusters are derived
        
        // the labels assigned to the first n points are the labels 
        // assigned to the rest pose

        std::vector<int> label_assignment_half;
        for( int f = 0 ; f < num_points ; f++ )
        {
            label_assignment_half.push_back( label_assignment.at(f) );
        }

        printf("\n\nlabel assignment half\n");
        for( int l : label_assignment_half )
        {
            printf(" %d   ",l);
        }
        
        std::vector<int> label_assignment_half_2;
        for( int f = num_points ; f < 2*num_points ; f++ )
        {
            label_assignment_half_2.push_back( label_assignment.at(f) );
        }

        printf("\n\nlabel assignment 2nd half\n");
        for( int l : label_assignment_half_2 )
        {
            printf(" %d   ",l);
        }
        // the rest pose points having the same label are put into
        // one cluster and the rest pose clusters are derived 

        pose_1_clusters.clear();
        std::vector<int> label_numbers;

        set_restpose_clusters( label_assignment_half, num_points, 
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
        
        reverse_pose_1_clusters.clear();
        std::vector<int> reverse_label_numbers;

        set_restpose_clusters( label_assignment_half_2, num_points, 
                               &reverse_pose_1_clusters, &reverse_label_numbers); 

        uu = 0;
        for( auto& rr : reverse_pose_1_clusters )
        {
            printf("\nreverse cluster label %d\n",reverse_label_numbers.at(uu));
            
            for( auto& oo : rr )
            {
                printf(" %d  ", oo );
            }
            uu++;
        }

        pose_2_clusters.clear();

        // given the restpose clusters the corresponding deformed pose 
        // clusters are computed
        compute_corresponding_clusters( rp_SamplePoints, 
                                        dfp_SamplePoints, 
                                        label_numbers, 
                                        Rs, Ts,
                                        pose_1_clusters, 
                                        &pose_2_clusters );
        
        for( auto& rr : pose_2_clusters )
        {
            printf("\ncluster 2\n");
            for( auto& oo : rr )
            {
                printf(" %d  ", oo );
            }
        }
        
        reverse_pose_2_clusters.clear();
        Rs_inv.clear();
        Ts_neg.clear();

        // based on the transformations from rest to deformed pose calculate
        // the trasformation the other way -->  R.inv * ( p_dfp - T ) = p_rp
        for( int l = 0 ; l < Rs.size() ;l++ )
        {
            Rs_inv.push_back( (Rs.at(l)).inverse() );
            Ts_neg.push_back( (Rs.at(l)).inverse() * (-1 * Ts.at(l)) );
        }

        // given the restpose clusters the corresponding deformed pose 
        // clusters are computed
        compute_corresponding_clusters( dfp_SamplePoints, 
                                        rp_SamplePoints, 
                                        reverse_label_numbers, 
                                        Rs_inv, Ts_neg,
                                        reverse_pose_1_clusters, 
                                        &reverse_pose_2_clusters );
        
        for( auto& rr : reverse_pose_2_clusters )
        {
            printf("\nreverse cluster 2\n");
            for( auto& oo : rr )
            {
                printf(" %d  ", oo );
            }
        }
        printf("\n\n");
      
        std::vector<Point> trans_points = transform_sample_points( rp_SamplePoints,
                                                                   Rs, Ts, 
                                                                   label_assignment_half); 

        std::vector<Point> reverse_trans_points = transform_sample_points( dfp_SamplePoints,
                                                                           Rs_inv, Ts_neg, 
                                                                           label_assignment_half_2); 
      
        std::vector<Point> rp_SP = scale_samples( rp_SamplePoints, scaling_p );
        std::vector<Point> rp_SP_tg = scale_samples( rp_SamplePoint_tg, scaling_p );
        std::vector<Point> dfp_SP = scale_samples( dfp_SamplePoints, scaling_p );
        std::vector<Point> dfp_SP_tg = scale_samples( dfp_SamplePoint_tg, scaling_p );
        trans_points = scale_samples( trans_points, scaling_p );
        reverse_trans_points = scale_samples( reverse_trans_points, scaling_p );
/*
        output_svg( rp_SamplePoints, trans_points, 
                    dfp_SamplePoints,  pose_1_clusters, 
                    pose_2_clusters,"forward", iter  );
        
        output_svg( dfp_SamplePoints, reverse_trans_points, 
                    rp_SamplePoints,  reverse_pose_1_clusters, 
                    reverse_pose_2_clusters,"reverse", iter  );
*/
        output_svg( rp_SP, trans_points, 
                    dfp_SP,  pose_1_clusters, 
                    pose_2_clusters,"forward", iter  );
        
        output_svg( dfp_SP, reverse_trans_points, 
                    rp_SP,  reverse_pose_1_clusters, 
                    reverse_pose_2_clusters,"reverse", iter  );
        // now having clusters in the rest pose and the corresponding 
        // clusters the transformations are recomputed
        
        Rs.clear();
        Ts.clear();

        Rs_inv.clear();
        Ts_neg.clear();

        recompute_RT(  rp_SamplePoints, dfp_SamplePoints, 
                       &Rs, &Ts, 
                       pose_1_clusters, 
                       pose_2_clusters );
        
        // make a vector of inverse Rs for the bidirectional ICP
        for( int l = 0 ; l < Rs.size() ;l++ )
        {
            Rs_inv.push_back( (Rs.at(l)).inverse() );
            Ts_neg.push_back( (Rs.at(l)).inverse() * (-1 * Ts.at(l)) );
        }

        // insert n random labels
        for( int i = 0 ; i < num_points ; i++ )
        {
            int random_index = rand() % (num_points * num_points);
            std::cout << " random number " << i << " is " << random_index << std::endl;
            Rs.push_back( all_possible_Rs.at(random_index) );
            Ts.push_back( all_possible_Ts.at(random_index) );

            Rs_inv.push_back( (all_possible_Rs.at(random_index)).inverse() );
            Ts_neg.push_back( (all_possible_Rs.at(random_index)).inverse() *
                                (-1 * all_possible_Ts.at(random_index)) );
        }
        
    }while( (++iter) < 5 );
    std::vector<Point> trans_points_final;
    //solver logic
    // now having the clusters in the rest pose and the corresponding clusters
    // in deformed pose for each pair we solve for the weights
    for(int i = 0 ; i < pose_1_clusters.size() ; i++ )
    {
        for(int j = 0 ; j < pose_1_clusters.at(i).size() ; j++)
        {
            int index1 = pose_1_clusters.at(i).at(j);
            int index2 = pose_2_clusters.at(i).at(j);
            Point p_i = rp_SamplePoints.at(index1);
            Point v_i = dfp_SamplePoints.at(index2);
            auto solution_w = QPSolve_i( Rs, Ts, p_i, v_i );
            trans_points_final.push_back( solution_w );
            std::cout << " ################################\n" << solution_w << std::endl;
        }
    }
        std::vector<Point> rp_SP = scale_samples( rp_SamplePoints, scaling_p );
        std::vector<Point> dfp_SP = scale_samples( dfp_SamplePoints, scaling_p );
        trans_points_final = scale_samples( trans_points_final, scaling_p );
        output_svg_2( rp_SP, dfp_SP, trans_points_final, "final_result", 0  );

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

std::vector<std::vector<int>> ssdr::bidirectional_correspondance( const std::vector<Point>& pose_1,
                                                                  const std::vector<Point>& pose_2,
                                                                  const std::vector<Eigen::Matrix2d>& Rs,
                                                                  const std::vector<Eigen::Vector2d>& Ts,
                                                                  const std::vector<int>& label_assignment,
                                                                  const std::vector<int>& label_assignment_2 )
{
    //stores the number of labels
    int num_labels = Rs.size();

    //number of sample points
    int num_samples = pose_1.size();
    
    std::vector<std::vector<int>> ret( 2 * num_samples ); 
    
    // For each vertex A in the rest position 
    for( int i = 0; i < num_samples; i++)
    {
        int l_i1 = label_assignment.at( i );
        int l_i2 = label_assignment_2.at( i );

        Point p = Rs.at(l_i1) * pose_1.at(i) + Ts.at(l_i1);
        int index = find_closest_point( p , pose_2 , nullptr );

        (ret.at(i)).push_back( index + num_samples );
        (ret.at(index + num_samples)).push_back( i );

        p = (Rs.at(l_i2)).inverse() * ( pose_1.at(i) + (-1 * Ts.at(l_i2)) );
        index = find_closest_point( p , pose_1 , nullptr );

        (ret.at( i + num_samples )).push_back( index );
        (ret.at( index )).push_back( i + num_samples );
    }
    
    for( int i = 0; i < 2*num_samples; i++)
    {    
        std::sort( ret[i].begin(), ret[i].end() );
        ret[i].erase( unique( ret[i].begin(), ret[i].end() ), ret[i].end() );

    }
    
    return ret;
}

/*******************************************************************************************************/
int* ssdr::combine_data_terms( int* data, int* data_reverse ,int num_points, int num_labels )
{
    int* ret = new int[2 * num_points * num_labels];
    int i = 0;

    for( i = 0 ; i < num_points ; i++ )
    {
        for( int label_index = 0; label_index < num_labels; label_index++ )
        {
            ret[ i * num_labels + label_index ] = data[ i * num_labels + label_index ];
        }
    }
    for( i = 0 ; i < num_points ; i++ )
    {
        for( int label_index = 0; label_index < num_labels; label_index++ )
        {
            ret[ ( i + num_points ) * num_labels + label_index ] = data_reverse[ i * num_labels + label_index ];
        }
    }
    return ret;

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
            smooth[l1+l2*num_labels] = (l1-l2)*(l1-l2) < 1  ? 0 : 10;
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
            //(*neighborsWeights)[i][j] = 1;
            (*neighborsWeights)[i][j] = 10;
        }
    }
}

/*******************************************************************************************************/
void ssdr::combine_neighborhood_info( int** numNeighbors, int*** neighborsIndexes, int*** neighborsWeights,
                                      int* numNeighbors_dfp, int** neighborsIndexes_dfp, int** neighborsWeights_dfp,
                                      const std::vector<std::vector<int>>& nbh_temp, int num_points )
{
    int *ret_numNeighbors = new int[ 2 * num_points ];
    int **ret_neighborsIndexes = new int*[ 2 * num_points ];
    int **ret_neighborsWeights = new int*[ 2 * num_points ];
    std::vector< std::vector<int> > temp( 2 * num_points );
    std::vector< std::vector<int> > temp_2( 2 * num_points );

    // combine number of neighbors
    for( int i = 0 ; i < num_points ; i++ )
    {
        // initially we set the number of neighbors of each vertex 
        // this does not include the new bidirectional neighbors 
        ret_numNeighbors[i] = (*numNeighbors)[i];
        ret_numNeighbors[ i + num_points ] = (numNeighbors_dfp)[i]; 
        
        // first push_back the neighbors in each direction separately
        for( int k = 0 ; k < ret_numNeighbors[i] ; k++ )
        {
            temp[i].push_back( (*neighborsIndexes)[i][k] );
        }
        
        for( int k = 0 ; k < ret_numNeighbors[ i + num_points ] ; k++ )
        {
            temp[ i + num_points ].push_back( (neighborsIndexes_dfp)[i][k] + num_points );
        }

        for( int k : nbh_temp.at(i) )
        {
            temp[i].push_back( k );
        }
        
        for( int k : nbh_temp.at( i + num_points ) )
        {
            temp[ i + num_points ].push_back( k );
        }

        std::sort( temp[ i + num_points ].begin(), temp[ i + num_points ].end() );
        temp[ i + num_points ].erase( unique( temp[ i + num_points ].begin(), temp[ i + num_points ].end() ), temp[ i + num_points ].end() );
        
        std::sort( temp[i].begin(), temp[i].end() );
        temp[i].erase( unique( temp[i].begin(), temp[i].end() ), temp[i].end() );
        
    }

    //combine neighbor Indexes and weights
    for( int i = 0 ; i < 2 * num_points ; i++ )
    {
        ret_neighborsIndexes[i] = new int[ (temp.at(i)).size() ];
        ret_neighborsWeights[i] = new int[ (temp.at(i)).size() ];

        ret_numNeighbors[i] = (temp.at(i)).size();

        for( int j = 0 ; j < ret_numNeighbors[i] ; j++ )
        {
            ret_neighborsIndexes[i][j] = (temp.at(i)).at(j);
            ret_neighborsWeights[i][j] = 10;
            for( int k = 0 ; k < (nbh_temp.at(i)).size() ; k++ )
            {
                if( ret_neighborsIndexes[i][j] == (nbh_temp.at(i)).at(k) )
                {
                    ret_neighborsWeights[i][j] = 1;
                }
            }
        }

    }
    
    *numNeighbors = ret_numNeighbors;
    *neighborsIndexes = ret_neighborsIndexes;
    *neighborsWeights = ret_neighborsWeights;

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
        
        printf("\nBefore optimization energy is %lld\n",gc->compute_energy());
        //gc->compute_energy();
		gc->expansion(2);
        gc->swap(2);
        
        printf("\nAfter optimization energy is %lld\n",gc->compute_energy());
        //gc->compute_energy();

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
std::vector<int> ssdr::assign_labels_2( int num_points,
                                        int num_labels,
                                        int *data_terms,
                                        GCoptimization::SmoothCostFnExtra fn,
                                        void* smooth_M,
                                        int *numNeighbors, 
                                        int **neighborsIndexes, 
                                        int** neighborsWeights )
{
    std::vector<int> result;
    
    try{

		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph( num_points , num_labels );
		gc->setDataCost( data_terms );
		gc->setSmoothCost( fn, smooth_M );

        //pass in all neighbor information at once
        gc->setAllNeighbors( numNeighbors , neighborsIndexes , neighborsWeights );
        
        printf("\nBefore optimization energy is %lld\n",gc->compute_energy());
        //gc->compute_energy();
		gc->expansion(2);
        gc->swap(2);
        
        printf("\nAfter optimization energy is %lld\n",gc->compute_energy());
        //gc->compute_energy();

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
   /* 
    //now we remove duplicate points in the corresponding cluster for instances were 2 
    // points mapped to the same point in the deformed pose
    for(auto &u : *pose_2_clusters )
    {
        std::sort( u.begin(), u.end() );
        u.erase( unique( u.begin(), u.end() ), u.end() );
    }
*/

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

// apply the appropriate labels to the points and 
std::vector<Point> ssdr::transform_sample_points( const std::vector<Point>& original_points, 
                                                  const std::vector<Eigen::Matrix2d>& Rs,
                                                  const std::vector<Eigen::Vector2d>& Ts,
                                                  const std::vector<int>& label_numbers)
{
    std::vector<Point> transformed_points;
    int i = 0;
    int index = 0;
    for( auto& p : original_points )
    {
        index = label_numbers.at(i);
        const Eigen::Matrix2d& R = Rs.at(index); 
        const Eigen::Vector2d& T = Ts.at(index);
        Point q = R * p + T;
        transformed_points.push_back( q );
      //  std::cout << "\nR\n" << R << "\nT\n" << T << "\npoint\n" << p << "\nother\n" << q << std::endl;
        i++;
        index++;
    }
    return transformed_points;
}

/*******************************************************************************************************/

void ssdr::output_svg( const std::vector<Point>& pose_1_samples, 
                       const std::vector<Point>& pose_2_samples,
                       const std::vector<Point>& pose_2_samples_B, 
                       const std::vector<std::vector<int>>& pose_1_clusters,
                       const std::vector<std::vector<int>>& pose_2_clusters,
                       char const *file_name,
                       int index )
{
    std::vector<Eigen::Vector3d> colors;
    Eigen::Vector3d temp = Eigen::Vector3d( 0.5, 0.5, 0.5 );
    colors.push_back( temp );
    temp = Eigen::Vector3d( 1, 0, 0 );
    colors.push_back( temp );
    //colors.push_back( Eigen::Vector3d( 0.5, 0 , 0) );
    //colors.push_back( Eigen::Vector3d( 1, 1, 0 ) );
    colors.push_back( Eigen::Vector3d( 0.5, 0.5, 0 ) );
    colors.push_back( Eigen::Vector3d( 0, 1, 0 ) );
    colors.push_back( Eigen::Vector3d( 0, 0.5, 0 ) );
    colors.push_back( Eigen::Vector3d( 0, 1, 1 ) );
    colors.push_back( Eigen::Vector3d( 0, 0, 1 ) );
    colors.push_back( Eigen::Vector3d( 1, 0, 1 ) );
    colors.push_back( Eigen::Vector3d( 0.5, 0, 0.5 ) );
    
    cairo_t *cr;
    cairo_surface_t *surface;
    cairo_pattern_t *pattern;
    int x,y;
    char buffer[16]; // make sure it's big enough
    snprintf(buffer, sizeof(buffer), "%s_%d.svg", file_name, index);
    
    surface =
      (cairo_surface_t *)cairo_svg_surface_create(buffer, 1200.0, 1200.0);
    cr = cairo_create(surface);
    
    // draw the original rest pose image curves in black ( 0 , 0 , 0 )
    
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width( cr, 5 );

    Point p0, p1, p2, p3;
    
    for( auto & curve : rp_Curves )
    {
        //end points
        p0 = rp_CurveEndPoints.at(curve[0]);
        p3 = rp_CurveEndPoints.at(curve[3]);
        
        //tangent control points
        p1 = rp_CurveMiddlePoints.at(curve[1]);
        p2 = rp_CurveMiddlePoints.at(curve[2]);
        
        cairo_move_to( cr, p0(0), p0(1) );
        cairo_curve_to( cr, p1(0), p1(1), p2(0), p2(1), p3(0), p3(1) );
        cairo_stroke( cr );

    }

    // draw the original deformed pose image curves in navy ( 0 , 0 , 0.5 )
    
    cairo_set_source_rgb(cr, 0, 0, 0.5);
    cairo_set_line_width( cr, 5 );
    
    for( auto & curve : dfp_Curves )
    {
        //end points
        p0 = dfp_CurveEndPoints.at(curve[0]);
        p3 = dfp_CurveEndPoints.at(curve[3]);
        
        //tangent control points
        p1 = dfp_CurveMiddlePoints.at(curve[1]);
        p2 = dfp_CurveMiddlePoints.at(curve[2]);
        
        cairo_move_to( cr, p0(0), p0(1) );
        cairo_curve_to( cr, p1(0), p1(1), p2(0), p2(1), p3(0), p3(1) );
        cairo_stroke( cr );

    }
    
    //draw the deformed pose samples
    cairo_set_source_rgb(cr, 0, 0, 0.5);
    cairo_set_line_width( cr, 3 );
    for( auto& p : pose_2_samples_B )
    {
        cairo_arc(cr, p(0), p(1), 5, 0, 2*M_PI);
        cairo_stroke( cr );
    }

    // draw the rest pose sample point's clusters
    
    cairo_set_line_width( cr, 3 );
    int j = 0;
    for( auto& cluster : pose_1_clusters )
    {
        Eigen::Vector3d c = colors.at( j );
        cairo_set_source_rgb(cr, c(0), c(1), c(2));

        for( auto& index : cluster )
        {
            Point p = pose_1_samples.at( index );
            cairo_arc(cr, p(0), p(1), 8, 0, 2*M_PI);
            cairo_stroke( cr );

        }
        j++;
    }

    // draw the rest pose sample point's clusters
    
    cairo_set_line_width( cr, 3 );
    j = 0;
    for( auto& cluster : pose_1_clusters )
    {
        Eigen::Vector3d c = colors.at( j );
        //cairo_set_source_rgb(cr, c(0), c(1), c(2));
        for( auto& index : cluster )
        {
            cairo_set_source_rgb( cr, 0 , 0, 0 );
            Point p = pose_2_samples.at( index );
            cairo_rectangle(cr, p(0), p(1), 12, 12 );
            cairo_stroke_preserve(cr);
            cairo_set_source_rgb(cr, c(0), c(1), c(2));
            cairo_fill( cr );

        }
        j++;
    }
    // draw the deformed pose sample point's clusters
    
    cairo_set_line_width( cr, 6 );
    j = 0;
    for( auto& cluster : pose_2_clusters )
    {
        Eigen::Vector3d c = colors.at( j );
        cairo_set_source_rgb(cr, c(0), c(1), c(2));

        for( auto& index : cluster )
        {
            Point p = pose_2_samples_B.at( index );
            cairo_arc(cr, p(0), p(1), 12, 0, 2*M_PI);
            cairo_stroke( cr );

        }
        j++;
    }

    /*cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width( cr, 5 );
    for( auto& p : points )
    {
        cairo_arc(cr, p(0), p(1), 5, 0, 2*M_PI);
        cairo_stroke( cr );
    }*/
    cairo_destroy (cr);
    cairo_surface_destroy (surface);

}

/*********************************************************************************************/

// iterates through the list of points and finds the max x annd max y values and returns them
// as a Point
Point ssdr::calculate_scaling_factor( const std::vector<Point>& point_coordinates )
{
    // find max x and max y
    int max_index_x = 0;
    int max_index_y = 0;
    int i = -1;
    Point max_point(0 , 0);
   
    for( const auto& p : point_coordinates )
    {
        if( max_point(0) < p(0) )
        {
            max_point(0) = p(0);
            max_index_x = i;
        }
        if( max_point(1) < p(1) )
        {
            max_point(1) = p(1);
            max_index_y = i;
        }
        i++;
    }
    return max_point;
}

/***************************************************************************************/

std::vector<Point> ssdr::scale_samples( const std::vector<Point>& point_coordinates, 
                                        const Point& scale_factors )
{
    std::vector<Point> norm_coord;
    for( const auto& p : point_coordinates )
    {
        Point q(p(0) * scale_factors(0), p(1) * scale_factors(1));
        
        norm_coord.push_back( q );
    }
    return norm_coord;

}

/******************************************************************************************/

int ssdr::smoothFn_unidirectional( int p1, int p2, int l1, int l2 )
{
    Point p1_prime = global_Rs.at(l1) * global_SamplePoints.at(p1) + global_Ts.at(l1); 
    Point p1_double_prime = global_Rs.at(l2) * global_SamplePoints.at(p1) + global_Ts.at(l2); 
    
    Point p2_prime = global_Rs.at(l2) * global_SamplePoints.at(p2) + global_Ts.at(l2); 
    Point p2_double_prime = global_Rs.at(l1) * global_SamplePoints.at(p2) + global_Ts.at(l1); 

    double d1 = (p1_prime - p1_double_prime).norm( );
    double d2 = (p2_prime - p2_double_prime).norm( );
    
    return static_cast<int>((d1 + d2)/2);

}

/******************************************************************************************/

int ssdr::smoothFn_bidirectional( int p1, int p2, int l1, int l2 )
{
    Point p1_prime; 
    Point p1_double_prime; 
    
    Point p2_prime; 
    Point p2_double_prime; 

    if( p1 < global_num_points )
    {
        p1_prime = global_Rs.at(l1) * rp_SamplePoints.at(p1) + global_Ts.at(l1); 
        p1_double_prime = global_Rs.at(l2) * rp_SamplePoints.at(p1) + global_Ts.at(l2); 
    }
    else
    {
        p1_prime = global_Rs_2.at(l1) * dfp_SamplePoints.at(p1 - global_num_points) 
                        + global_Ts_2.at(l1); 
        p1_double_prime = global_Rs_2.at(l2) * dfp_SamplePoints.at(p1 - global_num_points) 
                                + global_Ts_2.at(l2); 
    }
    if( p2 < global_num_points )
    {
        p2_prime = global_Rs.at(l1) * rp_SamplePoints.at(p2) + global_Ts.at(l1); 
        p2_double_prime = global_Rs.at(l2) * rp_SamplePoints.at(p2) + global_Ts.at(l2); 
    }
    else
    {
        p2_prime = global_Rs_2.at(l1) * dfp_SamplePoints.at(p2 - global_num_points) 
                        + global_Ts_2.at(l1); 
        p2_double_prime = global_Rs_2.at(l2) * dfp_SamplePoints.at(p2 - global_num_points) 
                                + global_Ts_2.at(l2); 
    }

    double d1 = (p1_prime - p1_double_prime).norm( );
    double d2 = (p2_prime - p2_double_prime).norm( );
    
    return static_cast<int>((d1 + d2)/2);

}

/******************************************************************************************/
void ssdr::set_nb_local(std::vector<std::vector<int>>& nb_local,
						std::vector<std::vector<double>>& nb_local_d,
						const std::vector<Point>& end_points,
						const std::vector<int>& p_curves,
						const Eigen::MatrixXd& p_neighbor_curves,
						double max_diagnal )
{
	
	for( int i = 0 ; i < end_points.size() ; i++ )
	{
		Point p1 = end_points.at(i);
		int c1 = p_curves.at(i);
		
		for( int j = i + 1 ; j < end_points.size() ; j++ )
		{
			Point p2 = end_points.at(j);
			int c2 = p_curves.at(j);
			double d = ( p2 - p1 ).norm();
				
			if( d < 0.2 * max_diagnal && 
			    p_neighbor_curves(c1,c2) != 1 )
			{
						nb_local.at(i).push_back(j);
						nb_local.at(j).push_back(i);
						nb_local_d.at(i).push_back(d);
						nb_local_d.at(j).push_back(d);
			}
			j++;
		}			
		i++;
		}				
}

/******************************************************************************************/

void ssdr::compute_neighborhood_info_2( const std::vector<std::vector<int>>& nbh_inf,
									  const std::vector<std::vector<double>>& nbh_d,
									  int **numNeighbors,
									  int ***neighborsIndexes,
									  int*** neighborsWeights,
									  double max_diagnal)
{
	 int num_points = nbh_inf.size();
	 int **nI = new int*[ num_points ];
	 int **nW = new int*[ num_points ];
	 double slope = -1 * ( 1 / ( max_diagnal * 0.2 ) );
	 
	 for( int i = 0 ; i < num_points; i++ )
     {
         std::cout << " i = " << i << "   " ;
     	nI[i] = new int[ (*numNeighbors)[i] + (nbh_inf.at(i)).size() ];
        nW[i] = new int[ (*numNeighbors)[i] + (nbh_inf.at(i)).size() ];
        
     	for( int j = 0 ; j < ((*numNeighbors)[i]) ; j++ )
     	{
     		nI[i][j] = (*neighborsIndexes)[i][j];
     		nW[i][j] = (*neighborsWeights)[i][j];
     	}
     	for( int j = (*numNeighbors)[i] ; j < ((*numNeighbors)[i] + (nbh_inf.at(i)).size()) ; j++ )
     	{
            int n = (*numNeighbors)[i];
     		nI[i][j] = nbh_inf.at(i).at(j - n);
     		double dist = nbh_d.at(i).at(j - n);
     		double r = slope * dist + 1 ;
     		nW[i][j] = r > 0 ? (static_cast<int>(r)) : 0;
     	}
     }
     
	 // add the quantity of the new neighbors to the list
	 for( int i = 0; i < num_points; i++ )
     {
         (*numNeighbors)[i] += nbh_inf.at(i).size();
     }
     
	 
     *neighborsIndexes = nI;
     *neighborsWeights = nW;
}

/*******************************************************************************************************/
void ssdr::combine_neighborhood_info_2( int** numNeighbors, int*** neighborsIndexes, int*** neighborsWeights,
                                      int* numNeighbors_dfp, int** neighborsIndexes_dfp, int** neighborsWeights_dfp,
                                      const std::vector<std::vector<int>>& nbh_temp, int num_points )
{
    int *ret_numNeighbors = new int[ 2 * num_points ];
    int **ret_neighborsIndexes = new int*[ 2 * num_points ];
    int **ret_neighborsWeights = new int*[ 2 * num_points ];
    std::vector< std::vector<int> > temp( 2 * num_points );
    std::vector< std::vector<int> > temp_2( 2 * num_points );

    // combine number of neighbors
    for( int i = 0 ; i < num_points ; i++ )
    {
        // initially we set the number of neighbors of each vertex 
        // this does not include the new bidirectional neighbors 
        ret_numNeighbors[i] = (*numNeighbors)[i];
        ret_numNeighbors[ i + num_points ] = (numNeighbors_dfp)[i]; 
        
        // first push_back the neighbors in each direction separately
        for( int k = 0 ; k < ret_numNeighbors[i] ; k++ )
        {
            temp[i].push_back( (*neighborsIndexes)[i][k] );
            temp_2[i].push_back( (*neighborsWeights)[i][k] );
        }
        
        for( int k = 0 ; k < ret_numNeighbors[ i + num_points ] ; k++ )
        {
            temp[ i + num_points ].push_back( (neighborsIndexes_dfp)[i][k] + num_points );
            temp_2[ i + num_points ].push_back( (neighborsWeights_dfp)[i][k] + num_points );
        }

        for( int k : nbh_temp.at(i) )
        {
            temp[i].push_back( k );
            temp_2[i].push_back( 1 );
        }
        
        for( int k : nbh_temp.at( i + num_points ) )
        {
            temp[ i + num_points ].push_back( k );
            temp_2[ i + num_points ].push_back( 1 );
        }

        std::sort( temp[ i + num_points ].begin(), temp[ i + num_points ].end() );
        temp[ i + num_points ].erase( unique( temp[ i + num_points ].begin(), temp[ i + num_points ].end() ), temp[ i + num_points ].end() );
        
        std::sort( temp[i].begin(), temp[i].end() );
        temp[i].erase( unique( temp[i].begin(), temp[i].end() ), temp[i].end() );
        
    }

    //combine neighbor Indexes and weights
    for( int i = 0 ; i < 2 * num_points ; i++ )
    {
        ret_neighborsIndexes[i] = new int[ (temp.at(i)).size() ];
        ret_neighborsWeights[i] = new int[ (temp.at(i)).size() ];

        ret_numNeighbors[i] = (temp.at(i)).size();

        for( int j = 0 ; j < ret_numNeighbors[i] ; j++ )
        {
            ret_neighborsIndexes[i][j] = (temp.at(i)).at(j);
            ret_neighborsWeights[i][j] = (temp_2.at(i)).at(j);
            for( int k = 0 ; k < (nbh_temp.at(i)).size() ; k++ )
            {
                if( ret_neighborsIndexes[i][j] == (nbh_temp.at(i)).at(k) )
                {
                    ret_neighborsWeights[i][j] = 1;
                }
            }
        }

    }
    
    *numNeighbors = ret_numNeighbors;
    *neighborsIndexes = ret_neighborsIndexes;
    *neighborsWeights = ret_neighborsWeights;

}

/*******************************************************************************************************/
/*******************************************************************************************************/

void ssdr::output_svg_2( const std::vector<Point>& pose_1_samples, 
                         const std::vector<Point>& pose_2_samples,
                         const std::vector<Point>& pose_2_samples_B, 
                         char const *file_name,
                         int index )
{
    std::vector<Eigen::Vector3d> colors;
    Eigen::Vector3d temp = Eigen::Vector3d( 0.5, 0.5, 0.5 );
    colors.push_back( temp );
    temp = Eigen::Vector3d( 1, 0, 0 );
    colors.push_back( temp );
    //colors.push_back( Eigen::Vector3d( 0.5, 0 , 0) );
    //colors.push_back( Eigen::Vector3d( 1, 1, 0 ) );
    colors.push_back( Eigen::Vector3d( 0.5, 0.5, 0 ) );
    colors.push_back( Eigen::Vector3d( 0, 1, 0 ) );
    colors.push_back( Eigen::Vector3d( 0, 0.5, 0 ) );
    colors.push_back( Eigen::Vector3d( 0, 1, 1 ) );
    colors.push_back( Eigen::Vector3d( 0, 0, 1 ) );
    colors.push_back( Eigen::Vector3d( 1, 0, 1 ) );
    colors.push_back( Eigen::Vector3d( 0.5, 0, 0.5 ) );
    
    cairo_t *cr;
    cairo_surface_t *surface;
    cairo_pattern_t *pattern;
    int x,y;
    char buffer[24]; // make sure it's big enough
    snprintf(buffer, sizeof(buffer), "%s_%d.svg", file_name, index);
    
    surface =
      (cairo_surface_t *)cairo_svg_surface_create(buffer, 1200.0, 1200.0);
    cr = cairo_create(surface);
    
    // draw the original rest pose image curves in black ( 0 , 0 , 0 )
    
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width( cr, 5 );

    Point p0, p1, p2, p3;
    
    for( auto & curve : rp_Curves )
    {
        //end points
        p0 = rp_CurveEndPoints.at(curve[0]);
        p3 = rp_CurveEndPoints.at(curve[3]);
        
        //tangent control points
        p1 = rp_CurveMiddlePoints.at(curve[1]);
        p2 = rp_CurveMiddlePoints.at(curve[2]);
        
        cairo_move_to( cr, p0(0), p0(1) );
        cairo_curve_to( cr, p1(0), p1(1), p2(0), p2(1), p3(0), p3(1) );
        cairo_stroke( cr );

    }

    // draw the original deformed pose image curves in navy ( 0 , 0 , 0.5 )
    
    cairo_set_source_rgb(cr, 0, 0, 0.5);
    cairo_set_line_width( cr, 5 );
    
    for( auto & curve : dfp_Curves )
    {
        //end points
        p0 = dfp_CurveEndPoints.at(curve[0]);
        p3 = dfp_CurveEndPoints.at(curve[3]);
        
        //tangent control points
        p1 = dfp_CurveMiddlePoints.at(curve[1]);
        p2 = dfp_CurveMiddlePoints.at(curve[2]);
        
        cairo_move_to( cr, p0(0), p0(1) );
        cairo_curve_to( cr, p1(0), p1(1), p2(0), p2(1), p3(0), p3(1) );
        cairo_stroke( cr );

    }
    
    //draw the deformed pose samples
    cairo_set_source_rgb(cr, 1, 0, 1);
    cairo_set_line_width( cr, 3 );
    for( auto& p : pose_2_samples_B )
    {
        cairo_arc(cr, p(0), p(1), 10, 0, 2*M_PI);
        cairo_stroke( cr );
    }

    //draw the rest pose samples
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width( cr, 3 );
    for( auto& p : pose_1_samples )
    {
        cairo_arc(cr, p(0), p(1), 5, 0, 2*M_PI);
        cairo_stroke( cr );
    }
    
    //draw the rest pose samples
    cairo_set_source_rgb(cr, 0, 0, 0.5);
    cairo_set_line_width( cr, 3 );
    for( auto& p : pose_2_samples )
    {
        cairo_arc(cr, p(0), p(1), 5, 0, 2*M_PI);
        cairo_stroke( cr );
    }
    cairo_destroy (cr);
    cairo_surface_destroy (surface);

}

/*********************************************************************************************/
}
