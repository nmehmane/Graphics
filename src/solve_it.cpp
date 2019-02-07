#include <vector>
#include <algorithm>
#include "solve_it.hpp"
#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>


Eigen::MatrixXd Calculate_A_i( const std::vector<Eigen::Matrix2d>& Rs, 
                               const std::vector<Eigen::Vector2d>& Ts, 
                               const Point& p_i )
{
    // each R is made into a column of the matrix by doing the following 
    // R = [ a  b ]                   T
    //     [ c  d ]  ===>  [ a b c d ]  
    //
    // Rs = [ a1  a2  ... ]
    //      [ b1  b2  ... ]
    //      [ c1  c2  ... ]
    //      [ d1  d2  ... ]
    
    // declare a 4 x #labels matrix
    
    int n = Rs.size();

    Eigen::MatrixXd Rs_m = Eigen::MatrixXd::Zero( 4, n );
    
    Eigen::MatrixXd Ts_m = Eigen::MatrixXd::Zero( 2, n );
    
    for( int i = 0 ; i < n ; i++ )
    {
        const Eigen::Matrix2d& R = Rs.at(i);
        const Eigen::Vector2d& T = Ts.at(i);
        
        Rs_m.col(i) << R(0,0) , R(0,1) , R(1,0) , R(1,1);
        Ts_m.col(i) << T(0) , T(1);

    }
    
    // kronecker product
    // [ x  y  0  0 ]
    // [ 0  0  x  y ]
    Eigen::MatrixXd kp_p_i = Eigen::MatrixXd::Zero( 2, 4 );
    kp_p_i << p_i(0) , p_i(1) , 0 , 0,
              0 , 0 , p_i(0) , p_i(1); 
    
    return ( kp_p_i * Rs_m + Ts_m );
} 

// the OSQP solves convex quadratic programs (QPs) of the form
//        T           T
// 1/2 * W * P * W + q * W 
// subject to 
// l <= BX <= u

Eigen::MatrixXd Calculate_P( const Eigen::MatrixXd& A )
{
    return ( A.transpose() * A * 2 );
}

Eigen::MatrixXd Calculate_q( const Eigen::MatrixXd& A, const Eigen::Vector2d v_i )
{
    return ( -2 * v_i.transpose() * A ).transpose();
}

// subject to 
// l <= BX <= u
// we have 2 conditions that need to be combined
// w1 + w2 + ... = 1
// 0 <= wj <= 1
//
// so we will combine them as follows
// [ 0 ]    [ 1  0   0   ...  0 ]        [ 1 ]
// [ 0 ]    [ 0  1   0   ...  0 ]        [ 1 ] 
// [ . ]    [ ... Identity ...  ]        [ 1 ]
// [ . ] <= [ ................  ] W  < = [ 1 ]
// [ 0 ]    [ 0  0   0   ...  1 ]        [ 1 ]
// [ 1 ]    [ 1  1   1   ...  1 ]        [ 1 ]

Eigen::MatrixXd Calculate_l( int num_labels )
{
    Eigen::MatrixXd l = Eigen::MatrixXd::Zero( num_labels + 1 , 1 );
    l(num_labels,0) = 1.0;
    return l;
}

Eigen::MatrixXd Calculate_u( int num_labels )
{
    Eigen::MatrixXd u = Eigen::MatrixXd::Ones( num_labels + 1 , 1 );
    return u;
}

Eigen::MatrixXd Calculate_B( int num_labels )
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Ones( num_labels + 1 , num_labels );
    B.topLeftCorner(num_labels,num_labels) = Eigen::MatrixXd::Identity(num_labels,num_labels);
    /*for( int i = 0 ; i < num_labels ; i++ )
    {
        B(num_labels, i ) = 1.0;
    }
    */
    return B;
}

//std::vectorXd QPSolve_i( const std::vector<Eigen::Matrix2d>& Rs, 
void QPSolve_i( const std::vector<Eigen::Matrix2d>& Rs, 
                         const std::vector<Eigen::Vector2d>& Ts, 
                         const Point& p_i, const Point& v_i )
{
    int num_labels = Rs.size();
    Eigen::MatrixXd l = Calculate_l( num_labels );
    Eigen::MatrixXd u = Calculate_u( num_labels );
    Eigen::MatrixXd B = Calculate_B( num_labels );


    Eigen::MatrixXd A = Calculate_A_i( Rs, Ts, p_i );
    Eigen::MatrixXd P = Calculate_P( A );
    Eigen::MatrixXd q = Calculate_q( A, v_i );

    std::cout << "\nSolver Stuff\nl = \n" << l << "\nu = \n" << u << "\nB = \n" << B << std::endl;
    std::cout << "\nA = \n" << A << "\nP = \n" << P << "\nq = \n" << q << std::endl; 
    
   /* int r = P.rows();
    int c = P.cols();

    c_float* P_x = (c_float*)c_malloc(sizeof(c_float) * r * c);
    for( int i = 0 ; i < r ; i++ )
    {
        for( int j = 0 ; j < c ; j++ )
        {
            P_x[i*c + j] = P(i,j);
        }
    }
    c_int P_nnz = r * c;
    */
}

