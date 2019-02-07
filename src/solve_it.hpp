#ifndef SOLVE_IT_HPP
#define SOLVE_IT_HPP

#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>

typedef Eigen::Vector2d Point;

Eigen::MatrixXd Calculate_A_i( const std::vector<Eigen::Matrix2d>& Rs, 
                               const std::vector<Eigen::Vector2d>& Ts, 
                               const Point& p_i );
Eigen::MatrixXd Calculate_P( const Eigen::MatrixXd& A );
Eigen::MatrixXd Calculate_q( const Eigen::MatrixXd& A, const Eigen::Vector2d v_i );
Eigen::MatrixXd Calculate_l( int num_labels );
Eigen::MatrixXd Calculate_u( int num_labels );
Eigen::MatrixXd Calculate_B( int num_labels );
//std::vectorXd QPSolve_i( const std::vector<Eigen::Matrix2d>& Rs, 
void QPSolve_i( const std::vector<Eigen::Matrix2d>& Rs, 
                         const std::vector<Eigen::Vector2d>& Ts, 
                         const Point& p_i, const Point& v_i );
#endif 
