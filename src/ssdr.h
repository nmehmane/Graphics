#ifndef SSDR_H
#define SSDR_H

#include <Eigen/Core>
#include <vector>

#define SPARSENESS 4
#define MAX_ITERATIONS 100

namespace ssdr_ns
{

typedef Eigen::Vector2d Point;

class ssdr {
     
  public:

    ssdr();
    ~ssdr();
    void perform_ssdr();
    std::vector<Point> sample_curves( const std::vector<Point>& end_points, 
                                      const std::vector<Point>& tangent_points, 
                                      const std::vector<Eigen::Vector4i>& cubic_curves,
                                      const int samples_per_curve,
                                      std::vector<Eigen::Vector2d> *p_samples_tg,
                                      std::vector<int>* p_curves);

    Point evaluate_curve( const Point& p0, const Point& p1,
                          const Point& p2, const Point& p3,
                          const double t );

    Point evaluate_curve_tg( const Point& p0, const Point& p1, 
                             const Point& p2, const Point& p3, 
                             const double t );

    std::vector<Point> normalize_tangents( const std::vector<Point>& sp_tg );
    
    std::vector<double> calculate_param(const std::vector<Point>& points,
                                        const std::vector<Point>& tg_points,
                                        const std::vector<Eigen::Vector4i>& cubic_curves,
                                        std::vector<int>* curve_numbers);

    std::vector<Point> compute_tangents( const std::vector<double>& points_param_t,
                                         const std::vector<int>& curve_numbers,
                                         const std::vector<Eigen::Vector4i>& cubic_curves,
                                         const std::vector<Point>& end_points,
                                         const std::vector<Point>& tangent_points);

    void compute_cluster_RT( std::vector<Point> pose_1,
                             std::vector<Point> pose_2,
                             Eigen::Matrix2d* R,
                             Eigen::Vector2d* T );

    
    void initialize_transformations( const std::vector<Point>& pose1_samples,
                                     const std::vector<Point>& pose2_samples,
                                     const std::vector<Point>& pose1_samples_tg,
                                     const std::vector<Point>& pose2_samples_tg,
                                     std::vector<Eigen::Matrix2d>* Rs,
                                     std::vector<Eigen::Vector2d>* Ts );

    int find_closest_point( const Point& p , const std::vector<Point>& points , double* dist );

    int* set_data_term( const std::vector<Point>& pose_1,
                              const std::vector<Point>& pose_2,
                              const std::vector<Eigen::Matrix2d>& Rs,
                              const std::vector<Eigen::Vector2d>& Ts );

    int* set_smoothness_term( int num_labels );

    Eigen::MatrixXd set_curve_neighborhood( const std::vector<Eigen::Vector4i>& pose_curves );

    void compute_neighborhood_info( const std::vector<int>& curve_numbers,
                                    const Eigen::MatrixXd& curve_nbh,
                                    int **numNeighbors, 
                                    int ***neighborsIndexes, 
                                    int*** neighborsWeights );

    std::vector<int> assign_labels( int num_points,
                                    int num_labels,
                                    int *data_terms,
                                    int *smoothness_term,
                                    int *numNeighbors, 
                                    int **neighborsIndexes, 
                                    int** neighborsWeights );


    void set_restpose_clusters( const std::vector<int>& label_assignment,
                                const int num_points,
                                std::vector<std::vector<int>>* pose_1_clusters,
                                std::vector<int>* label_numbers);

    void compute_corresponding_clusters( const std::vector<Point>& pose_1,
                                         const std::vector<Point>& pose_2, 
                                         const std::vector<int>& label_assignment,
                                         const std::vector<Eigen::Matrix2d>& Rs,
                                         const std::vector<Eigen::Vector2d>& Ts,
                                         std::vector<std::vector<int>>& pose_1_clusters,
                                         std::vector<std::vector<int>>* pose_2_clusters );

    void recompute_RT( const std::vector<Point>& pose_1,
                       const std::vector<Point>& pose_2, 
                       std::vector<Eigen::Matrix2d>* Rs,
                       std::vector<Eigen::Vector2d>* Ts,
                       const std::vector<std::vector<int>>& pose_1_clusters,
                       const std::vector<std::vector<int>>& pose_2_clusters );


    /*    Rest Pose    */
    // is set by svg parser to hold the end points of the curves
    std::vector<Point> rp_CurveEndPoints;
    // is set by svg parser to hold the tangent control points of the curves
    std::vector<Point> rp_CurveMiddlePoints;
    // is set by svg parser to hold the rest pose curves
    std::vector<Eigen::Vector4i> rp_Curves;
    
    /*    Deformed Pose    */
    // is set by svg parser to hold the end points of the curves
    std::vector<Point> dfp_CurveEndPoints;
    // is set by svg parser to hold the tangent control points of the curves
    std::vector<Point> dfp_CurveMiddlePoints;
    // is set by svg parser to hold the deformed pose curves
    std::vector<Eigen::Vector4i> dfp_Curves;
    
    /*    Rest Pose Samples    */
    std::vector<Point> rp_SamplePoints;
    std::vector<Point> rp_SamplePoint_tg;
    std::vector<int> rp_SamplePoint_Param_t;
    std::vector<int> rp_SamplePoint_Curves;

    /*    Deformed Pose Samples    */
    std::vector<Point> dfp_SamplePoints;
    std::vector<Point> dfp_SamplePoint_tg;
    std::vector<int> dfp_SamplePoint_Param_t;
    std::vector<int> dfp_SamplePoint_Curves;

    int num_handles;


};

}
#endif
