#ifndef SSDR_H
#define SSDR_H

#include <Eigen/Core>
#include <vector>

#define SPARSENESS 4

class ssdr {
    
  public:
    ssdr();
    ~ssdr();
    void init_bone_transforms(void);
  //make all fields public
    Eigen::MatrixXd rest_pose;
    std::vector<Eigen::MatrixXd> frame_poses;
    std::vector<Eigen::Matrix3d> handle_Transforms;
    Eigen::MatrixXd Weight;
    int num_handles;

    int num_frames;
    int num_vertices;

    int max_non_zero_weights;
    Eigen::MatrixXd weight_map;


};

#endif
