#ifndef SSDR_H
#define SSDR_H

#include <Eigen/Core>
#include <vector>

#define SPARSENESS 4
#define MAX_ITERATIONS 100

class ssdr {
    
  public:
    ssdr();
    ~ssdr();
    void init_bone_transforms(void);
    void init_bone_transforms_2(void);
    //double avg_dis(int i, int j);
    //void find_neighbors(int sj, int maxSelect, std::vector<int>& selected);
  //make all fields public
    Eigen::MatrixXd rest_pose;
    std::vector<Eigen::MatrixXd> frame_poses;
    std::vector<Eigen::Vector2d> handle_translations;
    std::vector<Eigen::Matrix2d> handle_rotations;
    Eigen::MatrixXd Weight;
    int num_handles;

    int num_frames;
    int num_vertices;

    int max_non_zero_weights;
    Eigen::MatrixXd weight_map;


};

#endif
