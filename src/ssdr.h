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
    void match_samples_to_curve(void);
    void normalize_tgs(void);
    void find_correspondence(void);
    void set_neighborhood(int *numNeighbors,int **a,int** neighborsWeights);

    //double avg_dis(int i, int j);
    //void find_neighbors(int sj, int maxSelect, std::vector<int>& selected);
  //make all fields public
    Eigen::MatrixXd rest_pose;
    Eigen::MatrixXd rp_tangents;
    std::vector<Eigen::Vector4i> rp_curves;

    std::vector<Eigen::MatrixXd> frame_poses;
    Eigen::MatrixXd df_tangents;
    std::vector<Eigen::Vector4i> df_curves;

    std::vector<Eigen::Vector2d> handle_translations;
    std::vector<Eigen::Matrix2d> handle_rotations;
    Eigen::MatrixXd Weight;
    int num_handles;

    int num_frames;
    int num_vertices;

    int max_non_zero_weights;
    Eigen::MatrixXd weight_map;

    Eigen::MatrixXd rest_pose_samples;
    Eigen::MatrixXd rp_tg_samples;
    int* rp_sample_curves;

    Eigen::MatrixXd deformed_pose_samples;
    Eigen::MatrixXd dp_tg_samples;

    // 2d structures for storing each pair's R and T 
    std::vector<std::vector<Eigen::Matrix2d>> pair_R;
    std::vector<std::vector<Eigen::Vector2d>> pair_T;


};

#endif
