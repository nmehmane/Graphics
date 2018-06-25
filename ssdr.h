
class ssdr {
  public:
    ssdr(Eigen::MatrixXd rest_pose);
    ~ssdr();
    add_pose(Eigen::MatrixXd new_pose);
  private:  //make all fields public
    Eigen::MatrixXd rest_pose;
    std::vector<Eigen::MatrixXd> frame_poses;
    std::vector<Eigen::Matrix3d> handle_Transforms;
    Eigen::MatrixXd Weight;
    int num_handles;
}

