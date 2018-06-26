#include "ssdr.h"

ssdr::ssdr()
{
}

ssdr::~ssdr()
{
}
/*
void ssdr::init_bone_transforms(void)
{
    //sets the number of deformed frames/poses provided
    num_frames = frame_poses.size();
    // sets the number of vertices
    num_vertices = rest_pose.rows();
    
    // initialize the Weight matrix (V x H)
    Weight = Eigen::MatrixXd::Zero(num_vertices ,num_handles);

    // Sparsness constraint / maximum number of allowed non-zero weights for a vertex
    ( num_handles > SPARSENESS ) ? ( max_non_zero_weights = SPARSENESS) :
                                   ( max_non_zero_weights = num_handles );
    
    weight_map = Eigen::MatrixXd::Zero(num_vertices , max_non_zero_weights);

    //


    
}
*/
