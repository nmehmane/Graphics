#include "ssdr.h"
#include <assert.h>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <gsl/gsl_linalg.h> // for svd
#include <Eigen/Dense>
#include <iostream>
#include "gco/GCoptimization.h"

ssdr::ssdr()
{
}

ssdr::~ssdr()
{
}

/* Given:
        rest_pose: an N-by-2 array; each row is a position of a vertex
        deformed_pose: an N-by-2 array; each row is a deformed position of a vertex
        H: A positive integer specifying the number of handles.
    Returns:
        transformations: a length-H sequence of pairs: (2-by-2 rotation, 2d translation vector)
        weights: an N-by-H array where element i,j contains the weights for vertex i and transformation j
    
    Implements [Le and Deng 2012] as described in Section 3.1 of their paper.
    The `weights` matrix is binary.
    Each 2-by-2 rotation is guaranteed to actually be a rotation (no scaling or skew or reflection).
    
    Optional parameters:
        max_iterations: The maximum number of iterations
    '''
*/
void ssdr::init_bone_transforms(void)
{
    // assume that we have the same number of vertices in both poses
    assert( rest_pose.rows() == frame_poses[0].rows() );
    //assert( num_handles > 1 );
    
    //sets the number of deformed frames/poses provided
    num_frames = frame_poses.size();
    // sets the number of vertices
    num_vertices = rest_pose.rows();
    
    /*
     * 1 Generate H random clusters of vertices (sets of indices).
     * 2 For each cluster, solve for the transformation that minimizes
     *   the difference between the rest and deformed pose positions.
     * 3 For each vertex, find the transformation that minimizes
     *   its reconstruction error.
     * Repeat steps 2 and 3 (in the paper, five times).
     * 4 Turn the clusters into weights and return them.
     */
    int iter = 0;
    /*      STEP 1       */

    // create an array of vertex indices and shuffle     
    std::vector<int> vertex_indices;
    for( int i = 0 ; i < num_vertices ; i++ )
    {
        vertex_indices.push_back(i);
    }
    int dsds = vertex_indices.size();
    //printf( " meh = %d " , dsds );
    // obtain a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::shuffle (vertex_indices.begin(), vertex_indices.end(), std::default_random_engine(seed));

    // randomly assigns clusters
    std::vector<std::vector<int>> clusters;
    
    std::vector<int> temp;
    temp.push_back( vertex_indices[0] );

    std::vector<std::vector<int>> clusters_copy;
    printf(" num handles = %d \n\n", num_handles );
    for( int i = 1 ; i < num_vertices ; i++ )
    {
        if( ((i+1)/(num_vertices / num_handles) ==  num_vertices/(num_vertices / num_handles)) && i == (num_vertices - 1) ) 
        {
            temp.push_back( vertex_indices[i] );
            clusters.push_back( temp );
            break;
        }        
        else if( (i/(num_vertices / num_handles) !=  num_vertices/(num_vertices / num_handles)) && i%(num_vertices / num_handles) == 0 )
        {
            clusters.push_back( temp );
            temp.clear();
        }

        temp.push_back( vertex_indices[i] ); 
    }

    /*        STEP 2       */
    do{

        printf( " rows = %d col = %d " , (int)clusters.size(), (int)clusters[0].size() );
        printf(" the cluster is ----> \n" );
        for( auto &v : clusters )
        {
            for( auto &p : v )
            {
                printf("%d  ", p );
            }
            printf("\n");
        }
    
        int row = 0;
        clusters_copy.clear();
        handle_translations.clear();
        handle_rotations.clear();
        //for each cluster
        for( auto &cluster : clusters )
        {

            std::vector<int> temp2;
            for( auto p : cluster )
            {
                temp2.push_back(p);
            }
            clusters_copy.push_back( temp2 );
            if(cluster.size() == 0 )
            {
                continue;
            }
            // get the rest pose and deformed coordinates
            Eigen::MatrixXd rest_positions = Eigen::MatrixXd::Zero( (int)cluster.size(), 2);
            Eigen::MatrixXd deformed_positions = Eigen::MatrixXd::Zero( (int)cluster.size(), 2);
            row = 0;

            Eigen::Vector2d p_bar(0, 0);
            Eigen::Vector2d q_bar(0, 0);
            
            for( auto &index : cluster )
            {
                rest_positions(row,0) = rest_pose(index,0);
                rest_positions(row,1) = rest_pose(index,1);
                
                // rest center
                p_bar(0) += rest_positions(row,0);
                p_bar(1) += rest_positions(row,1);

                deformed_positions(row,0) = frame_poses[0](index,0);
                deformed_positions(row,1) = frame_poses[0](index,1);
                
                // deformed center
                q_bar(0) += deformed_positions(row,0);
                q_bar(1) += deformed_positions(row,1);

                row++;
            }
            //std::cout << " cluster rest positions matrix" << std::endl << rest_positions << std::endl; 
            //std::cout << " cluster deformed position matrix" << std::endl << deformed_positions << std::endl; 
            
            p_bar/=(int)cluster.size();
            q_bar/=(int)cluster.size();
            
            std::cout << " p_bar" << std::endl << p_bar << std::endl; 
            std::cout << " q_bar" << std::endl << q_bar << std::endl; 

            // X and Y are the d × n matrices that have xi and yi as their columns
            // xi :=pi −p_bar, yi :=qi −q_bar
            Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2,(int)cluster.size());
            Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(2,(int)cluster.size());
            
            double xi_x = 0;
            double xi_y = 0;
            double yi_x = 0;
            double yi_y = 0;

            for(int l = 0; l < (int)cluster.size(); l++ )
            {
                xi_x = rest_positions(l,0) - p_bar(0);
                xi_y = rest_positions(l,1) - p_bar(1);

                X.col(l) << xi_x, xi_y;
                
                yi_x = deformed_positions(l,0) - q_bar(0);
                yi_y = deformed_positions(l,1) - q_bar(1);

                Y.col(l) << yi_x, yi_y;
            }
            std::cout << " X" << std::endl << X << std::endl; 
            std::cout << " Y" << std::endl << Y << std::endl; 
            
            // 2x2 covariance matrix
            Eigen::MatrixXd S_tmp = X * Y.transpose();
            std::cout << " 2x2 covariance matrix" << std::endl << S_tmp << std::endl; 
            // Compute SVD
            Eigen::JacobiSVD<Eigen::MatrixXd> svd;
            svd.compute( S_tmp, Eigen::ComputeThinU | Eigen::ComputeThinV);
            std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
            std::cout << "Its left singular vectors are the columns of the thin U matrix:" 
                    << std::endl << svd.matrixU() << std::endl;
            std::cout << "Its right singular vectors are the columns of the thin V matrix:" 
                    << std::endl << svd.matrixV() << std::endl;
            
            // The optimal rotation is V diag( 1,1,det(V*U.traspose) ) U.transpose

            Eigen::MatrixXd Ut = svd.matrixU().transpose();
            Eigen::MatrixXd diag_mat = Eigen::MatrixXd::Identity(2,2);
            diag_mat(1,1) = (svd.matrixV() * Ut).determinant();
            Eigen::MatrixXd R = svd.matrixV() * diag_mat * Ut;
            handle_rotations.push_back( R );

            // The optimal translation is q_bar - R p_bar
            handle_translations.push_back( q_bar - (R*p_bar) );
            
            std::cout << " R" << std::endl << R << std::endl; 
            std::cout << " T" << std::endl << ( q_bar - (R*p_bar) ) << std::endl; 

        }
        
        /*        STEP 3        */

        // Keep the last cluster_indices around for a termination test. (clusters_copy)

        // For each vertex, find the best transformation
        Eigen::Vector2d r_p(0, 0);
        Eigen::Vector2d d_p(0, 0);
        Eigen::Vector2d rest_transformed;
        double err = 0;
        double best_err = 0;
        int best_cluster = 0;

        // empty the clusters
        for(auto &w : clusters)
            w.clear();

        for(int i = 0 ; i < num_vertices ; i++)
        {
            printf("Calculating the best clustering!\n");

            //best_err = 1000000;
            best_cluster = 0;
            err = 0;
            
            // rest pose of vertex i
            r_p(0) = rest_pose(i,0);
            r_p(1) = rest_pose(i,1);

            // deformed pose of vertex i
            d_p(0) = (frame_poses[0])(i,0);
            d_p(1) = (frame_poses[0])(i,1);

            for(int j = 0 ; j < num_handles ; j++)
            {
                rest_transformed = handle_rotations[j]*r_p + handle_translations[j]; 
                err = (rest_transformed - d_p).dot(rest_transformed - d_p);
                printf(" err = %lf for index %d handle %d \n", err, i, j);
                
                if( j == 0 )
                {
                    best_err = err;
                    best_cluster = j;
                }
                else if( err < best_err )
                {
                    best_err = err;
                    best_cluster = j;
                }


            }
            printf("best cluster is %d for index %d \n", best_cluster, i );
            clusters[best_cluster].push_back( i );
            
        }
        printf( " rows = %d col = %d " , (int)clusters.size(), (int)clusters[0].size() );
        printf(" the best cluster is ----> \n" );
        for( auto &v : clusters )
        {
            for( auto &p : v )
            {
                printf("%d  ", p );
            }
            printf("\n");
        }
        printf( " Finished iteration: %d\n", iter );
        iter++;
        
        bool terminate_loop = true;
        // check for termination : if the clustering hasn't changed since last iteration then terminate
        if( clusters.size() != clusters_copy.size() )
        { 
            terminate_loop = false;
        }

        for( int n = 0 ; n < num_handles ; n++ )
        {
            std::sort( clusters_copy[n].begin(), clusters_copy[n].end() );
            std::sort( clusters[n].begin(), clusters[n].end() );

            if( clusters[n].size() != clusters_copy[n].size() )
            {
                terminate_loop = false;
            }
            for( int m = 0 ; m < clusters[n].size(); m++ )
            {
                if( clusters[n][m] != clusters_copy[n][m] )
                {
                    terminate_loop = false;

                }
            }
        }

        if( terminate_loop == true )
        {
            printf("No change in clusters! Terminating.\n");
            goto stop;
        }
        else
        {
            printf("START : ITERATION %d \n",iter );
        }

    }while(iter <= MAX_ITERATIONS);
    
    stop:
    printf( " rows = %d col = %d " , (int)clusters.size(), (int)clusters[0].size() );
    printf(" final clustering  ----> \n" );
    for( auto &v : clusters )
    {
        for( auto &p : v )
        {
            printf("%d  ", p );
        }
        printf("\n");
    }
    /*      STEP 4      */

    Weight = Eigen::MatrixXd::Zero( num_vertices, num_handles );
    
    for( int nn = 0 ; nn < num_handles ; nn++ )
    {
        if( clusters[nn].size() != 0 )
        {
            for(auto w : clusters[nn])
            {
                Weight(w,nn) = 1;
            }
        }
    }

    std::cout << " Weights  " << std::endl << Weight << std::endl;

    
}

     /*   gsl_matrix* S = gsl_matrix_alloc( 2, 2 );
        gsl_matrix_set( S, 0, 0, S_tmp(0,0) );
        gsl_matrix_set( S, 0, 1, S_tmp(0,1) );
        gsl_matrix_set( S, 1, 0, S_tmp(1,0) );
        gsl_matrix_set( S, 1, 1, S_tmp(1,1) );
        
        gsl_matrix* V = gsl_matrix_alloc( 2, 2 );
        gsl_vector* sigma = gsl_vector_alloc(2); 
        gsl_vector* work = gsl_vector_alloc(2);
        
        //This function factorizes the M-by-N matrix A into the singular value decomposition 
        //A = U S V^T for M >= N. On output the matrix A is replaced by U. The diagonal elements
        // of the singular value matrix S are stored in the vector S. The singular values are 
        //non-negative and form a non-increasing sequence from S_1 to S_N. The matrix V contains
        // the elements of V in untransposed form. To form the product U S V^T it is necessary to
        // take the transpose of V. A workspace of length N is required in work.
        //This routine uses the Golub-Reinsch SVD algorithm.

        gsl_linalg_SV_decomp( S , V, sigma, work );

        handle_rotations.push_back( ); */




    // initialize the Weight matrix (V x H)
   // Weight = Eigen::MatrixXd::Zero(num_vertices ,num_handles);


bool Compare_vectors(Eigen::Vector2d a,Eigen::Vector2d b) { return (a(1) < b(1)); }

void ssdr::init_bone_transforms_2(void)
{
    // assume that we have the same number of vertices in both poses
    assert( rest_pose.rows() == frame_poses[0].rows() );
    //assert( num_handles > 1 );
    
    //sets the number of deformed frames/poses provided
    num_frames = frame_poses.size();
    // sets the number of vertices
    num_vertices = rest_pose.rows();
    
    /*
     * 1 Use le and Deng 2014 to initialize clusters:
     *      1) put all vertices into one cluster.
     *      2) compute the rest pose centroid of the cluster. here we have one cluster containing all vertices. ( O)
     *      For every vertex in the cluster:
	 *          3) we compute the reconstruction error of that vertex( call it e(s) )
	 *          4) compute the distance of the vertex’s rest pose position to the centroid of the rest pose ( call it d( us, O ) )
	 *          5) compute the product of e(s)* d(us, O ) 
	 *          6) check to see if the product of this vertex is larger than all others, if it is then make this the seed vertex s.
     *      now that all vertices in the cluster have been processed we have the cluster seed vertex s which has the maximum product.
     *
     *      for every vertex in the cluster:
	 *      7) computer the euclidean distance from each vertex’s rest pose to the seed vertex’s rest pose
     *
     *      8) use the euclidean distance to evenly split the cluster into two clusters.
     *
     * 2 For each cluster, solve for the transformation that minimizes
     *   the difference between the rest and deformed pose positions.
     * 3 For each vertex, find the transformation that minimizes
     *   its reconstruction error and reassign vertices to clusters.
     * Repeat step one exept split each of the clusters into two again and then do  steps 2 and 3 till we have the number of clusters
     * we need 
     * 4 Turn the clusters into weights and return them.
     */

    
    /*  1) Put all vertices into one cluster    */

    // create an array of vertex indices and shuffle     
    std::vector<int> vertex_indices;
    for( int i = 0 ; i < num_vertices ; i++ )
    {
        vertex_indices.push_back(i);
    }
    int dsds = vertex_indices.size();
    //printf( " meh = %d " , dsds );

    std::vector<std::vector<int>> clusters;
    
    std::vector<int> temp;
    for( int i = 0 ; i < num_vertices ; i++ )
    {
        temp.push_back( vertex_indices[i] );
    }
    clusters.push_back( temp );

    int iter = 1; // quantity of current clusters

    printf(" num handles = %d \n\n", num_handles );

    // used to hold the index of the vertex that is the clusters seed. the index of this array 
    // corresponds to the cluster number
    std::vector<int> cluster_seeds;  
    // the distance of the vertex’s rest pose position to the centroid of the rest pose ( call it d( us, O ) )
    std::vector<double> euclidean_distances; 
    for( int i = 0 ; i < num_vertices ; i++ )
    {
        euclidean_distances.push_back( 0.0 );
    }

    do{

        
        // Displays clusters
        printf( " rows = %d col = %d " , (int)clusters.size(), (int)clusters[0].size() );
        printf(" the cluster is ----> \n" );
        for( auto &v : clusters )
        {
            for( auto &p : v )
            {
                printf("%d  ", p );
            }
            printf("\n");
        }
    
        int row = 0;
        
        handle_translations.clear();
        handle_rotations.clear();
        //for each cluster
        for( auto &cluster : clusters )
        {
            /*  2) compute the rest pose centroid (p_bar) of the cluster(s).Initially we have one cluster
             *     containing all vertices.        */

            // get the rest pose and deformed coordinates
            Eigen::MatrixXd rest_positions = Eigen::MatrixXd::Zero( (int)cluster.size(), 2);
            Eigen::MatrixXd deformed_positions = Eigen::MatrixXd::Zero( (int)cluster.size(), 2);
            row = 0;

            Eigen::Vector2d p_bar(0, 0);
            Eigen::Vector2d q_bar(0, 0);
            
            for( auto &index : cluster )
            {
                rest_positions(row,0) = rest_pose(index,0);
                rest_positions(row,1) = rest_pose(index,1);
                
                // rest center
                p_bar(0) += rest_positions(row,0);
                p_bar(1) += rest_positions(row,1);

                deformed_positions(row,0) = frame_poses[0](index,0);
                deformed_positions(row,1) = frame_poses[0](index,1);
                
                // deformed center
                q_bar(0) += deformed_positions(row,0);
                q_bar(1) += deformed_positions(row,1);

                row++;
            }
            //std::cout << " cluster rest positions matrix" << std::endl << rest_positions << std::endl; 
            //std::cout << " cluster deformed position matrix" << std::endl << deformed_positions << std::endl; 
           
            // these are the centriods of the rest and deformed cluster 
            p_bar/=(int)cluster.size();
            q_bar/=(int)cluster.size();

            /*  3)using the centroid of the cluster ( in both rest and deformed pose ) we calculate the 
             *    optimal initial translation and rotation of the cluster   */
            
            std::cout << " p_bar" << std::endl << p_bar << std::endl; 
            std::cout << " q_bar" << std::endl << q_bar << std::endl; 

            // X and Y are the d × n matrices that have xi and yi as their columns
            // xi :=pi −p_bar, yi :=qi −q_bar
            Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2,(int)cluster.size());
            Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(2,(int)cluster.size());
            
            double xi_x = 0;
            double xi_y = 0;
            double yi_x = 0;
            double yi_y = 0;

            for(int l = 0; l < (int)cluster.size(); l++ )
            {
                xi_x = rest_positions(l,0) - p_bar(0);
                xi_y = rest_positions(l,1) - p_bar(1);
                
                // assigning all euclidean distances vertex index = index of the array
                X.col(l) << xi_x, xi_y;
               
                euclidean_distances[ cluster[l] ] = sqrt( (xi_x * xi_x) + ( xi_y * xi_y ) );
                printf( " euclidean_distances = %lf \n", euclidean_distances[ cluster[l] ]);

                yi_x = deformed_positions(l,0) - q_bar(0);
                yi_y = deformed_positions(l,1) - q_bar(1);

                Y.col(l) << yi_x, yi_y;
            }
            std::cout << " X" << std::endl << X << std::endl; 
            std::cout << " Y" << std::endl << Y << std::endl; 
            
            // 2x2 covariance matrix
            Eigen::MatrixXd S_tmp = X * Y.transpose();
            std::cout << " 2x2 covariance matrix" << std::endl << S_tmp << std::endl; 
            // Compute SVD
            Eigen::JacobiSVD<Eigen::MatrixXd> svd;
            svd.compute( S_tmp, Eigen::ComputeThinU | Eigen::ComputeThinV);
            std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
            std::cout << "Its left singular vectors are the columns of the thin U matrix:" 
                    << std::endl << svd.matrixU() << std::endl;
            std::cout << "Its right singular vectors are the columns of the thin V matrix:" 
                    << std::endl << svd.matrixV() << std::endl;
            
            // The optimal rotation is V diag( 1,1,det(V*U.traspose) ) U.transpose

            Eigen::MatrixXd Ut = svd.matrixU().transpose();
            Eigen::MatrixXd diag_mat = Eigen::MatrixXd::Identity(2,2);
            diag_mat(1,1) = (svd.matrixV() * Ut).determinant();
            Eigen::MatrixXd R = svd.matrixV() * diag_mat * Ut;
            handle_rotations.push_back( R );

            // The optimal translation is q_bar - R p_bar
            handle_translations.push_back( q_bar - (R*p_bar) );
            
            std::cout << " R" << std::endl << R << std::endl; 
            std::cout << " T" << std::endl << ( q_bar - (R*p_bar) ) << std::endl; 

        }

        /*  4) now that we have the R and T s for each cluster, we compute the reconstruction error 
         *     for each vertex   */
        
        // remove the previous cluster seeds
        cluster_seeds.clear();
        
        int cluster_i = 0;
        Eigen::Vector2d r_p(0, 0);
        Eigen::Vector2d d_p(0, 0);
        Eigen::Vector2d rest_transformed;
        double err = 0;

        double max_product = 0.0;
        int max_index = 0;

        for( auto &cluster : clusters )
        {

            for( auto &index : cluster )
            {
                // rest pose of vertex i
                r_p(0) = rest_pose(index,0);
                r_p(1) = rest_pose(index,1);

                // deformed pose of vertex i
                d_p(0) = (frame_poses[0])(index,0);
                d_p(1) = (frame_poses[0])(index,1);
                
                // reconstuction error for each vertex
                rest_transformed = handle_rotations[cluster_i]*r_p + handle_translations[cluster_i]; 
                err = (rest_transformed - d_p).dot(rest_transformed - d_p);
                printf(" err = %lf for index %d handle %d \n", err, index, cluster_i);

                /*  5) compute the distance of the vertex’s rest pose position to the centroid of the
                 *     rest pose ( call it d( us, O ) ). */

                 // this distance has already been computed for every vertex 
                
                /*  6) compute the product of e(s)* d(us, O )   */

                double product_e_d = err * euclidean_distances[index] ;
                
                /*  7) check to see if the product of this vertex is larger than all others, 
                 *     if it is then make this the seed vertex s.   */
                if( max_product < product_e_d )
                {
                    max_product = product_e_d;
                    max_index = index;
                }

            }
            cluster_i++;
            printf( " max index is = %d\n", max_index);
            // assign the seed of this cluster
            cluster_seeds.push_back( max_index );
        }
        
        // now we have the seeds we need
        /*  for every vertex in the cluster:
	     *     8) computer the euclidean distance from each vertex’s rest pose to the seed vertex’s rest pose */
        
        int cluster_number = 0;
        std::vector<Eigen::Vector2d> euc_dist;
        std::vector<std::vector<int>> clusters_cpy;
        //for each cluster
        for( auto &cluster : clusters )
        {
            for( auto &index : cluster )
            {   
                double x_temp = rest_pose(index,0) - rest_pose(cluster_seeds[cluster_number],0);
                double y_temp = rest_pose(index,1) - rest_pose(cluster_seeds[cluster_number],1);

                Eigen::Vector2d tt( index , sqrt( (x_temp*x_temp) + (y_temp*y_temp) ) );
                printf(" e-d = %lf \n", tt(1) );
                euc_dist.push_back( tt );

            }
            // sort the euclidean distances
            std::sort( euc_dist.begin(), euc_dist.end(), Compare_vectors);
            
            /*  9) use the euclidean distance to evenly split the cluster into two clusters.   */

            // split the cluster in two 
            std::vector<int> temp1;
            std::vector<int> temp2;
            for( int i = 0 ; i < cluster.size() ; i++ )
            {
                if( i < ( (int)cluster.size() / 2 ) )
                {
                    temp1.push_back( euc_dist[i][0] );
                }
                else
                {
                    temp2.push_back( euc_dist[i][0] );
                }
            }
            clusters_cpy.push_back( temp1 );
            clusters_cpy.push_back( temp2 );
            iter++;
            cluster_number++;
            euc_dist.clear();
        }
        printf(" the cluster_cpy is ----> \n" );
        for( auto &v : clusters_cpy )
        {
            for( auto &p : v )
            {
                printf("%d  ", p );
            }
            printf("\n");
        }
        
        /*  10) compute the rest pose centroid of the clusters  */

        handle_translations.clear();
        handle_rotations.clear();
        //for each cluster
        for( auto &cluster : clusters_cpy )
        {

            // get the rest pose and deformed coordinates
            Eigen::MatrixXd rest_positions = Eigen::MatrixXd::Zero( (int)cluster.size(), 2);
            Eigen::MatrixXd deformed_positions = Eigen::MatrixXd::Zero( (int)cluster.size(), 2);
            row = 0;

            Eigen::Vector2d p_bar(0, 0);
            Eigen::Vector2d q_bar(0, 0);
            
            for( auto &index : cluster )
            {
                rest_positions(row,0) = rest_pose(index,0);
                rest_positions(row,1) = rest_pose(index,1);
                
                // rest center
                p_bar(0) += rest_positions(row,0);
                p_bar(1) += rest_positions(row,1);

                deformed_positions(row,0) = frame_poses[0](index,0);
                deformed_positions(row,1) = frame_poses[0](index,1);
                
                // deformed center
                q_bar(0) += deformed_positions(row,0);
                q_bar(1) += deformed_positions(row,1);

                row++;
            }
            //std::cout << " cluster rest positions matrix" << std::endl << rest_positions << std::endl; 
            //std::cout << " cluster deformed position matrix" << std::endl << deformed_positions << std::endl; 
            
            p_bar/=(int)cluster.size();
            q_bar/=(int)cluster.size();
            
            std::cout << " p_bar" << std::endl << p_bar << std::endl; 
            std::cout << " q_bar" << std::endl << q_bar << std::endl; 

            // X and Y are the d × n matrices that have xi and yi as their columns
            // xi :=pi −p_bar, yi :=qi −q_bar
            Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2,(int)cluster.size());
            Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(2,(int)cluster.size());
            
            double xi_x = 0;
            double xi_y = 0;
            double yi_x = 0;
            double yi_y = 0;

            for(int l = 0; l < (int)cluster.size(); l++ )
            {
                xi_x = rest_positions(l,0) - p_bar(0);
                xi_y = rest_positions(l,1) - p_bar(1);

                X.col(l) << xi_x, xi_y;
                
                yi_x = deformed_positions(l,0) - q_bar(0);
                yi_y = deformed_positions(l,1) - q_bar(1);

                Y.col(l) << yi_x, yi_y;
            }
            std::cout << " X" << std::endl << X << std::endl; 
            std::cout << " Y" << std::endl << Y << std::endl; 
            
            // 2x2 covariance matrix
            Eigen::MatrixXd S_tmp = X * Y.transpose();
            std::cout << " 2x2 covariance matrix" << std::endl << S_tmp << std::endl; 
            // Compute SVD
            Eigen::JacobiSVD<Eigen::MatrixXd> svd;
            svd.compute( S_tmp, Eigen::ComputeThinU | Eigen::ComputeThinV);
            std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
            std::cout << "Its left singular vectors are the columns of the thin U matrix:" 
                    << std::endl << svd.matrixU() << std::endl;
            std::cout << "Its right singular vectors are the columns of the thin V matrix:" 
                    << std::endl << svd.matrixV() << std::endl;
            
            // The optimal rotation is V diag( 1,1,det(V*U.traspose) ) U.transpose

            Eigen::MatrixXd Ut = svd.matrixU().transpose();
            Eigen::MatrixXd diag_mat = Eigen::MatrixXd::Identity(2,2);
            diag_mat(1,1) = (svd.matrixV() * Ut).determinant();
            Eigen::MatrixXd R = svd.matrixV() * diag_mat * Ut;
            handle_rotations.push_back( R );

            // The optimal translation is q_bar - R p_bar
            handle_translations.push_back( q_bar - (R*p_bar) );
            
            std::cout << " R" << std::endl << R << std::endl; 
            std::cout << " T" << std::endl << ( q_bar - (R*p_bar) ) << std::endl; 

        }
        
        /*        STEP 3        */

        // For each vertex, find the best transformation
        double best_err = 0;
        int best_cluster = 0;

        // empty the clusters
        for(auto &w : clusters)
            w.clear();
        clusters.clear();
        for( auto &w :clusters_cpy )
        {
            std::vector<int> temp3;
            clusters.push_back( temp3 );
        }
        

        for(int i = 0 ; i < num_vertices ; i++)
        {
            printf("Calculating the best clustering!\n");

            //best_err = 1000000;
            best_cluster = 0;
            err = 0;
            
            // rest pose of vertex i
            r_p(0) = rest_pose(i,0);
            r_p(1) = rest_pose(i,1);

            // deformed pose of vertex i
            d_p(0) = (frame_poses[0])(i,0);
            d_p(1) = (frame_poses[0])(i,1);

            for(int j = 0 ; j < iter ; j++)
            {
                rest_transformed = handle_rotations[j]*r_p + handle_translations[j]; 
                err = (rest_transformed - d_p).dot(rest_transformed - d_p);
                printf(" err = %lf for index %d handle %d \n", err, i, j);
                
                if( j == 0 )
                {
                    best_err = err;
                    best_cluster = j;
                }
                else if( err < best_err )
                {
                    best_err = err;
                    best_cluster = j;
                }


            }
            printf("best cluster is %d for index %d \n", best_cluster, i );
            clusters[best_cluster].push_back( i );
            
        }
        printf( " rows = %d col = %d " , (int)clusters.size(), (int)clusters[0].size() );
        printf(" the best cluster is ----> \n" );
        for( auto &v : clusters )
        {
            for( auto &p : v )
            {
                printf("%d  ", p );
            }
            printf("\n");
        }
        printf( " Finished iteration: %d\n", iter );
        
    }while( iter != num_handles );
    
    
}

void ssdr::match_samples_to_curve(void)
{
    //currently we are selecting samples in such a way that each sample corresponds to the t=0 endpoint of one of
    // the curves. so we simply look though the t=0 endpoint of all curves to find out which curve each sample 
    // belongs to. 
    
    rp_sample_curves = new int[rest_pose_samples.rows()];

    rp_tg_samples = Eigen::MatrixXd::Zero(rest_pose_samples.rows(),2);
    for(int i = 0; i < rest_pose_samples.rows(); i++)
    {
        for(int j = 0; j < rp_curves.size(); j++)
        {
            int c_i = (rp_curves.at(j))(0);
            double err = rest_pose_samples(i,0) - rest_pose(c_i,0);
            err *= err;
            err += ( ( rest_pose_samples(i,1) - rest_pose(c_i,1) ) * ( rest_pose_samples(i,1) - rest_pose(c_i,1) ) );
            //printf("err = %lf\n", err);

            if( err < 0.5 )
            {
                rp_sample_curves[i] = j;
                int tg_i = (rp_curves.at(j))(1);
                // set the tangent
                rp_tg_samples(i,0) = rp_tangents(tg_i,0) - rest_pose(c_i,0);
                rp_tg_samples(i,1) = rp_tangents(tg_i,1) - rest_pose(c_i,1);

            }
        }
    }
    std::cout << " rest pose sample tg \n" << rp_tg_samples << std::endl;

    Eigen::MatrixXd deformed_pose = frame_poses.at(0);
    dp_tg_samples = Eigen::MatrixXd::Zero(deformed_pose_samples.rows(),2);

    for(int i = 0; i < deformed_pose_samples.rows(); i++)
    {
        //printf("\n***************\n");
        for(int j = 0; j < df_curves.size(); j++)
        {
            int c_i = (df_curves.at(j))(0);
            double err = deformed_pose_samples(i,0) - deformed_pose(c_i,0);
            err *= err;
            err += ( ( deformed_pose_samples(i,1) - deformed_pose(c_i,1) ) * ( deformed_pose_samples(i,1) - deformed_pose(c_i,1) ) );
            //printf("err = %lf\n", err);

            if( err < 0.5 )
            {
                int tg_i = (rp_curves.at(j))(1);
                // set the tangent
                dp_tg_samples(i,0) = df_tangents(tg_i,0) - deformed_pose(c_i,0);
                dp_tg_samples(i,1) = df_tangents(tg_i,1) - deformed_pose(c_i,1);
            }
        }
    }
    std::cout << " deformed pose sample tg \n" << dp_tg_samples << std::endl;

}

void ssdr::normalize_tgs(void)
{
    double length = 0.0;
    for(int i = 0; i < dp_tg_samples.rows(); i++)
    {
        length = (dp_tg_samples(i,0)*dp_tg_samples(i,0)) + (dp_tg_samples(i,1)*dp_tg_samples(i,1));
        length = std::sqrt( length );
        printf(" length = %lf\n", length );
        dp_tg_samples(i,0) = dp_tg_samples(i,0) / length;
        dp_tg_samples(i,1) = dp_tg_samples(i,1) / length;
        printf(" normal df x = %lf y = %lf\n", dp_tg_samples(i,0),dp_tg_samples(i,1));
    }
    
    for(int i = 0; i < rp_tg_samples.rows(); i++)
    {
        length = (rp_tg_samples(i,0)*rp_tg_samples(i,0)) + (rp_tg_samples(i,1)*rp_tg_samples(i,1));
        length = std::sqrt( length );
        rp_tg_samples(i,0) = rp_tg_samples(i,0) / length;
        rp_tg_samples(i,1) = rp_tg_samples(i,1) / length;
        printf(" normal rp x = %lf y = %lf\n", rp_tg_samples(i,0),rp_tg_samples(i,1));
    }
}

void ssdr::find_correspondence(void)
{
    // for each sample in the rest-pose we pair it up with each sample from the 
    // deformed pose
    Eigen::Vector2d p_bar(0, 0);
    Eigen::Vector2d q_bar(0, 0);

    // holds the coordinates of the second rest pose and the second deformed pose
    Eigen::Vector2d r_point_2(0, 0);
    Eigen::Vector2d d_point_2(0, 0);

    for( int i = 0; i <rest_pose_samples.rows(); i++)
    {
        std::vector<Eigen::Matrix2d> temp_R;
        std::vector<Eigen::Vector2d> temp_T;

        // Add the normalized tangent to the coordinates of the vertex to get a second point
        r_point_2(0) = rest_pose_samples(i,0) + rp_tg_samples(i,0);
        r_point_2(1) = rest_pose_samples(i,1) + rp_tg_samples(i,1);

        for(int j = 0; j < deformed_pose_samples.rows(); j++)
        {
            // Add the normalized tangent to the coordinates of the vertex to get a second point
            d_point_2(0) = deformed_pose_samples(j,0) + dp_tg_samples(j,0);
            d_point_2(1) = deformed_pose_samples(j,1) + dp_tg_samples(j,1);
            
            //now we have a cluster containing two points in each pose we find the center of each cluster

            p_bar(0) = (rest_pose_samples(i,0) + r_point_2(0))/2;
            p_bar(1) = (rest_pose_samples(i,1) + r_point_2(1))/2;
            
            q_bar(0) = (deformed_pose_samples(j,0) + d_point_2(0))/2;
            q_bar(1) = (deformed_pose_samples(j,1) + d_point_2(1))/2;

            std::cout << " p_bar" << std::endl << p_bar << std::endl;
            std::cout << " q_bar" << std::endl << q_bar << std::endl;

            // using the centroid of the cluster ( in both rest and deformed pose ) we calculate the 
            // optimal initial translation and rotation of the cluster

            // X and Y are the d × (n=2) matrices that have xi and yi as their columns
            // xi = pi −p_bar, yi = qi −q_bar
            Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2,2);
            Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(2,2);

            double xi_x = 0;
            double xi_y = 0;
            double yi_x = 0;
            double yi_y = 0;
            
            //setting the first columns
            xi_x = rest_pose_samples(i,0) - p_bar(0);
            xi_y = rest_pose_samples(i,1) - p_bar(1);
            
            X.col(0) << xi_x, xi_y;
            
            yi_x = deformed_pose_samples(j,0) - q_bar(0);
            yi_y = deformed_pose_samples(j,1) - q_bar(1);

            Y.col(0) << yi_x, yi_y;

            //setting the second columns
            xi_x = r_point_2(0) - p_bar(0);
            xi_y = r_point_2(1) - p_bar(1);
            
            X.col(1) << xi_x, xi_y;
            
            yi_x = d_point_2(0) - q_bar(0);
            yi_y = d_point_2(1) - q_bar(1);

            Y.col(1) << yi_x, yi_y;

            std::cout << " X" << std::endl << X << std::endl;
            std::cout << " Y" << std::endl << Y << std::endl;

            // 2x2 covariance matrix
            Eigen::MatrixXd S_tmp = X * Y.transpose();
            std::cout << " 2x2 covariance matrix" << std::endl << S_tmp << std::endl;
            
            // Compute SVD
            Eigen::JacobiSVD<Eigen::MatrixXd> svd;
            svd.compute( S_tmp, Eigen::ComputeThinU | Eigen::ComputeThinV);
            std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
            std::cout << "Its left singular vectors are the columns of the thin U matrix:"
                      << std::endl << svd.matrixU() << std::endl;
            std::cout << "Its right singular vectors are the columns of the thin V matrix:"
                      << std::endl << svd.matrixV() << std::endl;

            // The optimal rotation is V diag( 1,1,det(V*U.traspose) ) U.transpose

            Eigen::MatrixXd Ut = svd.matrixU().transpose();
            Eigen::MatrixXd diag_mat = Eigen::MatrixXd::Identity(2,2);
            diag_mat(1,1) = (svd.matrixV() * Ut).determinant();
            Eigen::MatrixXd R = svd.matrixV() * diag_mat * Ut;
            Eigen::Vector2d T = q_bar - (R*p_bar);

            // store the pair's R and T
            temp_R.push_back( R );
            temp_T.push_back( T );

            std::cout << " R" << std::endl << R << std::endl;
            std::cout << " T" << std::endl << T << std::endl;

        }

        pair_R.push_back(temp_R);
        pair_T.push_back(temp_T);
    }

    // now we have n x n candidate R,T s for each sample
    /* For each vertex A in the rest position:
            For every R&T:
                Apply The current R and T to vertex A and get a position C.
                Find the closest vertex in the deformed pose to this point
                Calculate the distance and store it a a goodness metric
                
                Apply R and T to the Tangent of this point in rest pose and see if it points to the 
                same direction as the tangent of the point in the deformed pose
    */
    //holds the transformed sample point from the rest pose
    Eigen::Vector2d transformed_r_p(0.0, 0.0);
    
    //stores the number of labels
    int num_labels = rest_pose_samples.rows() * rest_pose_samples.rows();
    //number of sample points
    int num_samples = rest_pose_samples.rows();

    // setting up the array for data costs i.e. the cost of using Label l_i ( Ri,Ti ) for point p
    int *data = new int[num_samples * num_labels];

    // setting up the array for smooth costs
    int *smooth = new int[num_labels*num_labels];

    for( int i = 0; i <rest_pose_samples.rows(); i++)
    {
        int label_index = 0;
            printf("data terms for %d\n" , i );

        //for each sample in the rest_pose apply all possible R,T pairs one by one to get the transformed
        // point
        for( int ii = 0; ii < rest_pose_samples.rows(); ii++)
        {
            //the sample point at the rest pose is stored as a 2D vector
            Eigen::Vector2d rp_sample_p = (rest_pose_samples.row(i)).transpose();

            for(int jj = 0; jj < deformed_pose_samples.rows(); jj++)
            {
                transformed_r_p = ( pair_R[ii][jj] * rp_sample_p ) + pair_T[ii][jj];
                //std::cout << " transformed rest sample point" << std::endl << transformed_r_p << std::endl;

                double min_label_cost = 0.0;

                for(int kk = 0; kk < deformed_pose_samples.rows(); kk++)
                {
                    //compute euclidean distance squared between the transformed rest point and the deformed point
                    double e_distance = ( transformed_r_p - (deformed_pose_samples.row(kk)).transpose() ).squaredNorm( );

                    if( kk == 0 )
                    {
                        min_label_cost = e_distance;
                    }
                    // choose the vertex that gives the minimum
                    else if( e_distance < min_label_cost )
                    {
                        min_label_cost = e_distance;
                    }
                }

                data[i*num_labels + (label_index) ] = (int)min_label_cost;
                printf(" %d  " , data[i*num_labels + (label_index) ] );

                label_index++;
            }
        }

            printf("\n");

    }
    
    for ( int l1 = 0; l1 < num_labels; l1++ )
    {
        printf("smoothness term for %d\n", l1);
        for (int l2 = 0; l2 < num_labels; l2++ )
        {
            smooth[l1+l2*num_labels] = (l1-l2)*(l1-l2) < 1  ? 0 : 10000;
            printf(" %d ",smooth[l1+l2*num_labels]);
        }
        printf("\n");
    }

	int *result = new int[num_samples];
    try{
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_samples,num_labels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);

		// now set up a neighborhood system
        int** neighborsIndexes = new int*[rest_pose_samples.rows()];
        int* numNeighbors = new int[rest_pose_samples.rows()];
        int** neighborsWeights = new int*[rest_pose_samples.rows()];

        set_neighborhood(numNeighbors,neighborsIndexes,neighborsWeights);

        //pass in all neighbor information at once
        gc->setAllNeighbors(numNeighbors,neighborsIndexes,neighborsWeights);

		printf("\nBefore optimization energy is %lld",gc->compute_energy());
		gc->expansion(2);// run expansion for 2 iterations.
        gc->swap(2);
		printf("\nAfter optimization energy is %lld",gc->compute_energy());

		for ( int  i = 0; i < num_samples; i++ )
        {
			result[i] = gc->whatLabel(i);
            printf(" for point %d label %d was chosen\n",i,result[i]);
        }

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

	delete [] result;
	delete [] smooth;
	delete [] data;

}


void ssdr::set_neighborhood(int *numNeighbors,int **a,int** neighborsWeights)
{

    for( int i = 0; i <rest_pose_samples.rows(); i++)
    {
        std::vector<int> temp;

        for( int j = i; j < rest_pose_samples.rows(); j++)
        {
            if( rp_sample_curves[i] == rp_sample_curves[j] )
            {
                //then j is the t=1 end point of the same curve
                temp.push_back(j);
            }
            else if( (rp_sample_curves[i]-1) == rp_sample_curves[j] )
            {
                //then j is the t=0 end point of the curve whose t=1
                // end point is rest_pose_samples.rows(i)
                temp.push_back(j);
            }
        }
        numNeighbors[i] = temp.size();
        a[i] = new int[temp.size()];
        neighborsWeights[i] = new int[temp.size()];

        
        for(int k = 0 ; k < temp.size(); k++)
        {
            a[i][k] = temp.at(k);
            neighborsWeights[i][k] = 1;
        }
    }

}
