#include "ssdr.h"
#include <assert.h>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <gsl/gsl_linalg.h> // for svd
#include <Eigen/Dense>
#include <iostream>

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
    assert( num_handles > 1 );
    
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
        double best_cluster = 0;

        // empty the clusters
        for(auto &w : clusters)
            w.clear();

        for(int i = 0 ; i < num_vertices ; i++)
        {

            best_err = 1000000;
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
                if( err < best_err )
                {
                    best_err = err;
                    best_cluster = j;
                }


            }
            clusters[best_cluster].push_back( i );
            
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
    //w->resize(nV, nH);
	//nH2=min(SPARSENESS, nH);
	//boneID.resize(nV);
	//for (int i=0; i<nV; i++) boneID[i].resize(nH2); */
/*	am.resize(nFr*3, nH2);
	bm.resize(nFr*3);
	xm.resize(nH2);
	aW.resize(nFr*3, nH);
	bW.resize(nFr*3);
*/
	//v2.resize(nV*3);
	//transCache.resize(nH);
	//for (int i=0; i<nH; i++) {
	//	transCache[i].resize(nV);
	//	for (int j=0; j<nV; j++) transCache[i][j].resize(3);
    //	}
    
    // initializes each handle with a vertex. it chooses the vertices that 
    // have the largest minimum distances to other vertices. after this we will have H 
    // cluster centers/handles initialized with vertices. 
/*	std::vector<int> sel(num_handles);
	sel[0] = rand() % num_vertices;
	for (int j=1; j<num_handels; j++)  
    {
		double gMax=0.0;
		for (int i=0; i<nV; i++)  
        {
			double gMin=1e+100;
			for (int j2=0; j2<j; j2++)  
            {
				double dm=0.0;
                double tmp = rest_pose(i,0) - rest_pose(sel[j2],0);
                dm += tmp*tmp;
                tmp = rest_pose(i,1) - rest_pose(sel[j2],1);
                dm += tmp*tmp;
				gMin=min(dm, gMin);
			}
			if (gMin>gMax)  
            {
				gMax=gMin;
				sel[j]=i;
			}
		}
	}

	//Matrix identity=Identity(3);
	//hr->resize(nFr);
	//for (int t=0; t<num_frames; t++) 
    for(int t = 0 ; t < num_handles ; t++)
    {
        // initialize the rotation matrices to the identiity matrix
        handle_rotations.push_back( Eigen::Matrix2d.identity() );
        //hr->at(t)=vector<Matrix>(nH, identity);
    }

	//ht->resize(nFr, nH*3);

    //initialize the translation matrixes
	vector<int> nei;
	for (int j=0; j<num_handles; j++) 
    {
		find_neighbors(sel[j], 20, nei);
		initTrans(j, nei);
	}
*/
	/*if (UPDATE_H==0) {
		aH.resize(nV, nH*4);
		bH.resize(nV);
	}
	
	if (UPDATE_W==2) initQRSum1(nH);
    */
//}
/*
// computes the average distance between 2 vertices across poses

double ssdr::avg_dis(int i, int j) 
{
	double res=0.0;
	double d, tmp;
	for (int t=0; t<num_frames; t++) {
		d=0;
        tmp = (frame_poses[t])(i,0) - (frame_poses[t])(j,0);
        d+=tmp*tmp;
        tmp = (frame_poses[t])(i,1) - (frame_poses[t])(j,1);
        d+=tmp*tmp;
		//for (int k=0; k<3; k++) {
		//	tmp=v->at(t)[i*3+k]-v->at(t)[j*3+k];
		//	d+=tmp*tmp;
		//}
		res+=std::sqrt(d);
	}
	return res/num_frames;
}

// puts the maxSelected closest neighbors in the selected container
void ssdr::find_neighbors(int sj, int maxSelect, std::vector<int>& selected) 
{
	
    //clears the container
    selected.clear();
	std::vector<double> dis;
	for (int i = 0; i < num_vertices; i++) 
    {
        
		if ((int)selected.size()<=maxSelect) 
        {
			selected.push_back(0);
			dis.push_back(0);
		}
        
		selected.back()=i;
        //  avg_dis gets the average distance between this vertex and the handle inputed by the user
        // across Frames

		dis.back()=avgDis(i, sj);
		
        int p=(int)selected.size()-1;
        // sorts the vertices according to their average distance to the handle
		while ((p>0)&&(dis[p]<dis[p-1])) {
			std::swap(dis[p], dis[p-1]);
			std::swap(selected[p], selected[p-1]);
			p--;
		}
	}

	if ((int)selected.size()>maxSelect) {
		selected.pop_back();
		dis.pop_back();
	}
}

    
void init_trans(int j, vector<int>& selected) 
{
	int ns=(int)selected.size();
	
	vector<double> p, q;
	Eigen::Matrix3d pqT, mu, ms, mv;
	//pqT.resize(3, 3);

	for (int t=0; t<num_frames; t++) 
    {
		p = q = vector<double>(3, 0.0);
        //for all selected neighbors of the vertex
		for (int id=0; id<ns; id++) {
            //gets the id th selected neighbor vertex number
			int i = selected[id];
            // the sum of coordinated of all neighboring verteces in that frame is computed
            // the sum of  coordinates of all neighboring verteces  in the rest pose
			for (int k=0; k<3; k++) {
				p[k]+=v->at(t)[i*3+k];
				q[k]+=u->at(i*3+k);
			}
		}
        // then sum( xt-0 xt-1 ... xt-v-1 ) sum( yt-0 yt-1 ... yt-v-1 ) is averaged over number of selected vertices
		for (int k=0; k<3; k++) {
			p[k]/=ns;
			q[k]/=ns;
		}

		pqT=Zeros(3);
		for (int id=0; id<ns; id++) {
			int i=selected[id];
			for (int k=0; k<3; k++) 
				for (int l=0; l<3; l++)
					pqT[l][k]+=(v->at(t)[i*3+k]-p[k])*(u->at(i*3+l)-q[l]);
		}

		if ((UPDATE_H>=2)&&svd3x3(pqT, mu, ms, mv)) {
			if (det3x3(pqT)<0)	
				for (int k=0; k<3; k++) mv[k][2]=-mv[k][2];
			mul(mv, mu.transpose(), hr->at(t)[j]);
		} else hr->at(t)[j]=Identity(3);
		for (int k=0; k<3; k++) {
			ht->at(t)[j*3+k]=p[k];
			for (int l=0; l<3; l++) ht->at(t)[j*3+k]-=hr->at(t)[j][k][l]*q[l];
		}
	}
}
*/
