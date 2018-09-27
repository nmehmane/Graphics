#include "svg_parser.h"
#include "ssdr.h"

#include <iostream>
#include <sstream>
#include <Eigen/Core>
#include "rapidxml_ns/rapidxml_ns.hpp"
#include "rapidxml_ns/rapidxml_ns_utils.hpp"
#include "svgpp/policy/xml/rapidxml_ns.hpp"
#include "svgpp/svgpp.hpp"
#include "Graphics_Gems/NearestPoint.h"
#include "Graphics_Gems/GraphicsGems.h"

using namespace svgpp;


void Context::on_enter_element(tag::element::any)
{
      //std::cout<<" On_enter_element " << std::endl;
}

void Context::on_exit_element()
{
      //std::cout<<" On_exit_element " << std::endl;
}
  
// converts relative coordinates to absolute coordinates if (m)
void Context::path_move_to(double x, double y, tag::coordinate::absolute)
{
    // M --> absolute coordinates follow  m --> relative coordinates follow
    // If a relative moveto (m) appears as the first element of the path, then 
    // it is treated as a pair of absolute coordinates. In this case, subsequent
    // pairs of coordinates are treated as relative even though the initial 
    // moveto is interpreted as an absolute moveto.
    
    path_start = true;
    //std::cout<<" path_move_to " << std::endl;
    //printf(" x = %lf y = %lf \n",x,y);
    vertices.push_back(x);
    vertices.push_back(y);
    // index of the current vertex that was pushed back
    int i = (vertices.size()/2)-1;
    Eigen::Vector4i temp(i, 0, 0, 0);
    curves.push_back(temp);

    
}

void Context::path_line_to(double x, double y, tag::coordinate::absolute)
{
      //std::cout<<" path_line_to " << std::endl;
}

void Context::path_cubic_bezier_to(double x1, double y1,double x2, double y2,
                                   double x, double y, tag::coordinate::absolute)
{
    // converts relative coordinates to absolute coordinates if (c)
      //std::cout<<" path_cubic_bezier_to " << std::endl;
      //printf(" x = %lf y = %lf \n",x,y);
    vertices.push_back(x);
    vertices.push_back(y);
    
    //printf(" x = %lf y = %lf \n",x,y);
    //printf(" x1 = %lf y1 = %lf \n",x1,y1);
    //printf(" x2 = %lf y2 = %lf \n",x2,y2);
    // index of the current vertex that was pushed back
    int i = (vertices.size()/2)-1;

    if( path_start == false )
    {
        Eigen::Vector4i temp( i-1, 0, 0, 0 );
        curves.push_back(temp);

    }
    else
    {
        path_start = false;
    }
    // sets temp( i0 ,0 , 0 , i1=i )
    (curves.back())(3) = i; 
    
    //push back the tangents 
    v_tangents.push_back(x1);
    v_tangents.push_back(y1);
    v_tangents.push_back(x2);
    v_tangents.push_back(y2);

    i = (v_tangents.size()/2)-1;
    (curves.back())(2) = i;
    (curves.back())(1) = i-1;
    //std::cout << "curve = "<< curves.back() << std::endl;

    
}

void Context::path_quadratic_bezier_to(double x1, double y1, double x, double y,
                                       tag::coordinate::absolute)
{
     // std::cout<<" path_quadratic_bezier_to " << std::endl;
}

void Context::path_elliptical_arc_to(double rx, double ry, double x_axis_rotation,
                                     bool large_arc_flag, bool sweep_flag,
                                     double x, double y, tag::coordinate::absolute)
{

}

void Context::path_close_subpath()
{
    //std::cout<<" path_close_subpath " << std::endl;
}

void Context::path_exit()
{
    //std::cout<<" path_exit " << std::endl;
}

void loadSvg(xml_element_t xml_root_element, Eigen::MatrixXd& rest_pose, Eigen::MatrixXd& rp_tangents, std::vector<Eigen::Vector4i>& rp_curves)
{
    Context context;
    document_traversal<
      processed_elements<processed_elements_t>,
      processed_attributes<traits::shapes_attributes_by_element>
    >::load_document(xml_root_element, context);
    //  an N x 2 matrix
    rest_pose = Eigen::MatrixXd::Zero((context.vertices.size()/2) ,2);
    rp_tangents = Eigen::MatrixXd::Zero((context.v_tangents.size()/2) ,2);

    std::vector<double>::iterator it;
    int row = 0;
    for( it = context.vertices.begin(); it != context.vertices.end(); it++,row++  )
    {
        rest_pose(row,0) = *it;
        //printf(" index of vertex = %d\n", row );
        //printf(" x = %lf", rest_pose(row,0)); 
        rest_pose(row,1) = *(++it);
        //printf(" y = %lf\n", rest_pose(row,1)); 
    }

    row = 0;
    for( it = context.v_tangents.begin(); it != context.v_tangents.end(); it++,row++  )
    {
        rp_tangents(row,0) = *it;
        //printf(" index of tangent = %d\n", row );
        //printf(" x' = %lf",rp_tangents(row,0)); 
        rp_tangents(row,1) = *(++it);
        //printf(" y' = %lf\n", rp_tangents(row,1)); 
    }
    
    std::vector<Eigen::Vector4i>::iterator itt;
    row = 0;
    for( itt = context.curves.begin(); itt != context.curves.end(); itt++)
    {
        rp_curves.push_back(*itt);
        //printf("curve ( %d )\n", row );
        //std::cout<<  rp_curves.back() << std::endl;
        row++;
    }
}

int main(int argc, char** argv)
{
    // change args to ./exe <number of handles> <svg rest-pose> <number of frames> <svg pose 1> <svg pose2> ...
    // Minimum requirements: 1)number of handels 2)rest pose 3)number of frames ( minimum 2 ) 4,5)minimum of two frames
    if( argc < 7 )
    {
        std::cout<<"./exe <number of handles> <svg rest-pose> <number of frames> <svg pose 1> <svg pose2> ..."
                   " <rest-pose samples> <deformed-pose samples>"<<std::endl;
        return -1;
    }
    
    //create ssdr element 
    ssdr ssdr_elem;
    
    // set number of  handles
    std::istringstream s(argv[1]);
    if( !(s >> ssdr_elem.num_handles) )
    {
        std::cerr << "Invalid number of handles" << std::endl;
    }

    // Set the rest pose
    try
    {
        rapidxml_ns::file<> xml_file(argv[2]);  
        rapidxml_ns::xml_document<> doc;    // character type defaults to char
        doc.parse<rapidxml_ns::parse_no_string_terminators>(xml_file.data());  
        loadSvg(doc.first_node(), ssdr_elem.rest_pose, ssdr_elem.rp_tangents, ssdr_elem.rp_curves);

    }
    catch (std::exception const & e)
    {
        std::cerr << "Error loading rest Pose SVG: " << e.what() << std::endl;
        return 1;
    }

    // Set the other poses
    // need a minimum of two other poses
    std::istringstream ss( argv[3] );
    int num_poses;
    
    if(!(ss >> num_poses))
    {
        std::cerr << " Invalid number of poses " << std::endl;
        return 2;
    }
    if( (argc - 6) != num_poses )
    {
        std::cerr << " number of poses does not match number of files " << std::endl;
        return 2;
    }

    if( num_poses < 1 )
    {
        std::cerr << " Need a minimum of one deformed pose to Compute " << std::endl;
        return 2;
    }

    try
    {
        for( int i = 0 ; i < num_poses ; i++ )
        {
            Eigen::MatrixXd new_pose;
            rapidxml_ns::file<> xml_file(argv[4 + i]);  
            rapidxml_ns::xml_document<> doc;    // character type defaults to char
            doc.parse<rapidxml_ns::parse_no_string_terminators>(xml_file.data());  
            loadSvg(doc.first_node(), new_pose, ssdr_elem.df_tangents, ssdr_elem.df_curves);
            (ssdr_elem.frame_poses).push_back( new_pose );
        }

    }
    catch (std::exception const & e)
    {
        std::cerr << "Error loading other SVG poses: " << e.what() << std::endl;
        return 1;
    }

    // load the rest pose samples 
    try
    {
        rapidxml_ns::file<> xml_file(argv[4+num_poses]);  
        rapidxml_ns::xml_document<> doc;    // character type defaults to char
        Eigen::MatrixXd df_tangents_bs;
        std::vector<Eigen::Vector4i> df_curves_bs;

        doc.parse<rapidxml_ns::parse_no_string_terminators>(xml_file.data());  
        loadSvg(doc.first_node(), ssdr_elem.rest_pose_samples, df_tangents_bs ,df_curves_bs );

    }
    catch (std::exception const & e)
    {
        std::cerr << "Error loading rest Pose samples: " << e.what() << std::endl;
        return 1;
    }

    //load the deformed pose samples
    try
    {
        rapidxml_ns::file<> xml_file(argv[4+num_poses+1]);  
        rapidxml_ns::xml_document<> doc;    // character type defaults to char
        Eigen::MatrixXd df_tangents_bs;
        std::vector<Eigen::Vector4i> df_curves_bs;

        doc.parse<rapidxml_ns::parse_no_string_terminators>(xml_file.data());  
        loadSvg(doc.first_node(), ssdr_elem.deformed_pose_samples, df_tangents_bs ,df_curves_bs );

    }
    catch (std::exception const & e)
    {
        std::cerr << "Error loading Deformed Pose samples: " << e.what() << std::endl;
        return 1;
    }


    // testing
    int i = 1;
    for( auto elem : ssdr_elem.frame_poses )
    {
        // print the matrix
        std::cout << "\n\nPose\n\n" << elem << std::endl;
    }
     // print the rest pose matrix
     std::cout << "\n\nRest Pose\n\n" << ssdr_elem.rest_pose << std::endl;
     
     // print the rest pose samples matrix
     std::cout << "\n\nRest Pose Samples\n\n" << ssdr_elem.rest_pose_samples << std::endl;

     // print the deformed pose samples matrix
     std::cout << "\n\nDeformed Pose Samples\n\n" << ssdr_elem.deformed_pose_samples << std::endl;
    
    ssdr_elem.match_samples_to_curve( );
    ssdr_elem.normalize_tgs( );
    ssdr_elem.find_correspondence( );
    //ssdr_elem.init_bone_transforms( );
    //ssdr_elem.init_bone_transforms_2( );

    // testing shit

 static Point2 bezCurve[4] = {	/*  A cubic Bezier curve	*/
	{ 0.0, 0.0 },
	{ 1.0, 2.0 },
	{ 3.0, 3.0 },
	{ 4.0, 2.0 },
    };
    static Point2 arbPoint = { 3.5, 2.0 }; /*Some arbitrary point*/
    Point2	pointOnCurve;		 /*  Nearest point on the curve */
    double t_param = 0; 
    /*  Find the closest point */
    pointOnCurve = NearestPointOnCurve(arbPoint, bezCurve, &t_param);
    printf("pointOnCurve : (%4.4f, %4.4f)\n", pointOnCurve.x,
		pointOnCurve.y);
    printf("parameter of the point : %lf\n",t_param);

    return 0;
}
