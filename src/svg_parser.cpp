#include "svg_parser.h"
#include "ssdr.h"

#include <iostream>
#include <sstream>
#include <Eigen/Core>
#include "rapidxml_ns/rapidxml_ns.hpp"
#include "rapidxml_ns/rapidxml_ns_utils.hpp"
#include "svgpp/policy/xml/rapidxml_ns.hpp"
#include "svgpp/svgpp.hpp"

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
      
    //std::cout<<" path_move_to " << std::endl;
    //printf(" x = %lf y = %lf \n",x,y);
    vertices.push_back(x);
    vertices.push_back(y);
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

void loadSvg(xml_element_t xml_root_element, Eigen::MatrixXd& rest_pose)
{
    Context context;
    document_traversal<
      processed_elements<processed_elements_t>,
      processed_attributes<traits::shapes_attributes_by_element>
    >::load_document(xml_root_element, context);
    //  an N x 2 matrix
    rest_pose = Eigen::MatrixXd::Zero((context.vertices.size()/2) ,2);

    std::vector<double>::iterator it;
    int row = 0;
    for( it = context.vertices.begin(); it != context.vertices.end(); it++,row++  )
    {
        rest_pose(row,0) = *it;
        //printf(" x = %lf", rest_pose(row,0)); 
        rest_pose(row,1) = *(++it);
        //printf(" y = %lf\n", rest_pose(row,1)); 
    }
}

int main(int argc, char** argv)
{
    // change args to ./exe <number of handles> <svg rest-pose> <number of frames> <svg pose 1> <svg pose2> ...
    // Minimum requirements: 1)number of handels 2)rest pose 3)number of frames ( minimum 2 ) 4,5)minimum of two frames
    if( argc < 5 )
    {
        std::cout<<"./exe <number of handles> <svg rest-pose> <number of frames> <svg pose 1> <svg pose2> ..."<<std::endl;
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
        loadSvg(doc.first_node(), ssdr_elem.rest_pose);

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
    if( (argc - 4) != num_poses )
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
            loadSvg(doc.first_node(), new_pose);
            (ssdr_elem.frame_poses).push_back( new_pose );
        }

    }
    catch (std::exception const & e)
    {
        std::cerr << "Error loading other SVG poses: " << e.what() << std::endl;
        return 1;
    }

    // testing
    int i = 1;
    for( auto elem : ssdr_elem.frame_poses )
    {
        // print the matrix
        std::cout << "Pose\n" << elem << std::endl;
    }
     // print the rest pose matrix
     std::cout << "Rest Pose\n" << ssdr_elem.rest_pose << std::endl;

    //ssdr_elem.init_bone_transforms( );
    ssdr_elem.init_bone_transforms_2( );


    return 0;
}
