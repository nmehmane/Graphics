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
void Context::transform_matrix(const boost::array<double, 6> & matrix)
{

}
void Context::set_viewport(double viewport_x, double viewport_y, double viewport_width, double viewport_height)
{
    std::cout << " I'm in view port " << viewport_x << "  " << viewport_y << std::endl;
    std::cout << viewport_width << "   " << viewport_height << std::endl; 
    min_x = viewport_x;
    min_y = viewport_y;
    max_x = viewport_width;
    max_y = viewport_height;
}

void Context::set_viewbox_size(double viewbox_width, double viewbox_height)
{
    std::cout << " view box " << viewbox_width << "  " << viewbox_height << std::endl;
    box_width = viewbox_width;
    box_height = viewbox_height;

}
void Context::disable_rendering()
{

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

void loadSvg( xml_element_t xml_root_element, 
              std::vector<ssdr_ns::Point>& curveEndPoints, 
              std::vector<ssdr_ns::Point>& curveMiddlePoints, 
              std::vector<Eigen::Vector4i>& curves )
{
    Context context;
    document_traversal<
      processed_elements<processed_elements_t>,
//      processed_attributes<traits::shapes_attributes_by_element>
      processed_attributes<processed_attributes_t>,
      viewport_policy<policy::viewport::as_transform>
    >::load_document(xml_root_element, context);

    std::vector<double>::iterator it;
    int row = 0;
    for( it = context.vertices.begin(); it != context.vertices.end(); it++,row++  )
    {
        ssdr_ns::Point p(*it,*(++it));
         curveEndPoints.push_back(p);
    }

    row = 0;
    for( it = context.v_tangents.begin(); it != context.v_tangents.end(); it++,row++  )
    {
        ssdr_ns::Point p(*it,*(++it));
         curveMiddlePoints.push_back(p);
    }
    
    std::vector<Eigen::Vector4i>::iterator it2;
    row = 0;
    for( it2 = context.curves.begin(); it2 != context.curves.end(); it2++)
    {
        curves.push_back(*it2);
        row++;
    }
}

int main(int argc, char** argv)
{
    if( argc < 7 )
    {
        std::cout <<"./exe <number of handles> <svg rest-pose>" 
                  <<" <number of frames> <svg pose 1> <svg pose2> ..."
                  <<" <rest-pose samples> <deformed-pose samples>" 
                  << std::endl;
        return -1;
    }
    
    //create ssdr element 
    ssdr_ns::ssdr ssdr_elem;
    
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
        rapidxml_ns::xml_document<> doc;    
        doc.parse<rapidxml_ns::parse_no_string_terminators>(xml_file.data());  
        loadSvg(doc.first_node(), 
                ssdr_elem.rp_CurveEndPoints,
                ssdr_elem.rp_CurveMiddlePoints, 
                ssdr_elem.rp_Curves);

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
        rapidxml_ns::file<> xml_file(argv[4]);  
        rapidxml_ns::xml_document<> doc;
        doc.parse<rapidxml_ns::parse_no_string_terminators>(xml_file.data());  
        loadSvg(doc.first_node(),
                ssdr_elem.dfp_CurveEndPoints, 
                ssdr_elem.dfp_CurveMiddlePoints,
                ssdr_elem.dfp_Curves);

    }
    catch (std::exception const & e)
    {
        std::cerr << "Error loading other SVG poses: " << e.what() << std::endl;
        return 1;
    }

    // print the rest pose
    std::cout << "\n\nRest Pose Curve End Points\n\n" << std::endl;

    for( auto& p : ssdr_elem.rp_CurveEndPoints )
    {
        std::cout << p << std::endl;
    }
    std::cout << "\n\nRest Pose Curve Middle Points\n\n" << std::endl;

    for( auto& p : ssdr_elem.rp_CurveMiddlePoints )
    {
        std::cout << p << std::endl;
    }
    
    // print the deformed pose
    std::cout << "\n\nDeformed Pose End Points\n\n" << std::endl;

    for( auto& p : ssdr_elem.dfp_CurveEndPoints )
    {
        std::cout << p << std::endl;
    }
    std::cout << "\n\nDeformed Pose Curve Middle Points\n\n" << std::endl;

    for( auto& p : ssdr_elem.dfp_CurveMiddlePoints )
    {
        std::cout << p << std::endl;
    }
     
    
    ssdr_elem.perform_ssdr();
    // testing shit

    return 0;
}
