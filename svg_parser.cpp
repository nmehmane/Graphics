#include "svg_parser.h"

#include <iostream>
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
    rest_pose = Eigen::MatrixXd::Zero(context.vertices.size() ,2);

    std::vector<double>::iterator it;
    int row = 0;
    for( it = context.vertices.begin(); it != context.vertices.end(); it++,row++  )
    {
        rest_pose(row,0) = *it;
        printf(" x = %lf", rest_pose(row,0)); 
        rest_pose(row,1) = *(++it);
        printf(" y = %lf\n", rest_pose(row,1)); 
    }
}

int main(int argc, char** argv)
{
    // change args to ./exe <number of handles> <svg rest-pose> <number of frames> <svg pose 1> <svg pose2> ...
    if( argc != 2 )
    {
        std::cout<<"/executable <svg file name>"<<std::endl;
        return -1;
    }
    //create ssdr element 
    // fill num handles
    // fill the rest pose
    // fill the other poses
    Eigen::MatrixXd rest_pose;
    try
    {
        rapidxml_ns::file<> xml_file(argv[1]);  
        rapidxml_ns::xml_document<> doc;    // character type defaults to char
        doc.parse<rapidxml_ns::parse_no_string_terminators>(xml_file.data());  
        loadSvg(doc.first_node(), rest_pose);

    }
    catch (std::exception const & e)
    {
        std::cerr << "Error loading SVG: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
