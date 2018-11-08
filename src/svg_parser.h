#ifndef SVG_PARSER_H
#define SVG_PARSER_H

#include <iostream>
#include <Eigen/Core>
#include "rapidxml_ns/rapidxml_ns.hpp"
#include "rapidxml_ns/rapidxml_ns_utils.hpp"
#include "svgpp/policy/xml/rapidxml_ns.hpp"
#include "svgpp/svgpp.hpp"
#include "ssdr.h"

using namespace svgpp;

typedef rapidxml_ns::xml_node<> const * xml_element_t;

class Context
{
public:

  void on_enter_element(tag::element::any);
  void on_exit_element();
  // converts relative coordinates to absolute coordinates if (m)
  void path_move_to(double x, double y, tag::coordinate::absolute);
  void path_line_to(double x, double y, tag::coordinate::absolute);
  void path_cubic_bezier_to(double x1, double y1, double x2, double y2,
                            double x, double y, tag::coordinate::absolute);
  void path_quadratic_bezier_to(double x1, double y1, double x, double y,
                                tag::coordinate::absolute);

  void path_elliptical_arc_to(double rx, double ry, double x_axis_rotation,
                              bool large_arc_flag, bool sweep_flag,
                              double x, double y, tag::coordinate::absolute);

  void path_close_subpath();

  void path_exit();
public:
  std::vector<double> vertices; // will hold the vertices as [ x0 y0 x1 y1 ... ]
  std::vector<double> v_tangents; // will hold the tangents as [ x'0 y'0 ...]
  std::vector<Eigen::Vector4i> curves; // will hold the index of the 4 control point of the bezier curve as [i0 i'0 i'1 i1]
  bool path_start;
};

typedef 
  boost::mpl::set<
    // SVG Structural Elements
    tag::element::svg,
    tag::element::g,
    // SVG Shape Elements
    tag::element::circle,
    tag::element::ellipse,
    tag::element::line,
    tag::element::path,
    tag::element::polygon,
    tag::element::polyline,
    tag::element::rect
  >::type processed_elements_t;

void loadSvg(xml_element_t xml_root_element,
             std::vector<ssdr_ns::Point>& rest_pose,
             std::vector<ssdr_ns::Point>& rp_tangents, 
             std::vector<Eigen::Vector4i>& rp_curves);

#endif
