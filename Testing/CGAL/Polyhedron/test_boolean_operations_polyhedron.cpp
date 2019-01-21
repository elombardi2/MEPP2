// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/properties_polyhedron_3.h"

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include <CGAL/boost/graph/helpers.h>         // for CGAL::clear(mesh)

#include "FEVV/Filters/CGAL/Polyhedron/Boolean_Operations/boolean_operations.hpp"

int main(int argc, const char **argv)
{
  if(argc != 3)
  {
    std::cout
        << "Apply boolean union, intersection and subtraction on two meshes."
        << std::endl;
    std::cout << "Usage:  " << argv[0] << "  mesh_file_1  mesh_file_2"
              << std::endl;
    std::cout << "Example:  " << argv[0]
              << "  ../Testing/Data/CubeNonTriangleFaces.off  "
                 "../Testing/Data/tetra.off"
              << std::endl;
    return EXIT_FAILURE;
  }

  // input and output files
  std::string input_file_path_1 = argv[1];
  std::string input_file_path_2 = argv[2];
  std::string output_file_union_path =
      "test_boolean_operations_polyhedron.union.off";
  std::string output_file_inter_path =
      "test_boolean_operations_polyhedron.inter.off";
  std::string output_file_minus_path =
      "test_boolean_operations_polyhedron.minus.off";

  // read meshes from files
  FEVV::MeshPolyhedron m1;
  FEVV::PMapsContainer pmaps_bag_1;
  FEVV::Filters::read_mesh(input_file_path_1, m1, pmaps_bag_1);
  FEVV::MeshPolyhedron m2;
  FEVV::PMapsContainer pmaps_bag_2;
  FEVV::Filters::read_mesh(input_file_path_2, m2, pmaps_bag_2);

  // retrieve point maps
  auto pm1 = get(boost::vertex_point, m1);
  auto pm2 = get(boost::vertex_point, m2);

  // apply union filter
  {
    // create output mesh
    FEVV::MeshPolyhedron m_out;
    auto pm_out = get(boost::vertex_point, m_out);
    FEVV::PMapsContainer pmaps_bag_out;

    // apply union filter, result in m_out
    FEVV::Filters::boolean_union(m1, pm1, m2, pm2, m_out, pm_out);

    // write result to file
    FEVV::Filters::write_mesh(output_file_union_path, m_out, pmaps_bag_out);
  }

#if 0
  // apply intersection filter
  {
    // create output mesh
    FEVV::MeshPolyhedron m_out;
    auto pm_out = get(boost::vertex_point, m_out);
    FEVV::PMapsContainer pmaps_bag_out;

    // apply intersection filter, result in m_out
    FEVV::Filters::boolean_inter(m1, pm1, m2, pm2, m_out, pm_out);

    // write result to file
    FEVV::Filters::write_mesh(output_file_inter_path, m_out, pmaps_bag_out);
  }

  // apply subtraction filter
  {
    // create output mesh
    FEVV::MeshPolyhedron m_out;
    auto pm_out = get(boost::vertex_point, m_out);
    FEVV::PMapsContainer pmaps_bag_out;

    // apply subtraction filter, result in m_out
    FEVV::Filters::boolean_minus(m1, pm1, m2, pm2, m_out, pm_out);

    // write result to file
    FEVV::Filters::write_mesh(output_file_minus_path, m_out, pmaps_bag_out);
  }
#endif

  return 0;
}
