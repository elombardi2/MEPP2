// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include "FEVV/Filters/CGAL/Polyhedron/Boolean_Operations/boolpolyhedra.hpp"

namespace FEVV {
namespace Filters {


/**
 * TODO-elo-fix-this-header
 * \brief  Computes the union of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 * \param  g1      1st input mesh 
 * \param  pm1     point map of 1st mesh
 * \param  g2      2nd input mesh 
 * \param  pm2     point map of 2nd mesh
 * \param  g_out   output mesh 
 * \param  pm_out  point map of output mesh
 * \param  gt      the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 */
//TODO-elo-note: here HalfedgeGraph mus be a Polyhedron_3!
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_union(HalfedgeGraph &g1,
              PointMap      &pm1,
              HalfedgeGraph &g2,
              PointMap      &pm2,
              HalfedgeGraph &g_out,
              PointMap      &pm_out, //TODO-elo-really-necessary?
              const GeometryTraits &gt)
{
  BoolPolyhedra< HalfedgeGraph, PointMap >(
      &g1, &pm1, &g2, &pm2, &g_out, &pm_out, UNION);
}

/**
 * TODO-elo-fix-this-header
 * \brief  Computes the union of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  g1      1st input mesh 
 * \param  pm1     point map of 1st mesh
 * \param  g2      2nd input mesh 
 * \param  pm2     point map of 2nd mesh
 * \param  g_out   output mesh 
 * \param  pm_out  point map of output mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_union(HalfedgeGraph &g1,
              PointMap      &pm1,
              HalfedgeGraph &g2,
              PointMap      &pm2,
              HalfedgeGraph &g_out,
              PointMap      &pm_out //TODO-elo-really-necessary?
              )
{
  GeometryTraits gt(g1);
  boolean_union< HalfedgeGraph, PointMap, GeometryTraits >(
      g1, pm1, g2, pm2, g_out, pm_out, gt);
}


} // namespace Filters
} // namespace FEVV

