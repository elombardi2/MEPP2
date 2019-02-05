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

#include "boolops_polyhedra.hpp"

namespace FEVV {
namespace Filters {

//--------------------- UNION -------------------------

/**
 * \brief  Computes the union of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 * \param  gA      1st input mesh
 * \param  pmA     point map of 1st mesh
 * \param  gB      2nd input mesh
 * \param  pmB     point map of 2nd mesh
 * \param  g_out   output mesh
 * \param  pm_out  point map of output mesh
 * \param  gt      the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_union(HalfedgeGraph &gA,
              PointMap      &pmA,
                // for compliance with filter policy, not used
              HalfedgeGraph &gB,
              PointMap      &pmB,
                // for compliance with filter policy, not used
              HalfedgeGraph &g_out,
              PointMap      &pm_out,
                // for compliance with filter policy, not used
              const GeometryTraits &gt)
{
  BoolPolyhedra< HalfedgeGraph, PointMap >(
      gA, gB, g_out, UNION);
}

/**
 * \brief  Computes the union of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  gA      1st input mesh
 * \param  pmA     point map of 1st mesh
 * \param  gB      2nd input mesh
 * \param  pmB     point map of 2nd mesh
 * \param  g_out   output mesh
 * \param  pm_out  point map of output mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_union(HalfedgeGraph &gA,
              PointMap      &pmA,
                // for compliance with filter policy, not used
              HalfedgeGraph &gB,
              PointMap      &pmB,
                // for compliance with filter policy, not used
              HalfedgeGraph &g_out,
              PointMap      &pm_out
                // for compliance with filter policy, not used
              )
{
  GeometryTraits gt(gA);
  boolean_union< HalfedgeGraph, PointMap, GeometryTraits >(
      gA, pmA, gB, pmB, g_out, pm_out, gt);
}


//--------------------- INTERSECTION -------------------------

/**
 * \brief  Computes the intersection of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 * \param  gA      1st input mesh
 * \param  pmA     point map of 1st mesh
 * \param  gB      2nd input mesh
 * \param  pmB     point map of 2nd mesh
 * \param  g_out   output mesh
 * \param  pm_out  point map of output mesh
 * \param  gt      the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_inter(HalfedgeGraph &gA,
              PointMap      &pmA,
                // for compliance with filter policy, not used
              HalfedgeGraph &gB,
              PointMap      &pmB,
                // for compliance with filter policy, not used
              HalfedgeGraph &g_out,
              PointMap      &pm_out,
                // for compliance with filter policy, not used
              const GeometryTraits &gt)
{
  BoolPolyhedra< HalfedgeGraph, PointMap >(
      gA, gB, g_out, INTER);
}

/**
 * \brief  Computes the intersection of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  gA      1st input mesh
 * \param  pmA     point map of 1st mesh
 * \param  gB      2nd input mesh
 * \param  pmB     point map of 2nd mesh
 * \param  g_out   output mesh
 * \param  pm_out  point map of output mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_inter(HalfedgeGraph &gA,
              PointMap      &pmA,
                // for compliance with filter policy, not used
              HalfedgeGraph &gB,
              PointMap      &pmB,
                // for compliance with filter policy, not used
              HalfedgeGraph &g_out,
              PointMap      &pm_out
                // for compliance with filter policy, not used
              )
{
  GeometryTraits gt(gA);
  boolean_inter< HalfedgeGraph, PointMap, GeometryTraits >(
      gA, pmA, gB, pmB, g_out, pm_out, gt);
}


//--------------------- SUBTRACTION -------------------------

/**
 * \brief  Computes the subtraction of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 * \param  gA      1st input mesh
 * \param  pmA     point map of 1st mesh
 * \param  gB      2nd input mesh
 * \param  pmB     point map of 2nd mesh
 * \param  g_out   output mesh
 * \param  pm_out  point map of output mesh
 * \param  gt      the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_minus(HalfedgeGraph &gA,
              PointMap      &pmA,
                // for compliance with filter policy, not used
              HalfedgeGraph &gB,
              PointMap      &pmB,
                // for compliance with filter policy, not used
              HalfedgeGraph &g_out,
              PointMap      &pm_out,
                // for compliance with filter policy, not used
              const GeometryTraits &gt)
{
  BoolPolyhedra< HalfedgeGraph, PointMap >(
      gA, gB, g_out, MINUS);
}

/**
 * \brief  Computes the subtraction of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  gA      1st input mesh
 * \param  pmA     point map of 1st mesh
 * \param  gB      2nd input mesh
 * \param  pmB     point map of 2nd mesh
 * \param  g_out   output mesh
 * \param  pm_out  point map of output mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_minus(HalfedgeGraph &gA,
              PointMap      &pmA,
                // for compliance with filter policy, not used
              HalfedgeGraph &gB,
              PointMap      &pmB,
                // for compliance with filter policy, not used
              HalfedgeGraph &g_out,
              PointMap      &pm_out
                // for compliance with filter policy, not used
              )
{
  GeometryTraits gt(gA);
  boolean_minus< HalfedgeGraph, PointMap, GeometryTraits >(
      gA, pmA, gB, pmB, g_out, pm_out, gt);
}


} // namespace Filters
} // namespace FEVV

