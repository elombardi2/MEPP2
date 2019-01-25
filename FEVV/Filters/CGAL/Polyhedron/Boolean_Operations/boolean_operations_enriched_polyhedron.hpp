//TODO-elo-licence
#pragma once

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_3.h>

#include "boolean_operations_properties.h"


template <class Refs, class T, class P, class Norm>
class EnrichedVertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
public:
  // life cycle
  EnrichedVertex() {}

  // repeat mandatory constructors
  EnrichedVertex(const P &pt)
      : CGAL::HalfedgeDS_vertex_base< Refs, T, P >(pt)
  {
  }

  // vertex properties

  // subdivision:
  /*! \brief true if the vertex has been created*/
  //TODO-elo-rm-not-used?  bool Isnew;
  // boolean operations:
  /*! \brief An Id for the vertex*/
  VertexId Label;
};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class EnrichedHalfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
public:
  // life cycle
  EnrichedHalfedge() {}

  // halfedge properties

  // subdivision:
  /*! \brief true if the halfedge has been created or subdivided*/
  //TODO-elo-rm-not-used?  bool Isnew;
  // boolean operations:
  /*! \brief An Id for the halfedge*/
  HalfedgeId Label;
};


template <class Refs, class T, class Norm>
class EnrichedFacet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
public:
  // life cycle
  EnrichedFacet() {}

  // face properties

  // subdivision:
  /*! \brief true if the facet has been subdivided*/
  //TODO-elo-rm-not-used?  bool Issub;
  // operations booleennes :
  /*! \brief true if the facet belongs to the result*/
  bool IsExt;
  /*! \brief true if the facet has been processed*/
  bool IsOK;
  /*! \brief An Id for the facet*/
  FacetId Label;
};


struct EnrichedItems : public CGAL::Polyhedron_items_3
{
  // wrap vertex
  template< class Refs, class Traits >
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef EnrichedVertex< Refs, CGAL::Tag_true, Point, Normal > Vertex;
  };

  // wrap face
  template< class Refs, class Traits >
  struct Face_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef typename Traits::Plane_3 Plane;
    typedef EnrichedFacet< Refs, CGAL::Tag_true, Point, Normal, Plane > Face;
  };

  // wrap halfedge
  template< class Refs, class Traits >
  struct Halfedge_wrapper
  {
    typedef typename Traits::Vector_3 Normal;
    typedef EnrichedHalfedge< Refs,
                              CGAL::Tag_true,
                              CGAL::Tag_true,
                              CGAL::Tag_true,
                              Normal >
        Halfedge;
  };
};


#if 0 //TODO-elo-rm-section
using CGALKernel = CGAL::Cartesian< double >;
using MeshPolyhedron =
    CGAL::Polyhedron_3< CGALKernel, CGAL::Polyhedron_items_with_id_3 >;
typedef Enriched_polyhedron<Enriched_kernel, Enriched_items>	Polyhedron;
#endif

using CGALKernel = CGAL::Cartesian< double >;
using EnrichedPolyhedron =
    CGAL::Polyhedron_3< CGALKernel, EnrichedItems >;
