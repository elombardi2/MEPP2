//TODO-elo-licence
#pragma once

/*!
 * \file BoolPolyhedra.h
 * \author Cyril Leconte
 */

#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include "boolean_operations_definitions.hpp"
//TODO-elo-restore  #include "CPolyhedron_from_polygon_builder_3.h"
//TODO-elo-restore  #include "Boolean_Operations_triangulation.h"

#include <boost/graph/graph_traits.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"
#include <list>

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
#include "Time_measure.h"
#endif // BOOLEAN_OPERATIONS_DEBUG

/*! \typedef HDS
 * \brief Halfedge data structure
 */
//TODO-elo-rm  typedef Polyhedron::HalfedgeDS                HDS;

/*! \typedef AABB_Kernel
 * \brief Kernel used for the computations in a AABB-tree
 */
//typedef CGAL::Simple_cartesian<num_type>          AABB_Kernel;
typedef CGAL::Simple_cartesian<double>           AABB_Kernel;

/*! \class Enriched_Triangle
 *  \brief An enriched triangle
 */
template< typename HalfedgeGraph >
class Enriched_Triangle : public AABB_Kernel::Triangle_3
{
public:
  typedef boost::graph_traits< HalfedgeGraph >   GraphTraits;
  typedef typename GraphTraits::face_descriptor  Facet_handle;

  /*! \typedef Point_3
   *  \brief 3d point
   */
  typedef AABB_Kernel::Point_3  Point_3;

  /*! \typedef FT
   *  \brief number type used in the kernel
   */
  typedef AABB_Kernel::FT  FT;

  /*! \brief Constructor
   *  \brief creates a triangle a little bigger
   *  \param _f : handle of a triangular facet of a polyhedron
   */
  Enriched_Triangle(Facet_handle &_f)
      : AABB_Kernel::Triangle_3(
            to_K(_f->facet_begin()->vertex()->point() +
                  (_f->facet_begin()->vertex()->point() -
                  _f->facet_begin()->next()->vertex()->point()) /
                      1000),
            to_K(_f->facet_begin()->next()->vertex()->point() +
                  (_f->facet_begin()->next()->vertex()->point() -
                  _f->facet_begin()->next()->next()->vertex()->point()) /
                      1000),
            to_K(_f->facet_begin()->next()->next()->vertex()->point() +
                  (_f->facet_begin()->next()->next()->vertex()->point() -
                  _f->facet_begin()->vertex()->point()) /
                      1000)),
        f(_f)
  {
  }

  /*! \brief Accessor
   *  \return The handle of the facet used to build the triangle
   */
  Facet_handle facet() { return f; }

  /*! \brief convert any 3d point in the kernel used by the AABB-tree
   *  \param p : The point to convert
   *  \return The 3d point converted
   */
  template < typename Point3d >
  inline Point_3 to_K(Point3d p)
  {
    return Point_3((FT)p.x(), (FT)p.y(), (FT)p.z());
  } // MT: suppression référence //TODO-elo: restore reference ?

private:
  /*! \brief The handle of the facet used to build the triangle*/
  Facet_handle f;
};

/*! \class BoolPolyhedra
 * \brief The class that compute a Boolean operation*/
template< typename HalfedgeGraph, typename PointMap >
class BoolPolyhedra
{
private:
  typedef boost::graph_traits< HalfedgeGraph >      GraphTraits;
  typedef typename GraphTraits::face_descriptor     Facet_handle;
  typedef typename GraphTraits::halfedge_descriptor Halfedge_handle;
  //TODO-elo-why-not-working?  typedef typename GraphTraits::vertex_iterator   Vertex_iterator;
  typedef typename HalfedgeGraph::Vertex_iterator   Vertex_iterator;
  //TODO-elo-why-not-working?  typedef typename GraphTraits::face_iterator    Facet_iterator;
  typedef typename HalfedgeGraph::Facet_iterator    Facet_iterator;


  /*! \typedef Triangle
  * \brief A triangle enriched with a facet handle*/
  using Triangle = Enriched_Triangle< HalfedgeGraph >;

  /*! \typedef AABB_Primitive
  * \brief A primitive for an AABB-tree*/
  using TriangleListIterator = typename std::list< Triangle >::iterator;
  using AABB_Primitive =
      typename CGAL::AABB_triangle_primitive< AABB_Kernel,
                                              TriangleListIterator >;
  /*! \typedef AABB_Traits
  * \brief concept for AABB-tree*/
  typedef CGAL::AABB_traits< AABB_Kernel, AABB_Primitive > AABB_Traits;
  /*! \typedef AABB_Tree
  * \brief AABB-tree*/
  typedef CGAL::AABB_tree< AABB_Traits > AABB_Tree;

  /*!
  * \typedef VertexId
  * \brief Vertex Id
  */
  typedef unsigned long              VertexId;

  /*!
  * \typedef HalfedgeId
  * \brief Halfedge Id
  */
  typedef unsigned long              HalfedgeId;

  /*!
  * \typedef FacetId
  * \brief Facet Id
  */
  typedef unsigned long              FacetId;

  /*!
  * \typedef InterId
  * \brief Intersection Id
  */
  typedef unsigned long              InterId;

  /*! \struct Triangle_Cut
   *  \brief A structure containing informations about an intersected facet
   */
  struct Triangle_Cut
  {
    /*! \brief true if the facet belongs to the first polyhedron */
    bool Facet_from_A;
    /*! \brief An exact vector giving the direction of the normal */
    Vector_exact norm_dir;
    /*! \brief A list of segments (the intersections with the facets of the
     * other polyhedron) */
    std::vector< std::vector< InterId > > CutList;
    /*! \brief A list of points (when the intersection is a point) */
    std::set< InterId > PtList;
    /*! \brief The list of the intersections */
    std::map< HalfedgeId, InterId > RefInter;

    /*! \brief Default constructor */
    Triangle_Cut() {}
    /*! \brief Constructor
     \param V : The normal direction
     \param ffA : Must be true if the facet belongs to the first polyhedron*/
    Triangle_Cut(Vector_exact V, bool ffA)
    {
      norm_dir = V;
      Facet_from_A = ffA;
    } // MT
  };

  /*! \struct Info_Inter
   *  \brief Contains informations about an intersection between a facet and a
   *         halfedge
   */
  struct Info_Inter
  {
    /*! \brief The facet*/
    Facet_handle f;
    /*! \brief The halfedge*/
    Halfedge_handle he;
    /*! \brief true if the intersection is exactly on the vertex pointed by he*/
    bool IsOnVertex;
    /*! \brief A code for the location of the intersection :\n\n
     * 0 : Intersection is strictly in the facet\n
     * 1 : Intersection is on the first edge of the facet\n
     * 2 : Intersection is on the second edge of the facet\n
     * 3 : Intersection is exactly on the first vertex of the facet\n
     * 4 : Intersection is on the third edge of the facet\n
     * 5 : Intersection is exactly on the third vertex of the facet\n
     * 6 : Intersection is exactly on the second vertex of the facet\n
     * 7 : There is no intersection */
    unsigned short res;
    /*! \brief The intersection point (exact)*/
    Point3d_exact pt;
    /*! \brief The Id of the intersection point*/
    InterId Id;
  };

public:
  /*! \brief Constructor.
   * \brief Computes a boolean operation
   * \param pMA : The first polyhedron
   * \param pMB : The second polyhedron
   * \param pMout : The result polyhedron
   * \param BOOP : The Boolean operator. Must be UNION, INTER or MINUS*/
  BoolPolyhedra(HalfedgeGraph *pMA, //TODO-elo-rename to mA
                PointMap      *pmA,
                HalfedgeGraph *pMB, //TODO-elo-rename to mB
                PointMap      *pmB,
                HalfedgeGraph *pMout,
                PointMap      *pm_out,
                Bool_Op BOOP)
      : m_BOOP(BOOP)
  {

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
    std::ofstream ofstrMA("input_A_boolsum.off"); ofstrMA << *pMA;
    std::ofstream ofstrMB("input_B_boolsum.off"); ofstrMB << *pMB;
    Time_measure Timer_total, Timer;

    Timer_total.Start();
    Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

    Init(pMA, pmA, pMB, pmB);
    //
    /////////////////////////////////////////////////
    //                                             //
    //        WIP       WIP   WIP   WIPWIP         //
    //        WIP       WIP         WIP WIP        //
    //        WIP       WIP   WIP   WIPWIP         //
    //         WIP WIP WIP    WIP   WIP            //
    //           WIP WIP      WIP   WIP            //
    //                                             //
    /////////////////////////////////////////////////
    //
    //
    //
    //
#if 0 //TODO-elo-WIP

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
    duration_Init = Timer.GetDiff();
    Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

    FindCouples();

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
    duration_FindCouples = Timer.GetDiff();
    Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

    if(!m_Couples.empty())
    {
      ComputeIntersections();

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
      duration_ComputeIntersections = Timer.GetDiff();
      Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

      CutIntersectedFacets();

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
      duration_CutIntersectedFacets = Timer.GetDiff();
      Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

      PropagateFacets();

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
      duration_PropagateFacets = Timer.GetDiff();
      Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

      pMout->delegate(ppbuilder);

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG
      duration_delegate = Timer.GetDiff();
      duration_total = Timer_total.GetDiff();
      std::ofstream ofstrMS("output_boolsum.off"); ofstrMS << *pMout;
      ofstrtime.open("output_time.txt");
      WriteData(pMout);
      ColorType();
#endif // BOOLEAN_OPERATIONS_DEBUG

    }
#endif //TODO-elo-WIP
  }

  /*! \brief Destructor*/
  ~BoolPolyhedra() {}

private:
  /*! \brief Initialisation of the tags, and triangulation of the two input polyhedra
   * \param pMA : The first polyhedron
   * \param pMB : The second polyhedron*/
  void Init(HalfedgeGraph *pMA, //TODO-elo-rename to mA
            PointMap      *pmA,
            HalfedgeGraph *pMB, //TODO-elo-rename to mB
            PointMap      *pmB)
  {
    m_pA = pMA;
    m_pmA = pmA;
    m_pB = pMB;
    m_pmB = pmB;

    // initialize property maps
    m_vertex_Label_A =
        FEVV::make_vertex_property_map< HalfedgeGraph, VertexId >(*m_pA);
    m_vertex_Label_B =
        FEVV::make_vertex_property_map< HalfedgeGraph, VertexId >(*m_pB);

    m_halfedge_Label_A =
        FEVV::make_halfedge_property_map< HalfedgeGraph, HalfedgeId >(*m_pA);
    m_halfedge_Label_B =
        FEVV::make_halfedge_property_map< HalfedgeGraph, HalfedgeId >(*m_pB);

    m_face_Label_A =
        FEVV::make_face_property_map< HalfedgeGraph, FacetId >(*m_pA);
    m_face_Label_B =
        FEVV::make_face_property_map< HalfedgeGraph, FacetId >(*m_pB);
    m_face_IsExt_A =
        FEVV::make_face_property_map< HalfedgeGraph, bool >(*m_pA);
    m_face_IsExt_B =
        FEVV::make_face_property_map< HalfedgeGraph, bool >(*m_pB);
    m_face_IsOK_A =
        FEVV::make_face_property_map< HalfedgeGraph, bool >(*m_pA);
    m_face_IsOK_B =
        FEVV::make_face_property_map< HalfedgeGraph, bool >(*m_pB);

    // triangulation of the two input polyhedra
    // this is necessary for the AABB-tree, and simplify the computation of the
    // intersections
    if(!m_pA->is_pure_triangle())
      triangulate(m_pA);
    if(!m_pB->is_pure_triangle())
      triangulate(m_pB);

    // initialize the tags
    for(Vertex_iterator pVertex = m_pA->vertices_begin();
        pVertex != m_pA->vertices_end();
        ++pVertex)
    {
      // TODO-elo-rm  pVertex->Label = 0xFFFFFFFF;
      //TODO-elo-test-rm-beg
      //int i = *pVertex;
      //typename GraphTraits::vertex_descriptor vd;
      //int j = vd;
      //typename HalfedgeGraph::Vertex_handle vh;
      //int k = vh;
      //int l = pVertex;
      //typename GraphTraits::vertex_iterator vi;
      //int z = vi;
      // -> vd = same type as vh and pVertex!!!
      // -> vi = CGAL::Prevent_deref<GraphTraits::vertex_descriptor>
      //TODO-elo-test-rm-end
      put(m_vertex_Label_A, pVertex, 0xFFFFFFFF);
    }
    for(Vertex_iterator pVertex = m_pB->vertices_begin();
        pVertex != m_pB->vertices_end();
        ++pVertex)
    {
      // TODO-elo-rm  pVertex->Label = 0xFFFFFFFF;
      put(m_vertex_Label_B, pVertex, 0xFFFFFFFF);
    }
    for(Facet_iterator pFacet = m_pA->facets_begin();
        pFacet != m_pA->facets_end();
        ++pFacet)
    {
      // TODO-elo-rm  pFacet->Label = 0xFFFFFFFF;
      // TODO-elo-rm  pFacet->IsExt = false;
      // TODO-elo-rm  pFacet->IsOK = false;
      put(m_face_Label_A, pFacet, 0xFFFFFFFF);
      put(m_face_IsExt_A, pFacet, false);
      put(m_face_IsOK_A, pFacet, false);
    }
    for(Facet_iterator pFacet = m_pB->facets_begin();
        pFacet != m_pB->facets_end();
        ++pFacet)
    {
      //TODO-elo-rm  pFacet->Label = 0xFFFFFFFF;
      //TODO-elo-rm  pFacet->IsExt = false;
      //TODO-elo-rm  pFacet->IsOK = false;
      put(m_face_Label_B, pFacet, 0xFFFFFFFF);
      put(m_face_IsExt_B, pFacet, false);
      put(m_face_IsOK_B, pFacet, false);
    }

    { // TODO-elo-rm-dbg-display
      // init HE property maps to be able to compare
      // with Mepp1
      typename HalfedgeGraph::Halfedge_iterator pHe;
      for(pHe = m_pA->halfedges_begin(); pHe != m_pA->halfedges_end(); pHe++)
        put(m_halfedge_Label_A, pHe, 42424242);
      for(pHe = m_pB->halfedges_begin(); pHe != m_pB->halfedges_end(); pHe++)
        put(m_halfedge_Label_B, pHe, 42424242);
    }
  }

  void triangulate(HalfedgeGraph *m) //TODO-elo-make-a-separate-filter?
  {
    { //TODO-elo-dbg-rm
      std::cout << "Triangulating mesh..." << std::endl;
      std::cout << "before triangulating mesh:" << std::endl;
      std::cout << "  size_of_vertices = " << m->size_of_vertices() << std::endl;
      std::cout << "  size_of_halfedges = " << m->size_of_halfedges() << std::endl;
      std::cout << "  size_of_facets = " << m->size_of_facets() << std::endl;
    }

    Facet_iterator f = m->facets_begin();
    Facet_iterator f2 = m->facets_begin();
    do // for (; f != this->facets_end(); f++)
    {
      f = f2;
      if(f == m->facets_end())
      {
        break;
      }
      f2++;

      if(!(f->is_triangle()))
      {
        int num = (int)(f->facet_degree() - 3);
        Halfedge_handle h = f->halfedge();

        h = m->make_hole(h);

        Halfedge_handle g = h->next();
        g = g->next();
        Halfedge_handle new_he = m->add_facet_to_border(h, g);
        g = new_he;

        num--;
        while(num != 0)
        {
          g = g->opposite();
          g = g->next();
          Halfedge_handle new_he = m->add_facet_to_border(h, g);
          g = new_he;

          num--;
        }

        m->fill_hole(h);
      }

    } while(true);

    //TODO-elo-rm  this->compute_normals();
    //TODO-elo-rm  this->compute_type();
    { //TODO-elo-dbg-rm
      std::cout << "after triangulating mesh:" << std::endl;
      std::cout << "  size_of_vertices = " << m->size_of_vertices() << std::endl;
      std::cout << "  size_of_halfedges = " << m->size_of_halfedges() << std::endl;
      std::cout << "  size_of_facets = " << m->size_of_facets() << std::endl;
    }
  }


#if 0 //TODO-elo-WIP
  /*! \brief Finds every couple of facets between the two input polyhedra that intersects
   * \brief Each couple is stored in the member m_Couples*/
  void FindCouples()
  {
    //A AABB-tree is built on the facets of one of the polyhedra. A collision test is done with each facet of the other polyhedron.
    Facet_iterator pFacet =  NULL;
    std::list<typename AABB_Tree::Primitive_id> primitives;
    std::list<Triangle> triangles;

    HalfedgeId i = 0;
    FacetId j = 0;

    //The AABB-tree is built on the polyhedron with the less number of facets
    if(m_pA->size_of_facets() < m_pB->size_of_facets())
    {
      //Building the AABB-tree on the first polyhedron
      for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
        triangles.push_back(Triangle(pFacet));
      tree.rebuild(triangles.begin(),triangles.end());

      //collision test with each facet of the second polyhedron
      for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
      {
        //"primitives" is the list of the triangles intersected (as a list of
        //triangles)
        tree.all_intersected_primitives(Triangle(pFacet),
                                        std::back_inserter(primitives));
        // TODO-elo-note: as far as I understand, primitives contains
        //                all A faces that intersect pFacet B face
        if(primitives.size() !=0)
        {
          // TODO-elo-note: push B face to m_Facet_Handle
          m_Facet_Handle.push_back(pFacet);
          //update of the tags (the facet and the three incidents halfedges
          //TODO-elo-rm  pFacet->Label = j++;
          put(m_face_Label_B, pFacet, j++);
          //TODO-elo-rm  pFacet->facet_begin()->Label = i++;
          put(m_halfedge_Label_B, pFacet->facet_begin(), i++);
          //TODO-elo-rm  pFacet->facet_begin()->next()->Label = i++;
          put(m_halfedge_Label_B, pFacet->facet_begin()->next(), i++);
          //TODO-elo-rm  pFacet->facet_begin()->next()->next()->Label = i++;
          put(m_halfedge_Label_B, pFacet->facet_begin()->next()->next(), i++);
          // creation of a Triangle_Cut structure to store the informations
          // about the intersections
          m_Inter_tri.push_back(Triangle_Cut(
              Compute_Normal_direction(pFacet->facet_begin()), false));
          do
          {
            //TODO-elo-note: supposing primitives contains A mesh faces
            //same operations for the intersected primitives (only one time)
            //TODO-elo-rm  if(primitives.back()->facet()->Label == 0xFFFFFFFF)
            if(get(m_face_Label_A, primitives.back()->facet()) == 0xFFFFFFFF)
            {
              //TODO-elo-note: push A face to m_Facet_Handle
              m_Facet_Handle.push_back(primitives.back()->facet());
              //TODO-elo-rm  primitives.back()->facet()->Label = j++;
              put(m_face_Label_A, primitives.back()->facet(), j++);
              //TODO-elo-rm  primitives.back()->facet()->facet_begin()->Label = i++;
              put(m_halfedge_Label_A, primitives.back()->facet()->facet_begin(), i++);
              //TODO-elo-rm  primitives.back()->facet()->facet_begin()->next()->Label = i++;
              put(m_halfedge_Label_A, primitives.back()->facet()->facet_begin()->next(), i++);
              //TODO-elo-rm  primitives.back()->facet()->facet_begin()->next()->next()->Label = i++;
              put(m_halfedge_Label_A, primitives.back()->facet()->facet_begin()->next()->next(), i++);
              m_Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(primitives.back()->facet()->facet_begin()), true));
            }
            //store every couple of intersected facet
            //TODO-elo-rm  m_Couples[primitives.back()->facet()->Label].insert(pFacet->Label);
            m_Couples[get(m_face_Label_A, primitives.back()->facet())].insert(
                get(m_face_Label_B, pFacet));
            //TODO-elo-note: here m_Couples contains (label_A, label_B) couples
            primitives.pop_back();
          }
          while(primitives.size() != 0);
          //TODO-elo-note: primitives is empty here
        }
      }
    }
    else
    {
      //Building the AABB-tree on the second polyhedron
      for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
        triangles.push_back(Triangle(pFacet));
      tree.rebuild(triangles.begin(),triangles.end());

      //collision test with each facet of the first polyhedron
      for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
      {
        //"primitives" is the list of the triangles intersected (as a list of triangles)
        tree.all_intersected_primitives(Triangle(pFacet),
                                        std::back_inserter(primitives));
        // TODO-elo-note: as far as I understand, primitives contains
        //                all B faces that intersect pFacet A face
        if(primitives.size() !=0)
        {
          // TODO-elo-note: push A face to m_Facet_Handle
          m_Facet_Handle.push_back(pFacet);
          //update of the tags (the facet and the three incidents halfedges
          //TODO-elo-rm  pFacet->Label = j++;
          put(m_face_Label_A, pFacet, j++);
          //TODO-elo-rm  pFacet->facet_begin()->Label = i++;
          put(m_halfedge_Label_A, pFacet->facet_begin(), i++);
          //TODO-elo-rm  pFacet->facet_begin()->next()->Label = i++;
          put(m_halfedge_Label_A, pFacet->facet_begin()->next(), i++);
          //TODO-elo-rm  pFacet->facet_begin()->next()->next()->Label = i++;
          put(m_halfedge_Label_A, pFacet->facet_begin()->next()->next(), i++);
          // creation of a Triangle_Cut structure to store the informations
          // about the intersections
          m_Inter_tri.push_back(Triangle_Cut(
              Compute_Normal_direction(pFacet->facet_begin()), true));
          do
          {
            //TODO-elo-note: supposing primitives contains B mesh faces
            //same operations for the intersected primitives (only one time)
            //TODO-elo-rm  if(primitives.back()->facet()->Label == 0xFFFFFFFF)
            if(get(m_face_Label_B, primitives.back()->facet()) == 0xFFFFFFFF)
            {
              //TODO-elo-note: push B face to m_Facet_Handle
              m_Facet_Handle.push_back(primitives.back()->facet());
              //TODO-elo-rm  primitives.back()->facet()->Label = j++;
              put(m_face_Label_B, primitives.back()->facet(), j++);
              //TODO-elo-rm  primitives.back()->facet()->facet_begin()->Label = i++;
              put(m_halfedge_Label_B, primitives.back()->facet()->facet_begin(), i++);
              //TODO-elo-rm  primitives.back()->facet()->facet_begin()->next()->Label = i++;
              put(m_halfedge_Label_B, primitives.back()->facet()->facet_begin()->next(), i++);
              //TODO-elo-rm  primitives.back()->facet()->facet_begin()->next()->next()->Label = i++;
              put(m_halfedge_Label_B, primitives.back()->facet()->facet_begin()->next()->next(), i++);
              m_Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(primitives.back()->facet()->facet_begin()), false));
            }
            //store every couple of intersected facet
            //TODO-elo-rm  m_Couples[pFacet->Label].insert(primitives.back()->facet()->Label);
            m_Couples[get(m_face_Label_A, pFacet)].insert(
                get(m_face_Label_B, primitives.back()->facet()));
            //TODO-elo-note: here m_Couples contains (label_A, label_B) couples
            primitives.pop_back();
          }
          while(primitives.size() != 0);
          //TODO-elo-note: primitives is empty here
        }
      }
    }

    { //TODO-elo-rm-dbg-display
      std::cout << "end of FindCouples(), mesh A, face Label property:" << std::endl;
      for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
          std::cout << get(m_face_Label_A, pFacet) << std::endl;
      std::cout << "end of FindCouples(), mesh B, face Label property:" << std::endl;
      for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
          std::cout << get(m_face_Label_B, pFacet) << std::endl;
      std::cout << "end of FindCouples(), mesh A, halfedge Label property:" << std::endl;
      typename HalfedgeGraph::Halfedge_iterator pHe;
      for(pHe = m_pA->halfedges_begin(); pHe != m_pA->halfedges_end(); pHe++)
          std::cout << get(m_halfedge_Label_A, pHe) << std::endl;
      std::cout << "end of FindCouples(), mesh B, halfedge Label property:" << std::endl;
      for(pHe = m_pB->halfedges_begin(); pHe != m_pB->halfedges_end(); pHe++)
          std::cout << get(m_halfedge_Label_B, pHe) << std::endl;
    }
  }

  /*! \brief Compute the intersections*/
  void ComputeIntersections()
  {
    while(!m_Couples.empty())
    {
      FacetId fA, fB;
      fA = m_Couples.begin()->first;
      fB = *m_Couples[fA].begin();
      InterTriangleTriangle(fA, fB);
      rmCouple(fA, fB);
    }
  }


  /*! \brief Cuts the intersected facets and starts to build the result*/
  void CutIntersectedFacets()
  {
    Triangle_Cut TriCut;
    Halfedge_handle he;

    //every intersected facet is triangulated if at least one of the intersections is a segment
    for(FacetId Facet = 0 ; Facet != m_Inter_tri.size() ; ++Facet)
    {
      if(!m_Inter_tri[Facet].CutList.empty())
      {
        TriCut = m_Inter_tri[Facet];
        he = m_Facet_Handle[Facet]->facet_begin();
        bool IsExt[3];
        //creation of a triangulation
        Triangulation<Exact_Kernel> T(he, TriCut.norm_dir);
        //add the list of intersection points (only happens in case of intersection of two edges)
        for(std::set<InterId>::iterator i = TriCut.PtList.begin();i != TriCut.PtList.end();++i)
        {
                                        T.add_new_pt(m_InterPts[*i], (unsigned long &)*i);    // MT: ajout cast
        }
        //add the intersection segments
                                for(int i = 0;i!=(int)TriCut.CutList.size();++i)
        {
          T.add_segment(m_InterPts[TriCut.CutList[i][0]], m_InterPts[TriCut.CutList[i][1]], TriCut.CutList[i][0], TriCut.CutList[i][1]);
        }
        //get the triangles of the triangulation thay belong to the result
        //and determine if the three neighboring facets belongs to the result (using IsExt[3])
                                vector<vector<unsigned long> > Tri_set = T.get_triangles((m_BOOP == MINUS && !TriCut.Facet_from_A)?true:false, IsExt);
        //add these triangles to the result
        ppbuilder.add_triangle(Tri_set, he);

        //update the tags
        m_Facet_Handle[Facet]->IsOK = true;
        if(IsExt[0]) he->opposite()->facet()->IsExt = true;
        if(IsExt[1]) he->next()->opposite()->facet()->IsExt = true;
        if(IsExt[2]) he->next()->next()->opposite()->facet()->IsExt = true;
      }
    }
  }

  /*! \brief Complete the building of the result*/
  void PropagateFacets()
  {
    Facet_handle pFacet = NULL, f = NULL, nf = NULL;
    stack<Facet_handle> tmpTriangles;

    //add to a stack the intersected facets that have been cut during CutIntersectedFacets
    for (pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
    {
      if(pFacet->IsOK) tmpTriangles.push(pFacet);
    }

    //while the stack is not empty, we look the three neighboring facets
    //if these facets has not been validated (IsOK == false), the facet is validated and added to the stack
    //if this facet is taged as a part of the result (IsExt == true), the facet is added to the result
    while(!tmpTriangles.empty())
    {
      f = tmpTriangles.top();
      tmpTriangles.pop();
      nf = f->facet_begin()->opposite()->facet();
      if(!nf->IsOK)
      {
        nf->IsOK = true;
        tmpTriangles.push(nf);
        if(nf->IsExt) add_facet_to_solution(nf, true);
      }
      nf = f->facet_begin()->next()->opposite()->facet();
      if(!nf->IsOK)
      {
        nf->IsOK = true;
        tmpTriangles.push(nf);
        if(nf->IsExt) add_facet_to_solution(nf, true);
      }
      nf = f->facet_begin()->next()->next()->opposite()->facet();
      if(!nf->IsOK)
      {
        nf->IsOK = true;
        tmpTriangles.push(nf);
        if(nf->IsExt) add_facet_to_solution(nf, true);
      }
    }

    //same process for the second polyhedron
    for (pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
    {
      if(pFacet->IsOK) tmpTriangles.push(pFacet);
    }

    while(!tmpTriangles.empty())
    {
      f = tmpTriangles.top();
      tmpTriangles.pop();
      nf = f->facet_begin()->opposite()->facet();
      if(!nf->IsOK)
      {
        nf->IsOK = true;
        tmpTriangles.push(nf);
        if(nf->IsExt) add_facet_to_solution(nf, false);
      }
      nf = f->facet_begin()->next()->opposite()->facet();
      if(!nf->IsOK)
      {
        nf->IsOK = true;
        tmpTriangles.push(nf);
        if(nf->IsExt) add_facet_to_solution(nf, false);
      }
      nf = f->facet_begin()->next()->next()->opposite()->facet();
      if(!nf->IsOK)
      {
        nf->IsOK = true;
        tmpTriangles.push(nf);
        if(nf->IsExt) add_facet_to_solution(nf, false);
      }
    }
  }

  /*! \brief removes properly a couple from the list
   * \param A : Id of the first facet
   * \param B : Id of the second facet
   */
  void rmCouple(FacetId &A, FacetId &B)
  {
    if(m_Couples[A].count(B) != 0) m_Couples[A].erase(B);
    if(m_Couples[A].empty()) m_Couples.erase(A);
  }

  /*! \brief Compute the intersection between two facets
   * \param A : Facet Id of the first facet (from the first polyhedron)
   * \param B : Facet Id of the second facet (from the second polyhedron)*/
  void InterTriangleTriangle(FacetId &A, FacetId &B)
  {
    Vector_exact nA, nB;
    nA = m_Inter_tri[A].norm_dir;
    nB = m_Inter_tri[B].norm_dir;

    Facet_handle fA, fB, fA2, fB2;
    fA = m_Facet_Handle[A];
    fB = m_Facet_Handle[B];

    Halfedge_handle heA[3], heB[3];
    heA[0] = fA->facet_begin();
    heA[1] = heA[0]->next();
    heA[2] = heA[1]->next();
    heB[0] = fB->facet_begin();
    heB[1] = heB[0]->next();
    heB[2] = heB[1]->next();

    Point3d_exact ptA[3], ptB[3];
    ptA[0] = point_to_exact(heA[0]->vertex()->point());
    ptA[1] = point_to_exact(heA[1]->vertex()->point());
    ptA[2] = point_to_exact(heA[2]->vertex()->point());
    ptB[0] = point_to_exact(heB[0]->vertex()->point());
    ptB[1] = point_to_exact(heB[1]->vertex()->point());
    ptB[2] = point_to_exact(heB[2]->vertex()->point());

    //compute the position of the three vertices of each triangle regarding the plane of the other
    //positive if the vertex is above
    //negative if the vertex is under
    //zero if the vertex is exactly on the triangle
    num_type posA[3], posB[3];
    posA[0] = nB * (ptA[0] - ptB[0]);
    posA[1] = nB * (ptA[1] - ptB[0]);
    posA[2] = nB * (ptA[2] - ptB[0]);
    posB[0] = nA * (ptB[0] - ptA[0]);
    posB[1] = nA * (ptB[1] - ptA[0]);
    posB[2] = nA * (ptB[2] - ptA[0]);

    //a code is computed on 6 bits using these results (two bits for each point)
    //10 -> above ; 01 -> under ; 00 -> on the plane
    unsigned short posAbin, posBbin;
    posAbin =    ( (posA[0] > 0)? 32 : 0 )
          + ( (posA[0] < 0)? 16 : 0 )
          + ( (posA[1] > 0)? 8 : 0 )
          + ( (posA[1] < 0)? 4 : 0 )
          + ( (posA[2] > 0)? 2 : 0 )
          + ( (posA[2] < 0)? 1 : 0 );

    posBbin =    ( (posB[0] > 0)? 32 : 0 )
          + ( (posB[0] < 0)? 16 : 0 )
          + ( (posB[1] > 0)? 8 : 0 )
          + ( (posB[1] < 0)? 4 : 0 )
          + ( (posB[2] > 0)? 2 : 0 )
          + ( (posB[2] < 0)? 1 : 0 );

    //if the intersection is not a segment, the intersection is not computed
    //the triangles intersects on a point (one vertex on the plane and the two others under or above
    if(posAbin == 5 || posAbin == 10 || posAbin == 17 || posAbin == 34 || posAbin == 20 || posAbin == 40
    || posBbin == 5 || posBbin == 10 || posBbin == 17 || posBbin == 34 || posBbin == 20 || posBbin == 40) return;
    //no possible intersection (one of the triangle is completely under or above the other
    if(posAbin == 42 || posAbin == 21
    || posBbin == 42 || posBbin == 21) return;
    //the triangles are coplanar
    if(posAbin == 0) return;

    //if an edge of a triangle is on the plane of the other triangle, it is necessary to verify if the
    //two polyhedra are intersecting on these edges, or if it only is a contact to know if the intersection
    //between the triangles must be computed or not.
    //"edgeA" and "edgeB" are codes
    //0 : the first edge is on the plane
    //1 : the second edge is on the plane
    //2 : the third edge is on the plane
    //3 : there is no edge on the plane
    unsigned short edgeA = 3, edgeB = 3;
    if(     posAbin == 1  || posAbin == 2 ) edgeA = 1; //points 0 and 1 on the plane
    else if(posAbin == 16 || posAbin == 32) edgeA = 2; //points 1 and 2 on the plane
    else if(posAbin == 4  || posAbin == 8 ) edgeA = 0; //points 2 and 0 on the plane
    if(     posBbin == 1  || posBbin == 2 ) edgeB = 1; //points 0 and 1 on the plane
    else if(posBbin == 16 || posBbin == 32) edgeB = 2; //points 1 and 2 on the plane
    else if(posBbin == 4  || posBbin == 8 ) edgeB = 0; //points 2 and 0 on the plane

    Vector_exact nA2, nB2;
    num_type p;
    bool invert_direction = false;
    bool stop = false;

    //if an edge of the first triangle is on the plane
    if(edgeA != 3 && edgeB == 3)
    {
      fA2 = heA[edgeA]->opposite()->facet();
      //TODO-elo-rm  nA2 = m_Inter_tri[fA2->Label].norm_dir;
      nA2 = m_Inter_tri[get(m_face_Label_A, fA2)].norm_dir;
      p = CGAL::cross_product(nA, nB) * CGAL::cross_product(nA2, nB);
      //if p is negative, the two triangles of the first polyhedron (including edgeA) are on the same side
      //so there is no intersection
      if(p < 0) stop = true;
      //if p == 0, fA2 is coplanar with the plane of fB
      //in that case, it is necessary to consider the boolean
      //operator used to determine if there is a contact or not
      else if(p == 0)
      {
        switch(m_BOOP)
        {
        case UNION:
          if(posA[(edgeA+1)%3] * (nA2 * nB) > 0) stop = true;
          break;
        case INTER:
          if(posA[(edgeA+1)%3] > 0) stop = true;
          break;
        case MINUS:
          if(posA[(edgeA+1)%3] * (nA2 * nB) < 0) stop = true;
          break;
        }
      }
      //the intersection between fA2 and fB is the same so this couple is removed from the list
      //TODO-elo-rm  rmCouple(fA2->Label, fB->Label);
      rmCouple(get(m_face_Label_A, fA2), get(m_face_Label_B, fB));
    }
    //if an edge of the second triangle is on the plane
    else if(edgeA == 3 && edgeB != 3)
    {
      fB2 = heB[edgeB]->opposite()->facet();
      //TODO-elo-rm  nB2 = m_Inter_tri[fB2->Label].norm_dir;
      nB2 = m_Inter_tri[get(m_face_Label_B, fB2)].norm_dir;
      p = CGAL::cross_product(nA, nB) * CGAL::cross_product(nA, nB2);
      //if p is negative, the two triangles of the second polyhedron (including edgeB) are on the same side
      //so there is no intersection
      if(p < 0) stop = true;
      //if p == 0, fB2 is coplanar with the plane of fA
      //in that case, it is necessary to consider the boolean
      //operator used to determine if there is a contact or not
      else if(p == 0)
      {
        switch(m_BOOP)
        {
        case UNION:
          if(posB[(edgeB+1)%3] < 0) stop = true;
          break;
        case INTER:
          if(posB[(edgeB+1)%3] * (nB2 * nA) < 0) stop = true;
          break;
        case MINUS:
          if(posB[(edgeB+1)%3] > 0) stop = true;
          break;
        }
      }
      //the intersection between fA and fB2 is the same so this couple is removed from the list
      //TODO-elo-rm  rmCouple(fA->Label, fB2->Label);
      rmCouple(get(m_face_Label_A, fA), get(m_face_Label_B, fB2));
    }
    //if an edge of each triangle is on the plane of the other
    else if(edgeA != 3 && edgeB != 3)
    {
      //in this case, four triangles are concerned by the intersection
      //fA2 and fB2 are the two other concerned facets
      //we try to determine if fA and fA2 are inside or outside the second polyhedron, using fB and fB2
      bool Intersection = false;
      Vector_exact nAcnB2, nA2cnB;
      num_type nAnB2, nA2nB, nA2nB2;
      num_type posA2_A, posB_A, posB2_A, posB_B2, posA_B, posB2_B, posB_A2, posB2_A2, posA2_B, posA2_B2;
      Point3d_exact ptA2, ptB2;

      fA2 = heA[edgeA]->opposite()->facet();
      fB2 = heB[edgeB]->opposite()->facet();
      //TODO-elo-rm  nA2 = m_Inter_tri[fA2->Label].norm_dir;
      nA2 = m_Inter_tri[get(m_face_Label_A, fA2)].norm_dir;
      //TODO-elo-rm  nB2 = m_Inter_tri[fB2->Label].norm_dir;
      nB2 = m_Inter_tri[get(m_face_Label_B, fB2)].norm_dir;

      nAcnB2 = CGAL::cross_product(nA, nB2);
      nA2cnB = CGAL::cross_product(nA2, nB);

      nAnB2 = nA * nB2;
      nA2nB = nA2 * nB;
      nA2nB2 = nA2 * nB2;

      ptA2 = point_to_exact(heA[edgeA]->opposite()->next()->vertex()->point());
      ptB2 = point_to_exact(heB[edgeB]->opposite()->next()->vertex()->point());

      posA_B = posA[(edgeA+1)%3];
      posB_A = posB[(edgeB+1)%3];
      posB_A2 = nA2 * (ptB[(edgeB+1)%3] - ptA[edgeA]);
      posB_B2 = nB2 * (ptB[(edgeB+1)%3] - ptA[edgeA]);
      posA2_A = nA * (ptA2 - ptA[edgeA]);
      posA2_B = nB * (ptA2 - ptA[edgeA]);
      posA2_B2 = nB2 * (ptA2 - ptA[edgeA]);
      posB2_A = nA * (ptB2 - ptA[edgeA]);
      posB2_A2 = nA2 * (ptB2 - ptA[edgeA]);
      posB2_B = nB * (ptB2 - ptA[edgeA]);

      if(nAcnB2 == CGAL::NULL_VECTOR && nA2cnB == CGAL::NULL_VECTOR
        && nAnB2 * nA2nB > 0) stop = true;

      //firstly, we search the position of fA
      //if fA is inside the poyhedron, Intersection = true
      if(posB_A * posB2_A > 0) //fB and fB2 on the same side
      {
        if(posB_B2 > 0) Intersection = true;
      }
      else if(posB_A * posB2_A < 0) //fB and fB2 on opposite side
      {
        if(posA_B < 0) Intersection = true;
      }
      else  //fA and fB2 coplanar
      {
        if(posA_B * posB2_B < 0)
        {
          if(posB_B2 > 0) Intersection = true;
        }
        else
        {
          if(nAnB2 < 0)
          {
            if(m_BOOP == UNION) Intersection = true;
          }
          else
          {
            if(m_BOOP == MINUS) Intersection = true;
          }
        }
      }

      //secondly, we search the position of fA2
      //if fA2 is inside the poyhedron, "Intersection" is inverted
      if(posB_A2 * posB2_A2 > 0) //fB and fB2 on the same side
      {
        if(posB_B2 > 0) Intersection = !Intersection;
      }
      else if(posB_A2 * posB2_A2 < 0) //fB and fB2 on opposite side
      {
        if(posA2_B < 0) Intersection = !Intersection;
      }
      else if(posB2_A2 == 0) //fA2 and fB2 coplanar
      {
        if(posA2_B * posB2_B < 0)
        {
          if(posB_B2 > 0) Intersection = !Intersection;
        }
        else
        {
          if(nA2nB2 < 0)
          {
            if(m_BOOP == UNION) Intersection = !Intersection;
          }
          else
          {
            if(m_BOOP == MINUS) Intersection = !Intersection;
          }
        }
      }
      else //fA2 and fB coplanar
      {
        if(posA2_B2 * posB_B2 < 0)
        {
          if(posB_B2 > 0) Intersection = !Intersection;
        }
        else
        {
          if(nA2nB < 0)
          {
            if(m_BOOP == UNION) Intersection = !Intersection;
          }
          else
          {
            if(m_BOOP == MINUS) Intersection = !Intersection;
          }
        }
      }

      //if Intersection == false, fA and fA2 are both inside or outside the second polyhedron.
      if(!Intersection) stop = true;

      //the intersection between (fA, fB2), (fA2, fB) and (fA2, fB2) are the same so these couples are removed from the list
      //TODO-elo-rm  rmCouple(fA->Label, fB2->Label);
      rmCouple(get(m_face_Label_A, fA), get(m_face_Label_B, fB2));
      //TODO-elo-rm  rmCouple(fA2->Label, fB->Label);
      rmCouple(get(m_face_Label_A, fA2), get(m_face_Label_B, fB));
      //TODO-elo-rm  rmCouple(fA2->Label, fB2->Label);
      rmCouple(get(m_face_Label_A, fA2), get(m_face_Label_B, fB2));

      //it is possible that the direction of the intersection have to be inverted
      if(posB_A * posA2_A > 0 && posB_A * posB2_A >= 0 && posB2_B * posA_B > 0) invert_direction = true;
    }

    //if the intersection must not be compute
    if(stop) return;

    Info_Inter inter[4];
    inter[0].f = fA;
    inter[1].f = fA;
    inter[2].f = fB;
    inter[3].f = fB;
    //TODO-elo-note:
    //  face of inter[0] and inter[1] are from mesh A
    //  face of inter[2] and inter[3] are from mesh B

    //the two intersection points between the edges of a triangle and the
    //other triangle are computed for the two triangles
    switch(posBbin)
    {
    //common intersections : one point one one side of the plane and the two other points on the other side
    //TODO-elo-note: inside this switch,
    //               inter[0].f and inter[1].f are from mesh A
    //               and inter[0].he and inter[1].he are from mesh B
    case 26:
    case 37:
      inter[0].he = heB[0];
      InterTriangleSegment(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[1];
      InterTriangleSegment(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    case 25:
    case 38:
      inter[0].he = heB[1];
      InterTriangleSegment(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[2];
      InterTriangleSegment(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    case 22:
    case 41:
      inter[0].he = heB[2];
      InterTriangleSegment(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[0];
      InterTriangleSegment(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    //particular cases : one point on the plane, one point one one side and one point on the other side
    case 6:
    case 9:
      inter[0].he = heB[2];
      InterTriangleSegment(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[0];
      IsInTriangle(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    case 18:
    case 33:
      inter[0].he = heB[0];
      InterTriangleSegment(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[1];
      IsInTriangle(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    case 24:
    case 36:
      inter[0].he = heB[1];
      InterTriangleSegment(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[2];
      IsInTriangle(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    //particular case : two points on the plane
    case 1:
    case 2:
      inter[0].he = heB[0];
      IsInTriangle(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[2]->opposite();
      IsInTriangle(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    case 16:
    case 32:
      inter[0].he = heB[1];
      IsInTriangle(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[0]->opposite();
      IsInTriangle(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    case 4:
    case 8:
      inter[0].he = heB[2];
      IsInTriangle(&inter[0], m_face_Label_A, m_halfedge_Label_B);
      inter[1].he = heB[1]->opposite();
      IsInTriangle(&inter[1], m_face_Label_A, m_halfedge_Label_B);
      break;
    default:
      return;
    }

    switch(posAbin)
    {
    //common intersections : one point one one side of the plane and the two other points on the other side
    //TODO-elo-note: inside this switch,
    //               inter[2].f and inter[3].f are from mesh B
    //               and inter[2].he and inter[3].he are from mesh A
    case 26:
    case 37:
      inter[2].he = heA[0];
      InterTriangleSegment(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[1];
      InterTriangleSegment(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    case 25:
    case 38:
      inter[2].he = heA[1];
      InterTriangleSegment(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[2];
      InterTriangleSegment(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    case 22:
    case 41:
      inter[2].he = heA[2];
      InterTriangleSegment(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[0];
      InterTriangleSegment(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    //particular cases : one point on the plane, one point one one side and one point on the other side
    case 6:
    case 9:
      inter[2].he = heA[2];
      InterTriangleSegment(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[0];
      IsInTriangle(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    case 18:
    case 33:
      inter[2].he = heA[0];
      InterTriangleSegment(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[1];
      IsInTriangle(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    case 24:
    case 36:
      inter[2].he = heA[1];
      InterTriangleSegment(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[2];
      IsInTriangle(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    //particular case : two points on the plane
    case 1:
    case 2:
      inter[2].he = heA[0];
      IsInTriangle(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[2]->opposite();
      IsInTriangle(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    case 16:
    case 32:
      inter[2].he = heA[1];
      IsInTriangle(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[0]->opposite();
      IsInTriangle(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    case 4:
    case 8:
      inter[2].he = heA[2];
      IsInTriangle(&inter[2], m_face_Label_B, m_halfedge_Label_A);
      inter[3].he = heA[1]->opposite();
      IsInTriangle(&inter[3], m_face_Label_B, m_halfedge_Label_A);
      break;
    default:
      return;
    }

    //if two distincts points belongs to the two triangles
    if(IsSegment(inter))
    {
      //we get this segment in ptInter
      std::vector<InterId> ptInter;
      Get_Segment(inter, ptInter);
      //and we build the opposite segment in ptInterInv
      std::vector<InterId> ptInterInv;
      ptInterInv.push_back(ptInter[1]);
      ptInterInv.push_back(ptInter[0]);

      //the segments are stored in the concerned triangles, and oriented
      if(CGAL::cross_product(nA, nB) * (m_InterPts[ptInter[1]] - m_InterPts[ptInter[0]]) * ((invert_direction == true)?-1:1) > 0)
      {
        switch(m_BOOP)
        {
        case UNION:
          //TODO-elo-rm  m_Inter_tri[fA->Label].CutList.push_back(ptInter);
          m_Inter_tri[get(m_face_Label_A, fA)].CutList.push_back(ptInter);
          //TODO-elo-rm  if(edgeA != 3) m_Inter_tri[fA2->Label].CutList.push_back(ptInter);
          if(edgeA != 3) m_Inter_tri[get(m_face_Label_A, fA2)].CutList.push_back(ptInter);
          //TODO-elo-rm  m_Inter_tri[fB->Label].CutList.push_back(ptInterInv);
          m_Inter_tri[get(m_face_Label_B, fB)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  if(edgeB != 3) m_Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
          if(edgeB != 3) m_Inter_tri[get(m_face_Label_B, fB2)].CutList.push_back(ptInterInv);
          break;
        case INTER:
          //TODO-elo-rm  m_Inter_tri[fA->Label].CutList.push_back(ptInterInv);
          m_Inter_tri[get(m_face_Label_A, fA)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  if(edgeA != 3) m_Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
          if(edgeA != 3) m_Inter_tri[get(m_face_Label_A, fA2)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  m_Inter_tri[fB->Label].CutList.push_back(ptInter);
          m_Inter_tri[get(m_face_Label_B, fB)].CutList.push_back(ptInter);
          //TODO-elo-rm  if(edgeB != 3) m_Inter_tri[fB2->Label].CutList.push_back(ptInter);
          if(edgeB != 3) m_Inter_tri[get(m_face_Label_B, fB2)].CutList.push_back(ptInter);
          break;
        case MINUS:
          //TODO-elo-rm  m_Inter_tri[fA->Label].CutList.push_back(ptInter);
          m_Inter_tri[get(m_face_Label_A, fA)].CutList.push_back(ptInter);
          //TODO-elo-rm  if(edgeA != 3) m_Inter_tri[fA2->Label].CutList.push_back(ptInter);
          if(edgeA != 3) m_Inter_tri[get(m_face_Label_A, fA2)].CutList.push_back(ptInter);
          //TODO-elo-rm  m_Inter_tri[fB->Label].CutList.push_back(ptInter);
          m_Inter_tri[get(m_face_Label_B, fB)].CutList.push_back(ptInter);
          //TODO-elo-rm  if(edgeB != 3) m_Inter_tri[fB2->Label].CutList.push_back(ptInter);
          if(edgeB != 3) m_Inter_tri[get(m_face_Label_B, fB2)].CutList.push_back(ptInter);
          break;
        }
      }
      else
      {
        switch(m_BOOP)
        {
        case UNION:
          //TODO-elo-rm  m_Inter_tri[fA->Label].CutList.push_back(ptInterInv);
          m_Inter_tri[get(m_face_Label_A, fA)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  m_Inter_tri[fA->Label].CutList.push_back(ptInterInv);
          m_Inter_tri[get(m_face_Label_A, fA)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  if(edgeA != 3) m_Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
          if(edgeA != 3) m_Inter_tri[get(m_face_Label_A, fA2)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  m_Inter_tri[fB->Label].CutList.push_back(ptInter);
          m_Inter_tri[get(m_face_Label_B, fB)].CutList.push_back(ptInter);
          //TODO-elo-rm  if(edgeB != 3) m_Inter_tri[fB2->Label].CutList.push_back(ptInter);
          if(edgeB != 3) m_Inter_tri[get(m_face_Label_B, fB2)].CutList.push_back(ptInter);
          break;
        case INTER:
          //TODO-elo-rm  m_Inter_tri[fA->Label].CutList.push_back(ptInter);
          m_Inter_tri[get(m_face_Label_A, fA)].CutList.push_back(ptInter);
          //TODO-elo-rm  if(edgeA != 3) m_Inter_tri[fA2->Label].CutList.push_back(ptInter);
          if(edgeA != 3) m_Inter_tri[get(m_face_Label_A, fA2)].CutList.push_back(ptInter);
          //TODO-elo-rm  m_Inter_tri[fB->Label].CutList.push_back(ptInterInv);
          m_Inter_tri[get(m_face_Label_B, fB)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  if(edgeB != 3) m_Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
          if(edgeB != 3) m_Inter_tri[get(m_face_Label_B, fB2)].CutList.push_back(ptInterInv);
          break;
        case MINUS:
          //TODO-elo-rm  m_Inter_tri[fA->Label].CutList.push_back(ptInterInv);
          m_Inter_tri[get(m_face_Label_A, fA)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  if(edgeA != 3) m_Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
          if(edgeA != 3) m_Inter_tri[get(m_face_Label_A, fA2)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  m_Inter_tri[fB->Label].CutList.push_back(ptInterInv);
          m_Inter_tri[get(m_face_Label_B, fB)].CutList.push_back(ptInterInv);
          //TODO-elo-rm  if(edgeB != 3) m_Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
          if(edgeB != 3) m_Inter_tri[get(m_face_Label_B, fB2)].CutList.push_back(ptInterInv);
          break;
        }
      }
    }
  }

  /*! \brief Compute the intersection between a facet and a halfedge
   * \param inter : A pointer to an Info_Inter structure.*/
  void InterTriangleSegment(Info_Inter* inter, const typename FEVV::Face_pmap< HalfedgeGraph, FacetId > &face_Label, const typename FEVV::Halfedge_pmap< HalfedgeGraph, HalfedgeId > &halfedge_Label)
  {
    Facet_handle f = inter->f;
    Halfedge_handle he = inter->he;
    auto f_Label = get(face_Label, f); //ELO+
    auto he_Label = get(halfedge_Label, he); //ELO+
    //if the intersection has been computed, the function returns directly the Id of the intersection
    //TODO-elo-rm  if(m_Inter_tri[f->Label].RefInter.count(he->Label) != 0)
    if(m_Inter_tri[f_Label].RefInter.count(he_Label) != 0)
    {
      //TODO-elo-rm  inter->Id = m_Inter_tri[f->Label].RefInter[he->Label];
      inter->Id = m_Inter_tri[f_Label].RefInter[he_Label];
      return;
    }
    //else, the calculation is done

    //this method is called when the intersection is not on the vertex pointed by the halfedge
    inter->IsOnVertex = false;
    //the intersection does not have an Id. 0xFFFFFFFF is set (this value means "no Id")
    inter->Id = 0xFFFFFFFF;

    Vector_exact e1, e2, dir, p, s, q;
    num_type u, v, tmp;

    Point3d_exact s1 = point_to_exact(he->opposite()->vertex()->point());
    Point3d_exact s2 = point_to_exact(he->vertex()->point());
    Point3d_exact v0 = point_to_exact(f->facet_begin()->vertex()->point());
    Point3d_exact v1 = point_to_exact(f->facet_begin()->next()->vertex()->point());
    Point3d_exact v2 = point_to_exact(f->facet_begin()->next()->next()->vertex()->point());

    //computation of the intersection (exact numbers)
    e1 = v1 - v0;
    e2 = v2 - v0;
    dir = s2 - s1;
    p = CGAL::cross_product(dir, e2);
    tmp = (num_type)1/(p*e1);
    s = s1 - v0;
    u = tmp * s * p;
    if(u < 0 || u > 1)
    {
      //the intersection is not in the triangle
      inter->res = 7;
      return;
    }
    q = CGAL::cross_product(s, e1);
    v = tmp * dir * q;
    if(v < 0 || v > 1)
    {
      //the intersection is not in the triangle
      inter->res = 7;
      return;
    }
    if(u + v > 1)
    {
      //the intersection is not in the triangle
      inter->res = 7;
      return;
    }

    //the result is stored in inter->pt
    inter->pt = s1+(tmp*e2*q)*dir;

    //creation of the code for the location of the intersection
    inter->res = 0;
    if(u == 0) inter->res += 1;  //intersection on he(0)
    if(v == 0) inter->res += 2;  //intersection on he(1)
    if(u+v == 1) inter->res += 4;  //intersection on he(2)
  }


  /*! \brief Finds the position of a point in a 3d triangle
   * \param inter : A pointer to an Info_Inter structure*/
  void IsInTriangle(Info_Inter* inter, const typename FEVV::Face_pmap< HalfedgeGraph, FacetId > &face_Label, const typename FEVV::Halfedge_pmap< HalfedgeGraph, HalfedgeId > &halfedge_Label)
  {
    Facet_handle f = inter->f;
    Halfedge_handle he = inter->he;
    auto f_Label = get(face_Label, f); //ELO+
    auto he_Label = get(halfedge_Label, he); //ELO+
    //if the intersection has been computed, the function returns directly the Id of the intersection
    //TODO-elo-rm  if(m_Inter_tri[f->Label].RefInter.count(he->Label) != 0)
    if(m_Inter_tri[f_Label].RefInter.count(he_Label) != 0)
    {
      //TODO-elo-rm  inter->Id = m_Inter_tri[f->Label].RefInter[he->Label];
      inter->Id = m_Inter_tri[f_Label].RefInter[he_Label];
      return;
    }
    //else, the calculation is done

    //this method is called when the intersection is exactly on the vertex pointed by the halfedge
    inter->IsOnVertex = true;
    //the intersection does not have an Id. 0xFFFFFFFF is set (this value means "no Id")
    inter->Id = 0xFFFFFFFF;

    Point3d_exact p = point_to_exact(he->vertex()->point());
    Point3d_exact v0 = point_to_exact(f->facet_begin()->vertex()->point());
    Point3d_exact v1 = point_to_exact(f->facet_begin()->next()->vertex()->point());
    Point3d_exact v2 = point_to_exact(f->facet_begin()->next()->next()->vertex()->point());

    //TODO-elo-rm  Vector_exact N = m_Inter_tri[f->Label].norm_dir;
    Vector_exact N = m_Inter_tri[f_Label].norm_dir;
    num_type u, v, w;

    u = N * CGAL::cross_product(v0 - v2, p - v2);
    if(u < 0)
    {
      //the intersection is not in the triangle
      inter->res = 7;
      return;
    }
    v = N * CGAL::cross_product(v1 - v0, p - v0);
    if(v < 0)
    {
      //the intersection is not in the triangle
      inter->res = 7;
      return;
    }
    w = N * CGAL::cross_product(v2 - v1, p - v1);
    if(w < 0)
    {
      //the intersection is not in the triangle
      inter->res = 7;
      return;
    }

    //the point is in the triangle
    inter->pt = p;

    //creation of the code for the location of the intersection
    inter->res = 0;
    if(u == 0) inter->res += 1;  //intersection on he(0)
    if(v == 0) inter->res += 2;  //intersection on he(1)
    if(w == 0) inter->res += 4;  //intersection on he(2)
  }

  /*! \brief Verify that the intersection is a segment
   * \param inter : A pointer to four Info_Inter structures
   * \return true if two distinct points are found in the four intersections computed*/
  bool IsSegment(Info_Inter *inter)
  {
    bool point = false; //true if a point is founded
    Point3d_exact pt; //the point founded
    bool id = false; //true if an Id is founded
    unsigned long Id = 0; // the Id founded // MT

    //each intersection is checked separately.
    //first intersection
    if(inter[0].Id != 0xFFFFFFFF)
    {
      //an Id different than 0xFFFFFFFF is founded
      //this intersection has already been computed and is valid
      id = true;
      Id = inter[0].Id;
    }
    else if(inter[0].res != 7)
    {
      //the intersection have no Id (0xFFFFFFFF)
      //but the intersection is in the triangle
      point = true;
      pt = inter[0].pt;
    }
    //second intersection
    if(inter[1].Id != 0xFFFFFFFF)
    {
      //an Id different than 0xFFFFFFFF is founded
      //this intersection has already been computed and is valid

      //if a point or an Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
      //(it is not possible that the two first points are the same)
      if(point || id) return true;
      id = true;
      Id = inter[1].Id;
    }
    else if(inter[1].res != 7)
    {
      //the intersection have no Id (0xFFFFFFFF)
      //but the intersection is in the triangle

      //if a point or an Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
      //(it is not possible that the two first points are the same)
      if(point || id) return true;
      point = true;
      pt = inter[1].pt;
    }
    //third intersection
    if(inter[2].Id != 0xFFFFFFFF)
    {
      //an Id different than 0xFFFFFFFF is founded
      //this intersection has already been computed and is valid

      //if a point or a different Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
      //(it is not possible that the two first points are the same)
      if(point || (id && Id != inter[2].Id)) return true;
      id = true;
      Id = inter[2].Id;
    }
    else if(inter[2].res != 7)
    {
      //the intersection have no Id (0xFFFFFFFF)
      //but the intersection is in the triangle

      //if an Id or a different point has already be founded, we founded the two distinct valid points (the intersection is a segment)
      //(it is not possible that the two first points are the same)
      if((point && pt != inter[2].pt) || id) return true;
      point = true;
      pt = inter[2].pt;
    }
    //fourth intersection
    if(inter[3].Id != 0xFFFFFFFF)
    {
      //an Id different than 0xFFFFFFFF is founded
      //this intersection has already been computed and is valid

      //if a point or a different Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
      //(it is not possible that the two first points are the same)
      if(point || (id && Id != inter[3].Id)) return true;
    }
    else if(inter[3].res != 7)
    {
      //the intersection have no Id (0xFFFFFFFF)
      //but the intersection is in the triangle

      //if an Id or a different point has already be founded, we founded the two distinct valid points (the intersection is a segment)
      //(it is not possible that the two first points are the same)
      if((point && pt != inter[3].pt) || id) return true;
    }
    return false;
  }

  /*! \brief Extracts the segment from a set of four intersection points and store these points in the list of intersecion points
   * \n There must be two valid and distinct points in the set
   * \param inter : A pointer to four Info_Inter structure
   * \param I : A vector to store the Id of the two intersection points (the output segment)*/
  void Get_Segment(Info_Inter *inter, std::vector<InterId> &I)
  {
    for(unsigned int i = 0;i != 4;++i)
    {
      //if the point have an Id
      if(inter[i].Id != 0xFFFFFFFF)
      {
        //the Id is stored if it is not already done
        if(I.size() == 0 || I[0] != inter[i].Id) I.push_back(inter[i].Id);
      }
      //else if the point is valid
      else if(inter[i].res != 7)
      {
        //the intersection point is stored in the list of the intersection points
        //and its new Id is stored in the output segment
        if(I.size() == 0 || m_InterPts[I[0]] != inter[i].pt)
        {
          Store_Intersection(&inter[i]);
          I.push_back(inter[i].Id);
        }
      }
      //return if the two points are founded
      if(I.size() == 2) return;
    }
  }

  /*! \brief Store the intersection and memorize it for every couples of facet-halfedge
   * \param inter : A pointer to an Info_Inter structure*/
  void Store_Intersection(Info_Inter *inter)
  {
    Facet_handle f;
    Halfedge_handle he;
    f = inter->f;
    he = inter->he;
    InterId I;

    //store the point to the list of the intersections and store its new Id
    inter->Id = m_InterPts.size();
    I = inter->Id;
    m_InterPts.push_back(inter->pt);

    //add this point as a vertex of the result
    ppbuilder.add_vertex(point_to_double(inter->pt), inter->Id);

    //if the intersection is on the vertex pointed by the halfedge (he), we update the Id (Label) of this vertex
    if(inter->IsOnVertex) he->vertex()->Label = I;

    //the intersection is memorized for each possible couple of (facet, halfedge) concerned by the intersection
    //if the intersection is exactly on the vertex pointed by the halfedge (he), it is necessary to take account
    //of every halfedge pointing to this vertex
    switch(inter->res)
    {
    case 0: //intersection on the facet
      {
        if(!inter->IsOnVertex)
        {
          m_Inter_tri[f->Label].RefInter[he->Label] = I;
          m_Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
        }
        else
        {
          Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
          do {
            m_Inter_tri[f->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
            H_circ++;
          } while(H_circ != H_end);
        }
      }
      break;
    case 1: //Intersection on the first halfedge of the facet
      {
        m_Inter_tri[f->Label].PtList.insert(I);
        m_Inter_tri[f->facet_begin()->opposite()->facet()->Label].PtList.insert(I);
        if(!inter->IsOnVertex)
        {
          m_Inter_tri[he->facet()->Label].PtList.insert(I);
          m_Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);
          m_Inter_tri[f->Label].RefInter[he->Label] = I;
          m_Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
          m_Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[he->Label] = I;
          m_Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
          m_Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->Label] = I;
          m_Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
          m_Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->Label] = I;
          m_Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
        }
        else
        {
          Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
          do {
            m_Inter_tri[f->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
            H_circ++;
          } while(H_circ != H_end);
        }
      }
      break;
    case 2: //Intersection on the second halfedge of the facet
      {
        m_Inter_tri[f->Label].PtList.insert(I);
        m_Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].PtList.insert(I);
        if(!inter->IsOnVertex)
        {
          m_Inter_tri[he->facet()->Label].PtList.insert(I);
          m_Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);
          m_Inter_tri[f->Label].RefInter[he->Label] = I;
          m_Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
          m_Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[he->Label] = I;
          m_Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
          m_Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
          m_Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
          m_Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
          m_Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
        }
        else
        {
          Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
          do {
            m_Inter_tri[f->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
            H_circ++;
          } while(H_circ != H_end);
        }
      }
      break;
    case 3: //Intersection on the first and second halfedge of the facet (vertex pointed by the first halfedge)
      {
        //update the Id (Label) of the first vertex of the facet
        f->facet_begin()->vertex()->Label = I;
        if(!inter->IsOnVertex)
        {
          m_Inter_tri[he->facet()->Label].PtList.insert(I);
          m_Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);

          Halfedge_around_vertex_circulator  H_circ = f->facet_begin()->vertex_begin(),
                            H_end = f->facet_begin()->vertex_begin();
          do {
            m_Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
            m_Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            H_circ++;
          } while(H_circ != H_end);
        }
        else
        {
          Halfedge_around_vertex_circulator  H_circ = he->vertex_begin(),
                            H_end = he->vertex_begin();
          do {
            Halfedge_around_vertex_circulator  F_circ = f->facet_begin()->vertex_begin(),
                              F_end = f->facet_begin()->vertex_begin();
            do {
              m_Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
              m_Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
              m_Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
              m_Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
              F_circ++;
            } while(F_circ != F_end);
            H_circ++;
          } while(H_circ != H_end);
        }
      }
      break;
    case 4: //Intersection on the third halfedge of the facet
      {
        m_Inter_tri[f->Label].PtList.insert(I);
        m_Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].PtList.insert(I);
        if(!inter->IsOnVertex)
        {
          m_Inter_tri[he->facet()->Label].PtList.insert(I);
          m_Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);
          m_Inter_tri[f->Label].RefInter[he->Label] = I;
          m_Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
          m_Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[he->Label] = I;
          m_Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
          m_Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
          m_Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
          m_Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
          m_Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
        }
        else
        {
          Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
          do {
            m_Inter_tri[f->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
            H_circ++;
          } while(H_circ != H_end);
        }
      }
      break;
    case 5: //Intersection on the first and third halfedge of the facet (vertex pointed by the third halfedge)
      {
        //update the Id (Label) of the third vertex of the facet
        f->facet_begin()->next()->next()->vertex()->Label = I;
        if(!inter->IsOnVertex)
        {
          m_Inter_tri[he->facet()->Label].PtList.insert(I);
          m_Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);

          Halfedge_around_vertex_circulator  H_circ = f->facet_begin()->next()->next()->vertex_begin(),
                            H_end = f->facet_begin()->next()->next()->vertex_begin();
          do {
            m_Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
            m_Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            H_circ++;
          } while(H_circ != H_end);
        }
        else
        {
          Halfedge_around_vertex_circulator   H_circ = he->vertex_begin(),
                            H_end = he->vertex_begin();
          do {
            Halfedge_around_vertex_circulator   F_circ = f->facet_begin()->next()->next()->vertex_begin(),
                              F_end = f->facet_begin()->next()->next()->vertex_begin();
            do {
              m_Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
              m_Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
              m_Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
              m_Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
              F_circ++;
            } while(F_circ != F_end);
            H_circ++;
          } while(H_circ != H_end);
        }
      }
      break;
    case 6: //Intersection on the second and third halfedge of the facet (vertex pointed by the second halfedge)
      {
        //update the Id (Label) of the second vertex of the facet
        f->facet_begin()->next()->vertex()->Label = I;
        if(!inter->IsOnVertex)
        {
          m_Inter_tri[he->facet()->Label].PtList.insert(I);
          m_Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);

          Halfedge_around_vertex_circulator  H_circ = f->facet_begin()->next()->vertex_begin(),
                            H_end = f->facet_begin()->next()->vertex_begin();
          do {
            m_Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
            m_Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
            m_Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            m_Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
            m_Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
            H_circ++;
          } while(H_circ != H_end);
        }
        else
        {
          Halfedge_around_vertex_circulator  H_circ = he->vertex_begin(),
                            H_end = he->vertex_begin();
          do {
            Halfedge_around_vertex_circulator  F_circ = f->facet_begin()->next()->vertex_begin(),
                              F_end = f->facet_begin()->next()->vertex_begin();
            do {
              m_Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
              m_Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
              m_Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
              m_Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
              F_circ++;
            } while(F_circ != F_end);
            H_circ++;
          } while(H_circ != H_end);
        }
      }
      break;
    }
  }

  /*! \brief Add a facet to the result
   * \param pFacet : A handle to the facet to add
   * \param facet_from_A : must be true if the facet belongs to the first polyhedron*/
  void add_facet_to_solution(Facet_handle &pFacet, bool facet_from_A)
  {
    //if the facet contains an intersection point but no intersection segment, the facet must be triangulate before
    if(pFacet->Label < m_Inter_tri.size())
    {
      Triangle_Cut TriCut = m_Inter_tri[pFacet->Label];
      Halfedge_handle he = pFacet->facet_begin();
      //creation of the triangulation
      Triangulation<Exact_Kernel> T(he, TriCut.norm_dir);
      //add the intersection points to the triangulation
      for(std::set<InterId>::iterator i = TriCut.PtList.begin();i != TriCut.PtList.end();++i)
      {
                            T.add_new_pt(m_InterPts[*i], (unsigned long &)*i);    // MT: ajout cast
      }
      //get all the triangles of the triangulation
                        vector<vector<unsigned long> > Tri_set = T.get_all_triangles((m_BOOP == MINUS && !TriCut.Facet_from_A)?true:false);
      //add these triangles to the result
      ppbuilder.add_triangle(Tri_set, he);
    }
    else
    {
      //the facet is added to the result. If the facet belongs to the second polyhedron, and if the
      //Boolean operation is a Subtraction, it is necessary to invert the orientation of the facet.
      if(m_BOOP == MINUS && !facet_from_A) ppbuilder.add_triangle(pFacet, true);
      else ppbuilder.add_triangle(pFacet, false);
    }
    //the tag of the three neighboring facets is updated
    pFacet->facet_begin()->opposite()->facet()->IsExt = true;
    pFacet->facet_begin()->next()->opposite()->facet()->IsExt = true;
    pFacet->facet_begin()->next()->next()->opposite()->facet()->IsExt = true;
  }

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG

  /*! \brief Colors the facets of the input polyhedra
   * \n The intersected facets are red
   * \n The facets that belong to the result are in green
   * \n The facets that does not belong to the result are in blue*/
  void ColorType()
  {
    Facet_iterator pFacet =  NULL;
    for (pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
    {
      if(pFacet->Label < m_Inter_tri.size()) pFacet->color(1.0, 0.0 ,0.0);
      else if(pFacet->IsExt) pFacet->color(0.0, 1.0 ,0.0);
      else pFacet->color(0.0, 0.0 ,1.0);
    }
    for (pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
    {
      if(pFacet->Label < m_Inter_tri.size()) pFacet->color(1.0, 0.0 ,0.0);
      else if(pFacet->IsExt) pFacet->color(0.0, 1.0 ,0.0);
      else pFacet->color(0.0, 0.0 ,1.0);
    }
  }


  /*! \brief Writes a report containing the computation time of the diffrent parts of the algorithm
   * \param pMout : The result polyhedron*/
  void WriteData(PolyhedronPtr &pMout)
  {
    unsigned int N_IFA = 0;
    unsigned int N_FFA = 0;
    unsigned int N_LFFA = 0;
    unsigned int N_IFB = 0;
    unsigned int N_FFB = 0;
    unsigned int N_LFFB = 0;
    unsigned int N_CF;

    for(Facet_iterator pFacet = m_pA->facets_begin();pFacet != m_pA->facets_end();++pFacet)
    {
      if(pFacet->Label < m_Inter_tri.size()) N_IFA++;
      else if(pFacet->IsExt) N_FFA++;
      else N_LFFA++;
    }
    for(Facet_iterator pFacet = m_pB->facets_begin();pFacet != m_pB->facets_end();++pFacet)
    {
      if(pFacet->Label < m_Inter_tri.size()) N_IFB++;
      else if(pFacet->IsExt) N_FFB++;
      else N_LFFB++;
    }
    N_CF = pMout->size_of_facets() - N_FFA - N_FFB;

    ofstrtime << "Computation time :"                                << std::endl;
    ofstrtime                                            << std::endl;
    ofstrtime << "Algorithm"                                    << std::endl;
    ofstrtime << "   Initialization :             "      << tr(duration_Init)          << std::endl;
    ofstrtime << " + Finding the Intersections :  "      << tr(duration_FindCouples)        << std::endl;
    ofstrtime << " + Compute the Intersections :  "      << tr(duration_ComputeIntersections)  << std::endl;
    ofstrtime << " + Cut the Intersected Facets : "      << tr(duration_CutIntersectedFacets)  << std::endl;
    ofstrtime << " + Complete the result :        "      << tr(duration_PropagateFacets)      << std::endl;
    ofstrtime << " + Create the polyhedron :      "      << tr(duration_delegate)        << std::endl;
    ofstrtime << "---------------------------------------"                      << std::endl;
    ofstrtime << " Total :                        "      << tr(duration_total)          << std::endl;
    ofstrtime                                            << std::endl;
    ofstrtime                                            << std::endl;
    ofstrtime << "Details :"                                    << std::endl;
    ofstrtime                                            << std::endl;
    ofstrtime << "Polyedron A :"                                  << std::endl;
    ofstrtime << "Number of Facets :                   "  << m_pA->size_of_facets()        << std::endl;
    ofstrtime << "Number of Intersected Facets :       "  << N_IFA                << std::endl;
    ofstrtime << "Number of Facets not in the result : "  << N_LFFA                << std::endl;
    ofstrtime                                            << std::endl;
    ofstrtime << "Polyedron B :"                                  << std::endl;
    ofstrtime << "Number of Facets :                   "  << m_pB->size_of_facets()        << std::endl;
    ofstrtime << "Number of Intersected Facets :       "  << N_IFB                << std::endl;
    ofstrtime << "Number of Facets not in the result : "  << N_LFFB                << std::endl;
    ofstrtime                                            << std::endl;
    ofstrtime << "Result :"                                      << std::endl;
    ofstrtime << "Number of Facets :                   "  << pMout->size_of_facets()        << std::endl;
    ofstrtime << "Number of Facets from A :            "  << N_FFA                << std::endl;
    ofstrtime << "Number of Facets from B :            "  << N_FFB                << std::endl;
    ofstrtime << "Number of Created Facets :           "  << N_CF                  << std::endl;
  }

#endif // BOOLEAN_OPERATIONS_DEBUG
#endif //TODO-elo-WIP

  /*! \brief Boolean operation computed*/
  Bool_Op m_BOOP;

  /*! \brief The first input polyhedron*/
  HalfedgeGraph *m_pA; //TODO-elo-rm PolyhedronPtr m_pA;
  PointMap *m_pmA;
  /*! \brief The second input polyhedron*/
  HalfedgeGraph *m_pB; //TODO-elo-rm PolyhedronPtr m_pB;
  PointMap *m_pmB;

  /*! \brief  Property maps */
  typename FEVV::Vertex_pmap< HalfedgeGraph, VertexId >
      m_vertex_Label_A; // ELO replace pVertex->Label
  typename FEVV::Vertex_pmap< HalfedgeGraph, VertexId >
      m_vertex_Label_B; // ELO replace pVertex->Label

  typename FEVV::Halfedge_pmap< HalfedgeGraph, HalfedgeId >
      m_halfedge_Label_A; // ELO replace facet_begin()->Label
  typename FEVV::Halfedge_pmap< HalfedgeGraph, HalfedgeId >
      m_halfedge_Label_B; // ELO replace facet_begin()->Label

  typename FEVV::Face_pmap< HalfedgeGraph, FacetId >
      m_face_Label_A; // ELO replace pFacet->Label
  typename FEVV::Face_pmap< HalfedgeGraph, FacetId >
      m_face_Label_B; // ELO replace pFacet->Label
  typename FEVV::Face_pmap< HalfedgeGraph, bool >
      m_face_IsExt_A; // ELO replace pFacet->IsExt
  typename FEVV::Face_pmap< HalfedgeGraph, bool >
      m_face_IsExt_B; // ELO replace pFacet->IsExt
  typename FEVV::Face_pmap< HalfedgeGraph, bool >
      m_face_IsOK_A; // ELO replace pFacet->IsOK
  typename FEVV::Face_pmap< HalfedgeGraph, bool >
      m_face_IsOK_B; // ELO replace pFacet->IsOK

  /*! \brief The polyhedron builder*/
  //TODO-elo-rm  CPolyhedron_from_polygon_builder_3<HDS> ppbuilder;

  /*! \brief Lists the couples of facets that intersect*/
  std::map< FacetId, std::set< FacetId > > m_Couples;
  /*! \brief Lists the exact intersection points computed*/
  std::vector< Point3d_exact > m_InterPts;
  /*! \brief Informations about the intersected facets*/
  std::vector< Triangle_Cut > m_Inter_tri;
  /*! \brief Index to obtain the handle of a facet with its Id*/
  std::vector< Facet_handle > m_Facet_Handle;

  /*! \brief the AABB-tree*/
  AABB_Tree tree;

#if 0 //#ifdef BOOLEAN_OPERATIONS_DEBUG

  /*! \brief The output report file*/
  std::ofstream ofstrtime;
  /*! \brief Initialisation time*/
  double duration_Init;
  /*! \brief Time to find all the couples of facet that intersect*/
  double duration_FindCouples;
  /*! \brief Time to Compute all the intersections*/
  double duration_ComputeIntersections;
  /*! \brief Time to Cut the facets*/
  double duration_CutIntersectedFacets;
  /*! \brief Time to complete the result with the facets from the input polyhedra*/
  double duration_PropagateFacets;
  /*! \brief Time to create the result polyhedron*/
  double duration_delegate;
  /*! \brief Time to Compute a Boolean operation*/
  double duration_total;

#endif // BOOLEAN_OPERATIONS_DEBUG
};
