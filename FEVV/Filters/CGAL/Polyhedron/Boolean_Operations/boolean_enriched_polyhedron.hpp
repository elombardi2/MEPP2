//TODO-elo-licence
#pragma once

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_3.h>

template <class Refs, class T, class P, class Norm>
class MEPP_Common_Vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
	protected:
		// tag
		int m_tag;

		// normal
		Norm m_normal;

		// color
		float m_color[3];

		// texture coordinates
		float m_texture_coordinates[2];

	public:
		// life cycle
		MEPP_Common_Vertex()
		{
			color(0.5f, 0.5f, 0.5f);
			texture_coordinates(0.0f, 0.0f);
		}
		// repeat mandatory constructors
		MEPP_Common_Vertex(const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
		{
			color(0.5f, 0.5f, 0.5f);
			texture_coordinates(0.0f, 0.0f);
		}

		// color
		float color(int index) { return m_color[index]; }
		void color(float r, float g, float b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }

		// texture coordinates
		float texture_coordinates(int index) { return m_texture_coordinates[index]; }
		void texture_coordinates(float u, float v) { m_texture_coordinates[0] = u; m_texture_coordinates[1] = v; }

		// normal
		typedef Norm Normal_3;
  		//typedef Norm Vector;

		Normal_3& normal() { return m_normal; }
		const Normal_3& normal() const { return m_normal; }

		// tag
		int& tag() {  return m_tag; }
		const int& tag() const {  return m_tag; }
		void tag(const int& t)  { m_tag = t; }
};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class MEPP_Common_Halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
	protected:
		// tag
		int m_tag;

#if (0)
		// normal
		Norm m_normal;  // AJOUT C�line : halfedge normal (sum of 2 incident facet normals)
#endif

		// option for edge superimposing
		bool m_control_edge;

		// texture coordinates : AJOUT Laurent Chevalier
		float m_texture_coordinates[2];

	public:
		// life cycle
		MEPP_Common_Halfedge()
		{
			m_control_edge = true;

			// texture coordinates : AJOUT Laurent Chevalier
			texture_coordinates(0.0f, 0.0f);
		}

		// tag
		const int& tag() const { return m_tag;  }
		int& tag() { return m_tag;  }
		void tag(const int& t)  { m_tag = t; }

#if (0)
		// normal : AJOUT C�line
		typedef Norm Normal_3;
		Normal_3& normal() { return m_normal; }
		const Normal_3& normal() const { return m_normal; }
#endif

		// texture coordinates : AJOUT Laurent Chevalier
		float texture_coordinates(int index) { return m_texture_coordinates[index]; }
		void texture_coordinates(float u, float v) { m_texture_coordinates[0] = u; m_texture_coordinates[1] = v; }

		// control edge
		bool& control_edge()  { return m_control_edge; }
		const bool& control_edge()  const { return m_control_edge; }
		void control_edge(const bool& flag) { m_control_edge  = flag; }
};


template <class Refs, class T, class Norm>
class MEPP_Common_Facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
	protected:
		// tag
		int m_tag;

		// normal
		Norm m_normal;

		// color
		float m_color[3];

	public:
		// life cycle
		MEPP_Common_Facet()
		{
			color(0.5f, 0.5f, 0.5f);
		}

		// tag
		const int& tag() const { return m_tag;  }
		int& tag() { return m_tag;  }
		void tag(const int& t)  { m_tag = t; }

		// normal
		typedef Norm Normal_3;
		Normal_3& normal() { return m_normal; }
		const Normal_3& normal() const { return m_normal; }

		// color
		float color(int index) { return m_color[index]; }
		void color(float r, float g, float b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }
};


template <class Refs, class T, class P, class Norm>
class Enriched_vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:
		Enriched_vertex()
		{
		}

		// Le constructeur PAS par d�faut appelle les constructeurs PAR DEFAUT des classes ancetres
		// en plus de celles appelles appell�es explicitement.
		// On a besoin d'appeller explicitement le constructeur(pt) de base pour la cr�ation de polyhedre
		Enriched_vertex(const P& pt) : MEPP_Common_Vertex<Refs, T, P, Norm>(pt)
		{
			this->point() = pt;
		}

		// La creation du polyhedre implique un appel au constructeur par copie,
		// qui appelle le constructeur par d�faut de base et ne copie pas le point.
		// Il faut donc avoir un constructeur par copie explicite qui s'occupe du point.
		Enriched_vertex(const Enriched_vertex& v)
		{
			this->point() = v.point();
		}
};






struct Enriched_items : public CGAL::Polyhedron_items_3
{
    // wrap vertex
    template <class Refs, class Traits>
    struct Vertex_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_vertex<Refs,
                          CGAL::Tag_true,
                          Point,
                          Normal> Vertex;
    };

    // wrap face
    template <class Refs, class Traits>
    struct Face_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
		typedef typename Traits::Plane_3 Plane;
        typedef Enriched_facet<Refs,
                         CGAL::Tag_true,
                         Point,
                         Normal,
						 Plane> Face;
    };

    // wrap halfedge
    template <class Refs, class Traits>
    struct Halfedge_wrapper
    {
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_halfedge<Refs,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            Normal> Halfedge;
    };
};





using CGALKernel = CGAL::Cartesian< double >;
using MeshPolyhedron =
    CGAL::Polyhedron_3< CGALKernel, CGAL::Polyhedron_items_with_id_3 >;










typedef Enriched_polyhedron<Enriched_kernel, Enriched_items>	Polyhedron;