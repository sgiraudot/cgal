// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ 
// $Date$
// 
//
// Author(s)     : Aymeric PELLE <aymeric.pelle@sophia.inria.fr>

#include <CGAL/Periodic_3_Regular_triangulation_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Gmpz.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Periodic_3_Regular_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <set>


typedef CGAL::Epeck K;
typedef K::FT FT;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
typedef CGAL::Periodic_3_Regular_triangulation_traits_3<Regular_traits> Traits;

template class CGAL::Periodic_3_Regular_triangulation_3<Traits>;
typedef CGAL::Periodic_3_Regular_triangulation_3<Traits> P3RT3;
typedef P3RT3::Vertex_handle Vertex_handle;
typedef P3RT3::Cell_handle Cell_handle;

typedef Traits::Weighted_point Weighted_point;
typedef Traits::Bare_point Bare_point;
typedef Traits::Iso_cuboid_3 Iso_cuboid;


template <class PeriodicTriangulation>
void periodic_triangulation_to_medit_file (const PeriodicTriangulation& pt, std::ostream& stream)
{
//  typedef typename PeriodicTriangulation::Point Point;
  typedef typename PeriodicTriangulation::Triangle Triangle;
  typedef typename PeriodicTriangulation::Periodic_triangle_iterator Periodic_iterator;

  Periodic_iterator ps_b = pt.periodic_triangles_begin(PeriodicTriangulation::STORED);
  Periodic_iterator ps_e = pt.periodic_triangles_end(PeriodicTriangulation::STORED);
  std::size_t ps_dist = std::distance(ps_b, ps_e);

  stream << "MeshVersionFormatted 1\n"
            "Dimension 3\n"
            "Vertices\n"
         << (ps_dist * 3)
         << std::endl;

//  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
//  {
//    Segment s = pt.segment(*iter);
//    stream << s.source().x() << " " << s.source().y() << " " << s.source().z() << " 1" << '\n';
//    stream << s.target().x() << " " << s.target().y() << " " << s.target().z() << " 1" << '\n';
//  }

  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
  {
    Triangle t = pt.triangle(*iter);
    stream << t[0].x() << " " << t[0].y() << " " << t[0].z() << " 1\n";
    stream << t[1].x() << " " << t[1].y() << " " << t[1].z() << " 1\n";
    stream << t[2].x() << " " << t[2].y() << " " << t[2].z() << " 1\n";
  }
  unsigned count = 0;
  stream << "Triangles\n" << ps_dist << std::endl;
  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
  {
    stream << ++count;
    stream << " " << ++count;
    stream << " " << ++count;
    stream << " 2" << '\n';
  }
  stream.flush();
}

template <class PeriodicTriangulation>
void periodic_triangulation_to_medit_1_file (const PeriodicTriangulation& pt, std::ostream& stream)
{
//  typedef typename PeriodicTriangulation::Point Point;
  typedef typename PeriodicTriangulation::Triangle Triangle;
  typedef typename PeriodicTriangulation::Periodic_triangle_iterator Periodic_iterator;

  Periodic_iterator ps_b = pt.periodic_triangles_begin(PeriodicTriangulation::UNIQUE);
  Periodic_iterator ps_e = pt.periodic_triangles_end(PeriodicTriangulation::UNIQUE);
  std::size_t ps_dist = std::distance(ps_b, ps_e);

  stream << "MeshVersionFormatted 1\n"
            "Dimension 3\n"
            "Vertices\n"
         << (ps_dist * 3)
         << std::endl;

//  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
//  {
//    Segment s = pt.segment(*iter);
//    stream << s.source().x() << " " << s.source().y() << " " << s.source().z() << " 1" << '\n';
//    stream << s.target().x() << " " << s.target().y() << " " << s.target().z() << " 1" << '\n';
//  }

  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
  {
    Triangle t = pt.triangle(*iter);
    stream << t[0].x() << " " << t[0].y() << " " << t[0].z() << " 1\n";
    stream << t[1].x() << " " << t[1].y() << " " << t[1].z() << " 1\n";
    stream << t[2].x() << " " << t[2].y() << " " << t[2].z() << " 1\n";
  }
  unsigned count = 0;
  stream << "Triangles\n" << ps_dist << std::endl;
  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
  {
    stream << ++count;
    stream << " " << ++count;
    stream << " " << ++count;
    stream << " 2" << '\n';
  }
  stream.flush();
}

template <class PeriodicTriangulation>
void periodic_triangulation_to_medit_edges_file (const PeriodicTriangulation& pt, std::ostream& stream)
{
//  typedef typename PeriodicTriangulation::Point Point;
  typedef typename PeriodicTriangulation::Segment Segment;
  typedef typename PeriodicTriangulation::Periodic_segment_iterator Periodic_iterator;

  Periodic_iterator ps_b = pt.periodic_segments_begin(PeriodicTriangulation::STORED);
  Periodic_iterator ps_e = pt.periodic_segments_end(PeriodicTriangulation::STORED);
  std::size_t ps_dist = std::distance(ps_b, ps_e);

  stream << "MeshVersionFormatted 1\n"
            "Dimension 3\n"
            "Vertices\n"
         << (ps_dist * 2)
         << std::endl;

  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
  {
    Segment t = pt.segment(*iter);
    stream << t[1].x() << " " << t[1].y() << " " << t[1].z() << " 1\n";
    stream << t[2].x() << " " << t[2].y() << " " << t[2].z() << " 1\n";
  }
  unsigned count = 0;
  stream << "Edges\n" << ps_dist << std::endl;
  for (Periodic_iterator iter = ps_b, end_iter = ps_e; iter != end_iter; ++iter)
  {
    stream << ++count;
    stream << " " << ++count;
    stream << " 1" << '\n';
  }
  stream.flush();
}

Weighted_point read_wpoint (std::istream& stream)
{
  FT x = 0., y = 0., z = 0., w = 0.;
  stream >> x;
  assert(stream && !stream.eof());
  stream >> y;
  assert(stream && !stream.eof());
  stream >> z;
  assert(stream && !stream.eof());
  stream >> w;
  assert(stream);
  return Weighted_point(Bare_point(x, y, z), w);
}


void test_construction ()
{
  std::cout << "--- test_construction" << std::endl;

  P3RT3 p3rt3;
  assert(p3rt3.is_valid());
}

void test_insert_1 ()
{
  std::cout << "--- test_insert_1" << std::endl;

  P3RT3 p3rt3;

  Weighted_point p(Bare_point(0,0,0), 0.01);
  p3rt3.insert(p);

  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 1);
  assert(p3rt3.number_of_stored_vertices() == 27);

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
  assert(hidden_point_count == 0);
}

void test_insert_point ()
{
  std::cout << "--- test_insert_point" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  Vertex_handle vh;
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.9,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.9,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.9),0.01));
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  Weighted_point p(Bare_point(0.4, 0.4, 0.4), 0.001);
  vh = p3rt3.insert(p);
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 5);
  assert(p3rt3.number_of_stored_vertices() == 135);

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
  assert(hidden_point_count == 0);
}

void test_insert_hidden_point ()
{
  std::cout << "--- test_insert_hidden_point" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  Vertex_handle vh;
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.9,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.9,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.9),0.01));
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
  assert(hidden_point_count == 0);

  Weighted_point hidden_point(Bare_point(0.101, 0.101, 0.101), 0.001);
  vh = p3rt3.insert(hidden_point);
  assert(vh == Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  unsigned hidden_found_count = 0;
  hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
  {
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    hidden_found_count += static_cast<unsigned>(std::find(iter->hidden_points_begin(), iter->hidden_points_end(), hidden_point) != iter->hidden_points_end());
  }
  assert(hidden_point_count == 1);
  assert(hidden_found_count == 1);
}

void test_insert_hiding_point ()
{
  std::cout << "--- test_insert_hiding_point" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  Vertex_handle vh;
  vh = p3rt3.insert(Weighted_point(Bare_point(0.9,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.9,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.9),0.01));
  assert(vh != Vertex_handle());
  Weighted_point hidden_point(Bare_point(0.101, 0.101, 0.101), 0.001);
  vh = p3rt3.insert(hidden_point);
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
  assert(hidden_point_count == 0);

  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  unsigned hidden_found_count = 0;
  hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
  {
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    hidden_found_count += static_cast<unsigned>(std::find(iter->hidden_points_begin(), iter->hidden_points_end(), hidden_point) != iter->hidden_points_end());
  }
  assert(hidden_point_count == 1);
  assert(hidden_found_count == 1);
}

void test_insert_a_point_twice ()
{
  std::cout << "--- test_insert_a_point_twice" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  Vertex_handle vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 1);
  assert(p3rt3.number_of_stored_vertices() == 27);

  Vertex_handle vh2 = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.1),0.01));
  assert(vh2 == vh);
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 1);
  assert(p3rt3.number_of_stored_vertices() == 27);

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
  assert(hidden_point_count == 0);
}

void test_insert_two_points_with_the_same_position ()
{
  std::cout << "--- test_insert_two_points_with_the_same_position" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  Weighted_point hidden_point(Bare_point(0.1,0.1,0.1),0.01);
  Vertex_handle vh = p3rt3.insert(hidden_point);
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 1);
  assert(p3rt3.number_of_stored_vertices() == 27);

  Vertex_handle vh2 = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.1),0.015));
  assert(vh2 != Vertex_handle());
  assert(vh2 != vh);
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 1);
  assert(p3rt3.number_of_stored_vertices() == 27);

  unsigned hidden_found_count = 0;
  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
  {
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    hidden_found_count += static_cast<unsigned>(std::find(iter->hidden_points_begin(), iter->hidden_points_end(), hidden_point) != iter->hidden_points_end());
  }
  assert(hidden_point_count == 1);
  assert(hidden_found_count == 1);
}

void test_remove ()
{
  std::cout << "--- test_remove" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  Vertex_handle vh;
  vh = p3rt3.insert(Weighted_point(Bare_point(0.9,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.9,0.1),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.9),0.01));
  assert(vh != Vertex_handle());
  vh = p3rt3.insert(Weighted_point(Bare_point(0.1,0.1,0.1),0.01));
  assert(vh != Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  Weighted_point hidden_point(Bare_point(0.101, 0.101, 0.101), 0.001);
  Vertex_handle vhh = p3rt3.insert(hidden_point);
  assert(vhh == Vertex_handle());
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  p3rt3.remove(vh);
  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == 4);
  assert(p3rt3.number_of_stored_vertices() == 108);

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
  assert(hidden_point_count == 0);

  unsigned point_found_count = 0;
  for (P3RT3::Vertex_iterator iter = p3rt3.vertices_begin(), end_iter = p3rt3.vertices_end(); iter != end_iter; ++iter)
    point_found_count += (iter->point() == hidden_point);
  assert(point_found_count == 27);
}

void test_insert_rnd_as_delaunay (unsigned pt_count, double weight)
{
  std::cout << "--- test_insert_rnd_as_delaunay (" << pt_count << ',' << weight << ')' << std::endl;

  CGAL::Random random(7);
  typedef CGAL::Creator_uniform_3<double,Bare_point>  Creator;
  CGAL::Random_points_in_cube_3<Bare_point, Creator> in_cube(0.5, random);

  P3RT3::Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  P3RT3 p3rt3(iso_cuboid);

  for (unsigned cnt = 1; cnt <= pt_count; ++cnt)
  {
    Weighted_point p(*in_cube++, weight);
//    std::cout << cnt << " : " << p << std::endl;
    assert(iso_cuboid.has_on_bounded_side(p));
    assert(p.weight() < 0.015625);
    p3rt3.insert(p);
  }

  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_vertices() == pt_count);
  assert(p3rt3.number_of_sheets() == CGAL::make_array(3,3,3) ?
      p3rt3.number_of_stored_vertices() == 27 * pt_count
      : p3rt3.number_of_stored_vertices() == pt_count);

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
  assert(hidden_point_count == 0);
}

void test_insert_rnd_then_remove_all (unsigned pt_count, unsigned seed)
{
  std::cout << "--- test_insert_rnd (" << pt_count << ", " << seed << ')' << std::endl;

  CGAL::Random random(seed);
  typedef CGAL::Creator_uniform_3<double,Bare_point>  Creator;
  CGAL::Random_points_in_cube_3<Bare_point, Creator> in_cube(0.5, random);

  P3RT3::Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  P3RT3 p3rt3(iso_cuboid);

  std::ofstream stream("out_p3rt3_test");
  assert(stream);

  std::vector<Weighted_point> insert_set;
  insert_set.reserve(pt_count);
  std::vector<Weighted_point> remove_set;
  remove_set.reserve(pt_count);

  std::cout << "-- insert" << std::endl;
  for (unsigned cnt = 1; cnt <= pt_count; ++cnt)
  {
    Weighted_point p(*in_cube++, random.get_double(0., 0.015625));
//    std::cout << cnt << " : " << p << std::endl;
    assert(p.weight() < 0.015625);
    stream << p << std::endl;

    unsigned hidden_point_count = 0;
    for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());

    Vertex_handle vh = p3rt3.insert(p);

    unsigned hidden_point_count_2 = 0;
    for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count_2 += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count <= hidden_point_count_2);
    assert(hidden_point_count_2 + p3rt3.number_of_vertices() == cnt);

    std::cout << cnt << " - p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
    if (vh == Vertex_handle())
    {
      if (find(insert_set.begin(), insert_set.end(), p) == insert_set.end())
        insert_set.push_back(p);
    }
    else
      insert_set.push_back(p);
  }

  stream.close();

  assert(p3rt3.is_valid());

  std::cout << "-- remove" << std::endl;
  unsigned cnt = 1;
  for (; p3rt3.number_of_vertices() != 0; ++cnt)
  {
    P3RT3::Vertex_iterator iter = p3rt3.vertices_begin();
    for (unsigned j = random.get_int(0, p3rt3.number_of_vertices()); j; --j)
      ++iter;

    std::cout << cnt << " : " << iter->point() << std::endl;

    remove_set.push_back(iter->point());
    p3rt3.remove(iter);

    std::cout << "    p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;

    unsigned hidden_point_count = 0;
    for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());

    assert(hidden_point_count + cnt + p3rt3.number_of_vertices() == insert_set.size());
  }

  std::sort(insert_set.begin(), insert_set.end());
  std::sort(remove_set.begin(), remove_set.end());
  assert(insert_set == remove_set);

  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_stored_vertices() == 0);
}

void test_insert_from_file (const char* filename)
{
  std::cout << "--- test_insert_from_file" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5));

  std::ifstream stream(filename);
  assert(stream);

  unsigned cnt = 1;
  while (stream && !(stream.eof()))
  {
    Weighted_point p = read_wpoint(stream);
    std::cout << cnt << " : " << p << std::endl;
    assert(p.weight() < 0.015625);
    std::cout << "p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
    p3rt3.insert(p);
    std::cout << "p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
    ++cnt;
  }

  assert(p3rt3.is_valid());
}

void test_insert_rt3_pointset ()
{
  P3RT3 p3rt3(Iso_cuboid(Bare_point(-100,-100,-100),Bare_point(100,100,100)));

  for (int a=0;a!=10;a++)
    for (int b=0;b!=5;b++)
      for (int d=0;d!=5;d++)
      {
        Weighted_point p( Bare_point(a*b-d*a + (a-b)*10 +a , a-b+d +5*b, a*a-d*d+b),  a*b-a*d ); // TODO check weight
        std::cout << p << std::endl;
        p3rt3.insert(p);
      }
  assert(p3rt3.is_valid());
}

void test_27_to_1_sheeted_covering ()
{
  std::cout << "--- test_27_to_1_sheeted_covering" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0, 0, 0, 1, 1, 1));

  unsigned count = 1;
  for (unsigned i = 0; i < 6; ++i)
    for (unsigned j = 0; j < 6; ++j)
      for (unsigned k = 0; k < 8; ++k)
      {
        FT x = FT(i) / FT(6);
        if (k % 2)
          x += FT(1) / FT(12);
        FT y = FT(j) / FT(6);
        if (k % 2)
          y += FT(1) / FT(12);
        FT z = FT(k) / FT(8);
        std::cout << count++ << " - " << i << " " << j << " " << k << std::endl;
        Weighted_point point(Bare_point(x, y, z), 0);
        p3rt3.insert(point);
        if (CGAL::make_array(i,j,k) != CGAL::make_array<unsigned>(5,5,7))
        {
          assert(p3rt3.number_of_sheets() == CGAL::make_array(3,3,3));
        }
      }

  assert(p3rt3.number_of_sheets() == CGAL::make_array(1,1,1));
  assert(p3rt3.number_of_vertices() == 6*6*8);
  assert(p3rt3.is_valid());

  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());

  assert(hidden_point_count == 0);
}

void test_dummy_points ()
{
  std::cout << "--- test_dummy_points" << std::endl;

  P3RT3 p3rt3(P3RT3::Iso_cuboid(0, 0, 0, 1, 1, 1));

  std::vector<Vertex_handle> vertices;
  vertices.reserve(p3rt3.number_of_vertices());

  unsigned count = 1;
  for (unsigned i = 0; i < 6; ++i)
    for (unsigned j = 0; j < 6; ++j)
      for (unsigned k = 0; k < 8; ++k)
      {
        FT x = FT(i) / FT(6);
        if (k % 2)
          x += FT(1) / FT(12);
        FT y = FT(j) / FT(6);
        if (k % 2)
          y += FT(1) / FT(12);
        FT z = FT(k) / FT(8);
        std::cout << count++ << " - " << i << " " << j << " " << k << std::endl;
        Weighted_point point(Bare_point(x, y, z), 0);
        vertices.push_back(p3rt3.insert(point));
        if (CGAL::make_array(i,j,k) != CGAL::make_array<unsigned>(5,5,7))
        {
          assert(p3rt3.number_of_sheets() == CGAL::make_array(3,3,3));
        }
      }

  assert(p3rt3.number_of_sheets() == CGAL::make_array(1,1,1));
  assert(p3rt3.number_of_vertices() == 6*6*8);
  assert(p3rt3.is_valid());

  P3RT3 p3rt3_b(P3RT3::Iso_cuboid(0, 0, 0, 1, 1, 1));
  p3rt3_b.insert_dummy_points();

  assert(p3rt3_b.number_of_sheets() == CGAL::make_array(1,1,1));
  assert(p3rt3_b.number_of_vertices() == 6*6*8);
  assert(p3rt3_b.is_valid());

  assert(p3rt3 == p3rt3_b);
}

void test_insert_range (unsigned pt_count, unsigned seed)
{
  CGAL::Random random(seed);
  typedef CGAL::Creator_uniform_3<double,Bare_point>  Creator;
  CGAL::Random_points_in_cube_3<Bare_point, Creator> in_cube(0.5, random);

  P3RT3::Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  P3RT3 p3rt3(iso_cuboid);

  std::vector<Weighted_point> points;
  points.reserve(pt_count);

  while (points.size() != pt_count)
  {
    Weighted_point p(*in_cube++, random.get_double(0., 0.015625));
    points.push_back(p);
  }

  p3rt3.insert(points.begin(), points.end(), true);

  for (P3RT3::Vertex_iterator iter = p3rt3.vertices_begin(), end_iter = p3rt3.vertices_end(); iter != end_iter; ++iter)
  {
    std::vector<Weighted_point>::iterator it = std::find(points.begin(), points.end(), iter->point());
    assert(it != points.end());
  }
  unsigned hidden_point_count = 0;
  for (P3RT3::Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
  {
    for (P3RT3::Cell::Point_iterator it = iter->hidden_points_begin(), end_it = iter->hidden_points_end(); it != end_it; ++it)
    {
      assert(std::find(points.begin(), points.end(), *it) != points.end());
      ++hidden_point_count;
    }
  }
  assert(p3rt3.number_of_vertices() == 659);
  assert(p3rt3.number_of_vertices() + hidden_point_count == 800);

  assert(p3rt3.is_valid());
  assert(p3rt3.number_of_sheets() == CGAL::make_array(1,1,1));
}

int main (int argc, char** argv)
{
  std::cout << "TESTING ..." << std::endl;

//  test_dummy_points();
//  test_construction();
//  test_insert_1();
//  test_insert_point();
//  test_insert_hidden_point();
//  test_insert_hiding_point();
//  test_insert_a_point_twice();
//  test_insert_two_points_with_the_same_position();
//  test_remove();
//  test_27_to_1_sheeted_covering();
  test_insert_range(800, 7);
//    Iso_cuboid unitaire ->  0 <= weight < 0.015625
//  test_insert_rnd_as_delaunay(100, 0.);
//  test_insert_rnd_as_delaunay(100, 0.01);
//  test_insert_rnd_then_remove_all(5000, 7);
//  test_insert_rnd_then_remove_all(5000, 12);

  std::cout << "EXIT SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
