// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_POINT_SET_NEIGHBORHOOD_H
#define CGAL_CLASSIFICATION_POINT_SET_NEIGHBORHOOD_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/boost/iterator/counting_iterator.hpp>

#include <CGAL/Kernel/Dimension_utils.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/centroid.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/Classification/Image.h>


namespace CGAL {

namespace Classification {

/// \cond SKIP_IN_MANUAL
namespace Access {

template <typename GeomTraits, typename PointMap, typename DimensionTag>
struct Search_structures;

template <typename GeomTraits, typename PointMap>
struct Search_structures<GeomTraits, PointMap, Dimension_tag<2> >
{
  typedef typename GeomTraits::Point_2 Point;
  typedef Search_traits_2<GeomTraits> SearchTraits;
  typedef Search_traits_adapter <boost::uint32_t, PointMap, SearchTraits> Search_traits;
  typedef Sliding_midpoint<Search_traits> Splitter;
  typedef Distance_adapter<boost::uint32_t, PointMap, Euclidean_distance<SearchTraits> > Distance;
  typedef Kd_tree<Search_traits, Splitter, Tag_true> Tree;
  typedef Fuzzy_sphere<Search_traits> Sphere;
  typedef Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Tree> Knn;

  template <typename Map>
  static void voxelize_point_set (std::size_t nb_pts, std::vector<boost::uint32_t>& indices, Map point_map,
                                  double voxel_size)
  {
    std::map<Point, std::vector<boost::uint32_t> > grid;

    for (boost::uint32_t i = 0; i < nb_pts; ++ i)
    {
      const Point& p = get(point_map, i);
      Point ref (std::floor(p.x() / voxel_size),
                 std::floor(p.y() / voxel_size));
      typename std::map<Point, std::vector<boost::uint32_t> >::iterator it
        = grid.insert (std::make_pair (ref, std::vector<boost::uint32_t>())).first;
      it->second.push_back (i);
    }

    for (typename std::map<Point, std::vector<boost::uint32_t> >::iterator
           it = grid.begin(); it != grid.end(); ++ it)
    {
      const std::vector<boost::uint32_t>& pts = it->second;
      Point centroid = CGAL::centroid (boost::make_transform_iterator
                                       (pts.begin(),
                                        CGAL::Property_map_to_unary_function<Map>(point_map)),
                                       boost::make_transform_iterator
                                       (pts.end(),
                                        CGAL::Property_map_to_unary_function<Map>(point_map)));
      boost::uint32_t chosen = 0;
      double min_dist = (std::numeric_limits<double>::max)();
      for (std::size_t i = 0; i < pts.size(); ++ i)
      {
        double dist = double(CGAL::squared_distance(get(point_map, pts[i]), centroid));
        if (dist < min_dist)
        {
          min_dist = dist;
          chosen = pts[i];
        }
      }
      indices.push_back (chosen);
    }
  }
};

template <typename GeomTraits, typename PointMap>
struct Search_structures<GeomTraits, PointMap, Dimension_tag<3> >
{
  typedef typename GeomTraits::Point_3 Point;
  typedef Search_traits_3<GeomTraits> SearchTraits;
  typedef Search_traits_adapter <boost::uint32_t, PointMap, SearchTraits> Search_traits;
  typedef Sliding_midpoint<Search_traits> Splitter;
  typedef Distance_adapter<boost::uint32_t, PointMap, Euclidean_distance<SearchTraits> > Distance;
  typedef Kd_tree<Search_traits, Splitter, Tag_true> Tree;
  typedef Fuzzy_sphere<Search_traits> Sphere;
  typedef Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Tree> Knn;
  
  template <typename Map>
  static void voxelize_point_set (std::size_t nb_pts, std::vector<boost::uint32_t>& indices, Map point_map,
                           double voxel_size)
  {
    std::map<Point, std::vector<boost::uint32_t> > grid;

    for (boost::uint32_t i = 0; i < nb_pts; ++ i)
    {
      const Point& p = get(point_map, i);
      Point ref (std::floor(p.x() / voxel_size),
                 std::floor(p.y() / voxel_size),
                 std::floor(p.z() / voxel_size));
      typename std::map<Point, std::vector<boost::uint32_t> >::iterator it
        = grid.insert (std::make_pair (ref, std::vector<boost::uint32_t>())).first;
      it->second.push_back (i);
    }

    for (typename std::map<Point, std::vector<boost::uint32_t> >::iterator
           it = grid.begin(); it != grid.end(); ++ it)
    {
      const std::vector<boost::uint32_t>& pts = it->second;
      Point centroid = CGAL::centroid (boost::make_transform_iterator
                                       (pts.begin(),
                                        CGAL::Property_map_to_unary_function<Map>(point_map)),
                                       boost::make_transform_iterator
                                       (pts.end(),
                                        CGAL::Property_map_to_unary_function<Map>(point_map)));
      boost::uint32_t chosen = 0;
      double min_dist = (std::numeric_limits<double>::max)();
      for (std::size_t i = 0; i < pts.size(); ++ i)
      {
        double dist = double(CGAL::squared_distance(get(point_map, pts[i]), centroid));
        if (dist < min_dist)
        {
          min_dist = dist;
          chosen = pts[i];
        }
      }
      indices.push_back (chosen);
    }
  }
};

} // namespace Access
/// \endcond

  /*!
    \ingroup PkgClassificationPointSet

    \brief Class that precomputes spatial searching structures for an
    input point set and gives access to the local neighborhood of a
    point as a set of indices.

    It allows the user to generate models of `NeighborQuery` based on
    a fixed range neighborhood or on a fixed K number of neighbors. In
    addition, the spatial searching structures can be computed on a
    simplified version of the point set to allow for neighbor queries
    at a higher scale.

    \tparam GeomTraits is a model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.
  */
template <typename GeomTraits, typename PointRange, typename PointMap, typename DimensionTag = Dimension_tag<3> >
class Point_set_neighborhood
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename CGAL::Access::Point<GeomTraits, DimensionTag>::type Point;
  
  class My_point_property_map
  {
    const PointRange* input;
    PointMap point_map;
    
  public:
    typedef Point value_type;
    typedef const value_type& reference;
    typedef boost::uint32_t key_type;
    typedef boost::lvalue_property_map_tag category;
    My_point_property_map () { }
    My_point_property_map (const PointRange *input, PointMap point_map)
      : input (input), point_map (point_map) { }
    reference operator[] (key_type k) const { return reinterpret_cast<reference>(get(point_map, *(input->begin()+std::size_t(k)))); }
    friend inline reference get (const My_point_property_map& ppmap, key_type i)
    { return ppmap[i]; }
  };

  typedef Access::Search_structures<GeomTraits, My_point_property_map, DimensionTag> Search_structures;
  typedef typename Search_structures::Search_traits Search_traits;
  typedef typename Search_structures::Tree Tree;
  typedef typename Search_structures::Distance Distance;
  typedef typename Search_structures::Splitter Splitter;
  typedef typename Search_structures::Sphere Sphere;
  typedef typename Search_structures::Knn Knn;

  Tree* m_tree;
  Distance m_distance;
  
public:

  /*!
    Functor that computes the neighborhood of an input point with a
    fixed number of neighbors.

    \cgalModels CGAL::Classification::NeighborQuery

    \sa Point_set_neighborhood
  */
  class K_neighbor_query
  {
  public:
    typedef typename Point_set_neighborhood::Point_3 value_type; ///<
  private:
    const Point_set_neighborhood& neighborhood;
    unsigned int k;
  public:
    /*!
      \brief Constructs a K neighbor query object.
      \param neighborhood point set neighborhood object.
      \param k number of neighbors per query.
    */
    K_neighbor_query (const Point_set_neighborhood& neighborhood, unsigned int k)
      : neighborhood (neighborhood), k(k) { }

    /// \cond SKIP_IN_MANUAL
    template <typename OutputIterator>
    OutputIterator operator() (const value_type& query, OutputIterator output) const
    {
      neighborhood.k_neighbors (query, k, output);
      return output;
    }
    /// \endcond
  };

  /*!
    Functor that computes the neighborhood of an input point defined
    as the points lying in a sphere of fixed radius centered at the
    input point.

    \cgalModels CGAL::Classification::NeighborQuery

    \sa Point_set_neighborhood
  */
  class Sphere_neighbor_query
  {
  public:
    typedef typename Point_set_neighborhood::Point_3 value_type; ///<
  private:
    const Point_set_neighborhood& neighborhood;
    double radius;
  public:
    /*!
      \brief Constructs a range neighbor query object.
      \param neighborhood point set neighborhood object.
      \param radius radius of the neighbor query sphere.
    */
    Sphere_neighbor_query (const Point_set_neighborhood& neighborhood, double radius)
      : neighborhood (neighborhood), radius(radius) { }

    /// \cond SKIP_IN_MANUAL
    template <typename OutputIterator>
    OutputIterator operator() (const value_type& query, OutputIterator output) const
    {
      neighborhood.sphere_neighbors (query, radius, output);
      return output;
    }
    /// \endcond
  };

  /// \cond SKIP_IN_MANUAL
  friend class K_neighbor_query;
  friend class Sphere_neighbor_query;

  Point_set_neighborhood () : m_tree (NULL) { }
  /// \endcond

  /// \name Constructors
  /// @{

  /*!
    \brief Constructs a neighborhood object based on the input range.

    \param input point range.
    \param point_map property map to access the input points.
  */
  Point_set_neighborhood (const PointRange& input,
                          PointMap point_map)
    : m_tree (NULL)
  {
    My_point_property_map pmap (&input, point_map);
    m_tree = new Tree (boost::counting_iterator<boost::uint32_t> (0),
                       boost::counting_iterator<boost::uint32_t> (boost::uint32_t(input.size())),
                       Splitter(),
                       Search_traits (pmap));
    m_distance = Distance (pmap);
    m_tree->build();
  }

  /*!
    \brief Constructs a simplified neighborhood object based on the input range.

    This method first computes a simplified version of the input point
    set by voxelization: a 3D grid is defined and for each subset
    present in one cell, only the point closest to the centroid of
    this subset is used.

    \param input input range.
    \param point_map property map to access the input points.
    \param voxel_size size of the cells of the 3D grid used for simplification.
  */
  Point_set_neighborhood (const PointRange& input,
                          PointMap point_map,
                          double voxel_size)
    : m_tree (NULL)
  {
    // First, simplify
    std::vector<boost::uint32_t> indices;
    My_point_property_map pmap (&input, point_map);
    Search_structures::voxelize_point_set(input.size(), indices, pmap, voxel_size);
    
    m_tree = new Tree (indices.begin(), indices.end(),
                       Splitter(),
                       Search_traits (pmap));
    
    m_distance = Distance (pmap);
    m_tree->build();
  }

  /// @}
  
  /// \cond SKIP_IN_MANUAL
  ~Point_set_neighborhood ()
  {
    if (m_tree != NULL)
      delete m_tree;
  }
  /// \endcond

  /// \name Queries
  /// @{

  /*!
    \brief Returns a neighbor query object with fixed number of neighbors `k`.
  */
  K_neighbor_query k_neighbor_query (const unsigned int k) const
  {
    return K_neighbor_query (*this, k);
  }

  /*!
    \brief Returns a neighbor query object with fixed radius `radius`.
  */
  Sphere_neighbor_query sphere_neighbor_query (const double radius) const
  {
    return Sphere_neighbor_query (*this, radius);
  }

  /// @}

private:

  template <typename OutputIterator>
  void sphere_neighbors (const Point_3& query, const FT radius_neighbors, OutputIterator output) const
  {
    CGAL_assertion (m_tree != NULL);
    Sphere fs (reinterpret_cast<const Point&>(query), radius_neighbors, 0, m_tree->traits());
    m_tree->search (output, fs);
  }

  template <typename OutputIterator>
  void k_neighbors (const Point_3& query, const unsigned int k, OutputIterator output) const
  {
    CGAL_assertion (m_tree != NULL);
    Knn search (*m_tree, reinterpret_cast<const Point&>(query), k, 0, true, m_distance);
    for (typename Knn::iterator it = search.begin(); it != search.end(); ++ it)
      *(output ++) = it->first;
  }

};
  

}
  
}


#endif // CGAL_CLASSIFICATION_POINT_SET_POINT_SET_NEIGHBORHOOD_H
