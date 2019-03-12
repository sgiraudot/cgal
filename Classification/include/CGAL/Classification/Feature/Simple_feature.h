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

#ifndef CGAL_CLASSIFICATION_SIMPLE_FEATURE_H
#define CGAL_CLASSIFICATION_SIMPLE_FEATURE_H

#include <CGAL/license/Classification.h>

#include <vector>
#include <CGAL/Classification/Feature_base.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on a user-defined scalar field.

    \tparam InputRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PropertyMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is statically castable to `double`.
  */
template <typename InputRange, typename PropertyMap>
class Simple_feature : public Feature_base
{
  const InputRange& m_input;
  PropertyMap m_pmap;
  
public:
  /*!
    \brief Constructs the feature using an input range and a property map.

    \param input point range.
    \param property_map property map to access scalar field.
    \param name name of the property (no default name is given).
  */
  Simple_feature (const InputRange& input,
                  PropertyMap property_map,
                  const std::string& name)
    : m_input (input), m_pmap (property_map)
  {
    this->set_name (name);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return static_cast<double> (get (m_pmap, *(m_input.begin()+pt_index)));
  }
  /// \endcond
};

class Mean_of_feature : public Feature_base
{
  std::vector<internal_float> m_values;
  
public:
  
  template <typename InputRange, typename ItemMap, typename PropertyMap, typename NeighborQuery>
  Mean_of_feature (const InputRange& input,
                   ItemMap item_map,
                   PropertyMap property_map,
                   const NeighborQuery& neighbor_query,
                   const std::string& name)
  {
    this->set_name ("mean_" + name);

    m_values.reserve(input.size());

    std::vector<std::size_t> neighborhood;
    for (typename InputRange::const_iterator it = input.begin();
         it != input.end(); ++ it)
    {
      neighbor_query (get (item_map, *it), std::back_inserter (neighborhood));

      double mean = 0.;
      for (std::size_t i = 0; i < neighborhood.size(); ++ i)
        mean += get (property_map, *(input.begin() + neighborhood[i]));
      m_values.push_back (mean / neighborhood.size());

      neighborhood.clear();
    }

  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return static_cast<double> (m_values[pt_index]);
  }
  /// \endcond
};

class Variance_of_feature : public Feature_base
{
  std::vector<internal_float> m_values;
  
public:
  
  template <typename InputRange, typename ItemMap, typename PropertyMap, typename NeighborQuery>
  Variance_of_feature (const InputRange& input,
                       const ItemMap& item_map,
                       PropertyMap property_map,
                       const NeighborQuery& neighbor_query,
                       Feature_handle mean_feature,
                       const std::string& name)
  {
    this->set_name ("variance_" + name);

    m_values.reserve(input.size());

    std::vector<std::size_t> neighborhood;
    std::size_t idx = 0;
    for (typename InputRange::const_iterator it = input.begin();
         it != input.end(); ++ it, ++ idx)
    {
      neighbor_query (get (item_map, *it), std::back_inserter (neighborhood));

      double mean = mean_feature->value(idx);
      double variance = 0.;
      for (std::size_t i = 0; i < neighborhood.size(); ++ i)
      {
        double v = get (property_map, *(input.begin() + neighborhood[i]));
        variance += (v - mean) * (v - mean);
      }
      m_values.push_back (variance / neighborhood.size());

      neighborhood.clear();
    }

  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return static_cast<double> (m_values[pt_index]);
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_SIMPLE_FEATURE_H
