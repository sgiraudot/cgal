// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
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

#ifndef CLASSIFICATION_K_MEANS_CLASSIFIER_H
#define CLASSIFICATION_K_MEANS_CLASSIFIER_H

#include <CGAL/license/Classification.h>

#include <CGAL/Random.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/Classification/internal/verbosity.h>
#include <CGAL/tags.h>
#include <CGAL/algorithm.h>

#include <iostream>

namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationClassifiers

  \brief %Classifier 

  \cgalModels `CGAL::Classification::Classifier`
*/
class K_means_classifier
{
  struct Cluster
  {
    std::vector<double> mean;
  };
  
  const Label_set& m_labels;
  const Feature_set& m_features;
  std::vector<double> m_feature_means;
  std::vector<double> m_feature_sd;
  std::vector<Cluster> m_clusters;
  
public:

  /// \name Constructor
  /// @{
  
/*!

  \brief Instantiate the classifier using the sets of `labels` and `features`.

  \note If the label set of the feature set are modified after
  instantiating this object (addition of removal of a label and/or of
  a feature), another classifier object should be instantiated as the
  internal data structures of this one are invalidated.
*/
  K_means_classifier (const Label_set& labels,
                      const Feature_set& features)
    : m_labels (labels), m_features (features)
  {
  }

  /// @}

  /// \cond SKIP_IN_MANUAL
  void operator() (std::size_t item_index,
                   std::vector<float>& out) const
  {
    out.resize (m_labels.size());
    float total = 0.f;
    for (std::size_t l = 0; l < m_labels.size(); ++ l)
    {
      out[l] = std::exp(-distance_to_cluster_mean (item_index, l));
      total += out[l];
    }
    
    for (std::size_t l = 0; l < m_labels.size(); ++ l)
      out[l] /= total;
  }
  
  void compute_normalization_coefficients (std::size_t number_of_points)
  {
    m_feature_means.clear();
    m_feature_means.resize (m_features.size(), 0.f);
    
    for (std::size_t i = 0; i < number_of_points; ++ i)
      for (std::size_t f = 0; f < m_features.size(); ++ f)
        m_feature_means[f] += m_features[f]->value(i);
    for (std::size_t f = 0; f < m_features.size(); ++ f)
      m_feature_means[f] /= number_of_points;

    m_feature_sd.clear();
    m_feature_sd.resize (m_features.size(), 0.f);

    for (std::size_t i = 0; i < number_of_points; ++ i)
      for (std::size_t f = 0; f < m_features.size(); ++ f)
        m_feature_sd[f] +=
          (m_features[f]->value(i) - m_feature_means[f]) *
          (m_features[f]->value(i) - m_feature_means[f]);
    for (std::size_t f = 0; f < m_features.size(); ++ f)
    {
      m_feature_sd[f] = std::sqrt (m_feature_sd[f] / number_of_points);
      if (m_feature_sd[f] == 0.f)
        m_feature_sd[f] = 1.f;

      if (std::isnan(m_feature_means[f]))
        m_feature_means[f] = 0.f;
      if (std::isnan(m_feature_sd[f]))
        m_feature_sd[f] = 1.f;
    }
  }

  float normalized (std::size_t feature_index, std::size_t item_index) const
  {
    return float((m_features[feature_index]->value(item_index) - m_feature_means[feature_index])
                 / m_feature_sd[feature_index]);
  }
  /// \endcond

  /// \name Training
  /// @{

  /*!
    \brief Runs the training algorithm.
  */
  void train (std::size_t number_of_points)
  {
    CGAL::Random& random = CGAL::get_default_random();

    compute_normalization_coefficients (number_of_points);

    std::vector<std::vector<float> > features(m_features.size());
    for (std::size_t i = 0; i < features.size(); ++ i)
    {
      features[i].reserve (number_of_points);
      for (std::size_t j = 0; j < number_of_points; ++ j)
        features[i].push_back(normalized(i,j));
    }

    // Initialize
    m_clusters.clear();
    m_clusters.resize(m_labels.size());

    std::set<std::size_t> done;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      std::size_t idx;
      do
      {
        idx = std::size_t(random.get_int (0, number_of_points));
      } while (!done.insert(idx).second);
      
      for (std::size_t j = 0; j < m_features.size(); ++ j)
        m_clusters[i].mean.push_back (features[j][idx]);
    }


    std::cerr << "Clusters initialized with:" << std::endl;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      std::cerr << " * " << m_labels[i]->name() << ":";
      for (std::size_t j = 0; j < m_features.size(); ++ j)
        std::cerr << " " << m_clusters[i].mean[j];
      std::cerr << std::endl;
    }

    std::vector<int> assigned_labels (number_of_points, -1);

    std::size_t iter = 0;
    std::size_t stuck_time = 0;
    std::size_t nb_changed_before = 0;
    
    while (true)
    {
      std::size_t nb_changed_points = 0;
      ++ iter;
      
      // Redistribute points
      for (std::size_t i = 0; i < number_of_points; ++ i)
      {
        int label_before = assigned_labels[i];

        float d_min = std::numeric_limits<float>::max();
        int label_after = -1;
        
        for (std::size_t j = 0; j < m_labels.size(); ++ j)
        {
          float d = distance_to_cluster_mean (i, j, features);
          if (d < d_min)
          {
            d_min = d;
            label_after = j;
          }
        }

        if (label_after != label_before)
          nb_changed_points ++;
        assigned_labels[i] = label_after;
      }

      std::cerr << " * Iteration " << iter << ": " << nb_changed_points << " point(s) changed, ";
      
      // Stop at convergence
      if (nb_changed_points == 0)
        break;

      if (nb_changed_points == nb_changed_before)
        ++ stuck_time;
      nb_changed_before = nb_changed_points;

      if (stuck_time == 10)
        break;

      std::vector<std::size_t> nb_inliers(m_labels.size(), 0);
      // Else, recompute means
      for (std::size_t i = 0; i < m_labels.size(); ++ i)
        for (std::size_t j = 0; j < m_features.size(); ++ j)
          m_clusters[i].mean[j] = 0.f;
      
      for (std::size_t i = 0; i < number_of_points; ++ i)
      {
        std::size_t label = std::size_t(assigned_labels[i]);
        nb_inliers[label] ++;

        for (std::size_t j = 0; j < m_features.size(); ++ j)
          m_clusters[label].mean[j] += features[j][i];
      }
      
      for (std::size_t i = 0; i < m_labels.size(); ++ i)
        for (std::size_t j = 0; j < m_features.size(); ++ j)
          m_clusters[i].mean[j] /= nb_inliers[i];

      std::cerr << "distribution = [ ";
      for (std::size_t i = 0; i < nb_inliers.size(); ++ i)
        std::cerr << nb_inliers[i] << " ";
      std::cerr << "]" << std::endl;
    }
    
    std::cerr << "Clusters finalized as:" << std::endl;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      std::cerr << " * " << m_labels[i]->name() << ":";
      for (std::size_t j = 0; j < m_features.size(); ++ j)
        std::cerr << " " << m_clusters[i].mean[j];
      std::cerr << std::endl;
    }

  }

  /// @}

  /// \endcond

  /// \name Input/Output
  /// @{
  
  /*!
    \brief Saves the current configuration in the stream `output`.
  */
  void save_configuration (std::ostream&)
  {
  }
  
  /*!
    \brief Loads a configuration from the stream `input`.  */
  bool load_configuration (std::istream&)
  {
    return false;
  }

  /// @}

private:

  float distance_to_cluster_mean (std::size_t item_index, std::size_t cluster_index) const
  {
    float out = 0.f;

    for (std::size_t i = 0; i < m_features.size(); ++ i)
    {
      float v = normalized(i,item_index);
      out += (v - m_clusters[cluster_index].mean[i]) * (v - m_clusters[cluster_index].mean[i]);
    }

    out /= m_features.size();
    return out;
  }

  float distance_to_cluster_mean (std::size_t item_index, std::size_t cluster_index,
                                  const std::vector<std::vector<float> >& features) const
  {
    float out = 0.f;

    for (std::size_t i = 0; i < m_features.size(); ++ i)
    {
      float v = features[i][item_index];
      out += (v - m_clusters[cluster_index].mean[i]) * (v - m_clusters[cluster_index].mean[i]);
    }

    out /= m_features.size();
    return out;
  }

};

}

}

#endif //  CLASSIFICATION_K_MEANS_CLASSIFIER_H
