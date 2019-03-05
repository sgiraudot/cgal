#ifndef CGAL_CLASSIFICATION_JOINT_CLASSIFICATION_WRAPPER_H
#define CGAL_CLASSIFICATION_JOINT_CLASSIFICATION_WRAPPER_H

namespace CGAL {

namespace Classification {

class Joint_classification_wrapper
{
public:
  
  typedef boost::uint32_t Index_type;
  typedef std::pair<Index_type, Index_type> Pair;
  typedef int unspecified_type;
    
  class Joint_feature : public Feature_base
  {
    Feature_handle m_feature;
    const std::vector<Pair>& m_neighborhood;
    bool m_self;
    
  public:
    Joint_feature (Feature_handle feature, const std::vector<Pair>& neighborhood, bool self)
      : m_feature (feature), m_neighborhood (neighborhood), m_self (self)
    {
      if (m_self)
        this->set_name (feature->name() + std::string("_self"));
      else
        this->set_name (feature->name() + std::string("_other"));
    }

    double value (std::size_t idx)
    {
      if (m_self)
        return m_feature->value(m_neighborhood[idx].first);
      return m_feature->value(m_neighborhood[idx].second);
    }
  };

  template <typename Classifier>
  class Joint_classifier
  {
    const std::vector<Pair>& m_neighborhood;
    const std::vector<Index_type>& m_neighborhood_index;
    std::size_t m_nb_labels;
    const Classifier& m_classifier;
    bool m_complete;

  public:

    Joint_classifier (const std::vector<Pair>& neighborhood,
                      const std::vector<Index_type>& neighborhood_index,
                      std::size_t nb_labels,
                      const Classifier& classifier,
                      bool complete)
      : m_neighborhood (neighborhood)
      , m_neighborhood_index (neighborhood_index)
      , m_nb_labels (nb_labels)
      , m_classifier (classifier)
      , m_complete (complete)
    { }

    void operator() (std::size_t item_index, std::vector<double>& out) const
    {
      std::size_t first = m_neighborhood_index[item_index];
      std::size_t last = m_neighborhood_index[item_index+1];

      out.resize(m_nb_labels, 0.f);
      
      for (std::size_t i = first; i < last; ++ i)
      {
        std::vector<double> local_out;
        m_classifier(i, local_out);
        if (m_complete)
          for (std::size_t j = 0; j < out.size(); ++ j)
            for (std::size_t k = 0; k < m_nb_labels; ++ k)
              out[j] += local_out[j * m_nb_labels + k];
        else
          for (std::size_t j = 0; j < out.size(); ++ j)
            out[j] += local_out[j] + local_out[j + m_nb_labels];
      }
      for (std::size_t i = 0; i < out.size(); ++ i)
        out[i] /= (last - first);
    }

    bool complete() const { return m_complete; }

    std::size_t first (std::size_t item_index) const
    {
      return m_neighborhood_index[item_index];
    }
    std::size_t last (std::size_t item_index) const
    {
      return m_neighborhood_index[item_index + 1];
    }

    const Classifier& base() const { return m_classifier; }

    std::size_t neighbor (std::size_t index) const
    {
      return m_neighborhood[index].second;
    }
  };

private:

  std::vector<Pair> m_neighborhood;
  std::vector<Index_type> m_neighborhood_index;
  
  Label_set m_labels;
  Feature_set m_features;
  std::vector<int> m_ground_truth;

  bool m_complete;
  
public:

  template <typename ItemRange, typename ItemMap, typename NeighborQuery>
  Joint_classification_wrapper (const ItemRange& input,
                                const ItemMap item_map,
                                const NeighborQuery& neighbor_query,
                                bool complete = false,
                                std::size_t allocation_hint = 0)
    : m_complete (complete)
  {
    std::vector<std::size_t> neighbors;
    Index_type idx = 0;
    m_neighborhood_index.reserve (input.size() + 1);
    m_neighborhood.reserve (allocation_hint);
    
    for (typename ItemRange::const_iterator it = input.begin(); it != input.end(); ++ it, ++ idx)
    {
      neighbors.clear();
      neighbor_query (get (item_map, *it), std::back_inserter (neighbors));
      
      m_neighborhood_index.push_back(m_neighborhood.size());
      for (std::size_t i = 0; i < neighbors.size(); ++ i)
        if (idx != Index_type(neighbors[i]))
          m_neighborhood.push_back (std::make_pair (idx, Index_type(neighbors[i])));
    }
    m_neighborhood_index.push_back(m_neighborhood.size());
  }

  const Label_set& labels (const Label_set& l)
  {
    m_labels.clear();

    if (m_complete)
      for (std::size_t i = 0; i < l.size(); ++ i)
        for (std::size_t j = 0; j < l.size(); ++ j)
        {
          std::string lname = l[i]->name() + "_next_to_" + l[j]->name();
          m_labels.add (lname.c_str());
        }
    else
    {
      for (std::size_t i = 0; i < l.size(); ++ i)
      {
        std::string lname = l[i]->name() + "_same";
        m_labels.add (lname.c_str());
      }
      for (std::size_t i = 0; i < l.size(); ++ i)
      {
        std::string lname = l[i]->name() + "_different";
        m_labels.add (lname.c_str());
      }
    }

    return m_labels;
  }

  const Feature_set& features (const Feature_set& f)
  {
    m_features.clear();
    for (std::size_t i = 0; i < f.size(); ++ i)
      m_features.add<Joint_feature>(f[i], m_neighborhood, true);
    for (std::size_t i = 0; i < f.size(); ++ i)
      m_features.add<Joint_feature>(f[i], m_neighborhood, false);
    return m_features;
  }

  template <typename LabelIndexRange>
  const std::vector<int>& ground_truth (const LabelIndexRange& gt)
  {
    std::size_t label_size = m_labels.size() / 2;
    if (m_complete)
      label_size = std::sqrt (m_labels.size());
    
    m_ground_truth.clear();
    m_ground_truth.reserve(m_neighborhood.size());
    for (std::size_t i = 0; i < gt.size(); ++ i)
    {
      std::size_t first = m_neighborhood_index[i];
      std::size_t last = m_neighborhood_index[i+1];
      if (m_complete)
      {
        for (std::size_t j = first; j < last; ++ j)
          if (int(gt[i]) == -1 || gt[m_neighborhood[j].second] == -1)
            m_ground_truth.push_back (-1);
          else 
            m_ground_truth.push_back(gt[i] * label_size + gt[m_neighborhood[j].second]);
      }
      else
      {
        for (std::size_t j = first; j < last; ++ j)
          if (int(gt[i]) == -1)
            m_ground_truth.push_back (-1);
          else if (gt[i] == gt[m_neighborhood[j].second])
            m_ground_truth.push_back(gt[i]);
          else
            m_ground_truth.push_back(gt[i] + label_size);
      }
    }
    return m_ground_truth;
  }

  template <typename Classifier>
  Joint_classifier<Classifier> classifier (const Classifier& classifier)
  {
    if (m_complete)
      return Joint_classifier<Classifier> (m_neighborhood, m_neighborhood_index, std::sqrt(m_labels.size()), classifier, m_complete);

    return Joint_classifier<Classifier> (m_neighborhood, m_neighborhood_index, m_labels.size() / 2, classifier, m_complete);
  }

  unspecified_type graphcut_neighbor_query() const
  {
    return 0;
  }

};

/// \cond SKIP_IN_MANUAL
namespace internal {

  template <typename ItemRange, typename ItemMap,
            typename Classifier,
            typename LabelIndexRange>
  class Classify_functor_joint
  {
    typedef Joint_classification_wrapper::Joint_classifier<Classifier> Joint_classifier;
    const ItemRange& m_input;
    ItemMap m_item_map;
    const Label_set& m_labels;
    const Joint_classifier& m_classifier;
    double m_strength;
    const std::vector<std::vector<std::size_t> >& m_indices;
    const std::vector<std::pair<std::size_t, std::size_t> >& m_input_to_indices;
    LabelIndexRange& m_out;

#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
    typedef CGAL::internal::Alpha_expansion_graph_cut_boost             Alpha_expansion;
#else
    typedef CGAL::internal::Alpha_expansion_graph_cut_boykov_kolmogorov Alpha_expansion;
#endif
    
  public:

    Classify_functor_joint (const ItemRange& input,
                            ItemMap item_map,
                            const Label_set& labels,
                            const Joint_classifier& classifier,
                            double strength,
                            const std::vector<std::vector<std::size_t> >& indices,
                            const std::vector<std::pair<std::size_t, std::size_t> >& input_to_indices,
                            LabelIndexRange& out)
    : m_input (input), m_item_map (item_map), m_labels (labels),
      m_classifier (classifier),
      m_strength (strength), m_indices (indices), m_input_to_indices (input_to_indices), m_out (out)
    { }

#ifdef CGAL_LINKED_WITH_TBB
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        apply(s);
    }
#endif // CGAL_LINKED_WITH_TBB
    
    inline void apply (std::size_t sub) const
    {
      if (m_indices[sub].empty())
        return;
        
      std::vector<std::pair<std::size_t, std::size_t> > edges;
      std::vector<double> edge_weights;
      std::vector<std::vector<double> > probability_matrix
        (m_labels.size(), std::vector<double>(m_indices[sub].size(), 0.));
      std::vector<std::size_t> assigned_label (m_indices[sub].size());

      for (std::size_t j = 0; j < m_indices[sub].size(); ++ j)
      {
        std::size_t s = m_indices[sub][j];
        
        std::vector<double> values(m_labels.size(), 0.f);

        std::size_t first = m_classifier.first(s);
        std::size_t last = m_classifier.last(s);
        
        if (first == last)
          continue;
        
        for (std::size_t i = first; i < last; ++ i)
        {
          std::size_t neighbor_s = m_classifier.neighbor(i);
          
          std::vector<double> values_local;
          
          m_classifier.base()(i, values_local);
          
          if (m_classifier.complete())
            for (std::size_t k = 0; k < values.size(); ++ k)
              for (std::size_t l = 0; l < m_labels.size(); ++ l)
                values[k] += values_local[k * m_labels.size() + l];
          else
            for (std::size_t k = 0; k < values.size(); ++ k)
              values[k] += values_local[k] + values_local[k + values.size()];
          
          if (sub == m_input_to_indices[neighbor_s].first
              && j != m_input_to_indices[neighbor_s].second)
          {
            edges.push_back (std::make_pair (j, m_input_to_indices[neighbor_s].second));

            double weight = 0.;
            if (m_classifier.complete())
              for (std::size_t k = 0; k < m_labels.size(); ++ k)
                weight += values_local[k * m_labels.size() + k];
            else
              for (std::size_t k = 0; k < m_labels.size(); ++ k)
                weight += values_local[k];

            edge_weights.push_back (weight * m_strength);
          }
        }
        
        std::size_t nb_class_best = 0;
        double val_class_best = 0.f;
        for(std::size_t k = 0; k < m_labels.size(); ++ k)
        {
          double value = values[k] / (last - first);
          probability_matrix[k][j] = -std::log(value);
            
          if(val_class_best < value)
          {
            val_class_best = value;
            nb_class_best = k;
          }
        }

        assigned_label[j] = nb_class_best;
      }
    
      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, probability_matrix, assigned_label);

      for (std::size_t i = 0; i < assigned_label.size(); ++ i)
        m_out[m_indices[sub][i]] = static_cast<typename LabelIndexRange::iterator::value_type>(assigned_label[i]);
    }

  };

} // namespace internal
  
/// \endcond
  

  /*! 
    \ingroup PkgClassificationMain

    \brief Runs the classification algorithm with a global
    regularization based on a graph cut.

    The computed classification energy is globally regularized through
    an alpha-expansion algorithm. This method is slow but provides
    the user with good quality results.

    To speed up computation, the input domain can be subdivided into
    smaller subsets such that several smaller graph cuts are applied
    instead of a big one. The computation of these smaller graph cuts can
    be done in parallel. Increasing the number of subsets allows for
    faster computation times but can also reduce the quality of the
    results.

    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Parallel_tag` or `Sequential_tag`.
    \tparam ItemRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`.
    \tparam ItemMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `ItemRange` and value type
    is the type of item to classify (for example, `CGAL::Point_3`).
    \tparam NeighborQuery model of `NeighborQuery`.
    \tparam Classifier model of `Classifier`.
    \tparam Model of `Range` with random access iterators whose value
    type is an integer type.

    \param input input range.
    \param item_map property map to access the input items.
    \param labels set of input labels.
    \param classifier input classifier.
    \param neighbor_query used to access neighborhoods of items.
    \param strength strength of the regularization with respect to the
    classification energy. Higher values produce more regularized
    output but may result in a loss of details.
    \param min_number_of_subdivisions minimum number of subdivisions
    (for parallel processing to be efficient, this should be at least
    the number of cores of the processor).
    \param output where to store the result. It is stored as a sequence,
    ordered like the input range, containing for each point the index
    (in the `Label_set`) of the assigned label.

  */
  template <typename ConcurrencyTag,
            typename ItemRange,
            typename ItemMap,
            typename Classifier,
            typename LabelIndexRange>
  void classify_with_graphcut (const ItemRange& input,
                               const ItemMap item_map,
                               const Label_set& labels,
                               const Joint_classification_wrapper::Joint_classifier<Classifier>& classifier,
                               const Joint_classification_wrapper::unspecified_type&,
                               const double strength,
                               const std::size_t min_number_of_subdivisions,
                               LabelIndexRange& output)
  {
    CGAL::Bbox_3 bbox = CGAL::bbox_3
      (boost::make_transform_iterator (input.begin(), CGAL::Property_map_to_unary_function<ItemMap>(item_map)),
       boost::make_transform_iterator (input.end(), CGAL::Property_map_to_unary_function<ItemMap>(item_map)));

    double Dx = double(bbox.xmax() - bbox.xmin());
    double Dy = double(bbox.ymax() - bbox.ymin());
    double A = Dx * Dy;
    double a = A / min_number_of_subdivisions;
    double l = std::sqrt(a);
    std::size_t nb_x = std::size_t(Dx / l) + 1;
    std::size_t nb_y = std::size_t((A / nb_x) / a) + 1;
    std::size_t nb = nb_x * nb_y;

    std::vector<CGAL::Bbox_3> bboxes;
    bboxes.reserve(nb);
    for (std::size_t x = 0; x < nb_x; ++ x)
      for (std::size_t y = 0; y < nb_y; ++ y)
      {
        bboxes.push_back
          (CGAL::Bbox_3 (bbox.xmin() + Dx * (x / double(nb_x)),
                         bbox.ymin() + Dy * (y / double(nb_y)),
                         bbox.zmin(),
                         bbox.xmin() + Dx * ((x+1) / double(nb_x)),
                         bbox.ymin() + Dy * ((y+1) / double(nb_y)),
                         bbox.zmax()));
      }

#ifdef CGAL_CLASSIFICATION_VERBOSE
    std::cerr << "Number of divisions = " << nb_x * nb_y << std::endl;
    std::cerr << " -> Size of division: " << Dx / nb_x << " " << Dy / nb_y << std::endl;
#endif

    std::vector<std::vector<std::size_t> > indices (nb);
    std::vector<std::pair<std::size_t, std::size_t> > input_to_indices(input.size());
    
    for (std::size_t s = 0; s < input.size(); ++ s)
    {
      CGAL::Bbox_3 b = get(item_map, *(input.begin() + s)).bbox();
        
      for (std::size_t i = 0; i < bboxes.size(); ++ i)
        if (CGAL::do_overlap (b, bboxes[i]))
        {
          input_to_indices[s] = std::make_pair (i, indices[i].size());
          indices[i].push_back (s);
          break;
        }
    }

    internal::Classify_functor_joint<ItemRange, ItemMap, Classifier, LabelIndexRange>
      f (input, item_map, labels, classifier, strength, indices, input_to_indices, output);
    
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, indices.size ()), f);
    }
    else
#endif
    {
      for (std::size_t sub = 0; sub < indices.size(); ++ sub)
        f.apply (sub);
    }
  }

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_JOINT_CLASSIFICATION_WRAPPER_H
