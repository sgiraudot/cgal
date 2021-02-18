#ifndef CGAL_ARR_LIGHTWEIGHT_POLYLINE_TRAITS_2_H
#define CGAL_ARR_LIGHTWEIGHT_POLYLINE_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <fstream>

#include <boost/variant.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_lightweight_polyline_subtraits_2.h>

namespace CGAL {

template <typename Kernel_, typename PointIterator>
class Arr_lightweight_polyline_traits_2
  : public Arr_polycurve_basic_traits_2
    <Arr_lightweight_polyline_subtraits_2<Kernel_, PointIterator>,
     internal::Lightweight_polyline_2<Kernel_, PointIterator> >
{
public:

  using Kernel = Kernel_;
  using Point_iterator = PointIterator;
  using Subcurve_traits_2 = Arr_lightweight_polyline_subtraits_2<Kernel, Point_iterator>;
  using Curve_2 = internal::Lightweight_polyline_2<Kernel, Point_iterator>;

private:
  using Base = Arr_polycurve_basic_traits_2<Subcurve_traits_2, Curve_2>;
  using Self = Arr_lightweight_polyline_traits_2<Kernel, Point_iterator>;

public:

  using Has_left_category = typename Base::Has_left_category;
  using Has_do_intersect_category = typename Base::Has_do_intersect_category;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  using Are_all_sides_oblivious_tag = typename Base::Are_all_sides_oblivious_tag;

  using X_monotone_subcurve_2 = typename Base::X_monotone_subcurve_2;
  using Size = typename Base::Size;
  using size_type = typename Base::size_type;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  using Compare_x_2 = typename Base::Compare_x_2;
  using Compare_xy_2 = typename Base::Compare_xy_2;
  using Construct_min_vertex_2 = typename Base::Construct_min_vertex_2;
  using Construct_max_vertex_2 = typename Base::Construct_max_vertex_2;
  using Is_vertical_2 = typename Base::Is_vertical_2;
  using Compare_y_at_x_2 = typename Base::Compare_y_at_x_2;
  using Compare_y_at_x_left_2 = typename Base::Compare_y_at_x_left_2;
  using Compare_y_at_x_right_2 = typename Base::Compare_y_at_x_right_2;
  using Equal_2 = typename Base::Equal_2;
  using Compare_endpoints_xy_2 = typename Base::Compare_endpoints_xy_2;
  using Construct_opposite_2 = typename Base::Construct_opposite_2;
  using Approximate_2 = typename Base::Approximate_2;
  using Construct_x_monotone_curve_2 = typename Base::Construct_x_monotone_curve_2;
  using Parameter_space_in_x_2 = typename Base::Parameter_space_in_x_2;
  using Parameter_space_in_y_2 = typename Base::Parameter_space_in_y_2;
  using Compare_x_on_boundary_2 = typename Base::Compare_x_on_boundary_2;
  using Compare_x_at_limit_2 = typename Base::Compare_x_at_limit_2;
  using Compare_x_near_boundary_2 = typename Base::Compare_x_near_boundary_2;
  using Compare_x_near_limit_2 = typename Base::Compare_x_near_limit_2;
  using Compare_y_on_boundary_2 = typename Base::Compare_y_on_boundary_2;
  using Compare_y_near_boundary_2 = typename Base::Compare_y_near_boundary_2;
  using Is_on_y_identification_2 = typename Base::Is_on_y_identification_2;
  using Is_on_x_identification_2 = typename Base::Is_on_x_identification_2;
  using Trim_2 = typename Base::Trim_2;

  using Has_merge_category = typename Subcurve_traits_2::Has_merge_category;
  using Multiplicity = typename Subcurve_traits_2::Multiplicity;
  using Subcurve_2 = typename Subcurve_traits_2::Curve_2;


  Arr_lightweight_polyline_traits_2() : Base() { }
  Arr_lightweight_polyline_traits_2(const Subcurve_traits_2* geom_traits) : Base(geom_traits) { }

  class Make_x_monotone_2
  {
  protected:
    using Traits = Arr_lightweight_polyline_traits_2<Kernel, Point_iterator>;
    const Traits& m_traits;
    Make_x_monotone_2(const Traits& traits) : m_traits(traits) {}
    friend class Arr_lightweight_polyline_traits_2<Kernel, Point_iterator>;
  public:
    using Curve_iterator = typename Curve_2::iterator;

    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      using Make_x_monotone_result = boost::variant<Point_2, X_monotone_curve_2>;

      auto compare_x_2 = m_traits.compare_x_2_object();

      Curve_iterator start = cv.points_begin();
      Curve_iterator it_prev = start + 1;
      Curve_iterator it_curr = it_prev + 1;

      Comparison_result previous_comp = compare_x_2 (*start, *it_prev);
      for (; it_curr != cv.points_end(); ++ it_prev, ++ it_curr)
      {
        Comparison_result current_comp = compare_x_2 (*it_prev, *it_curr);
        if (current_comp != previous_comp)
        {
          previous_comp = current_comp;
          *oi ++ = Make_x_monotone_result(X_monotone_curve_2 (start, it_curr));
          start = it_prev;
        }
      }
      *oi ++ = Make_x_monotone_result(X_monotone_curve_2(start, cv.points_end()));

      return oi;
    }
  };
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this); }

  class Split_2
  {
  protected:
    using Traits = Arr_lightweight_polyline_traits_2<Kernel, Point_iterator>;
    const Traits& m_traits;
    Split_2(const Traits& traits) : m_traits(traits) {}
    friend class Arr_lightweight_polyline_traits_2<Kernel, Point_iterator>;
  public:
    void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const
    {
      const Subcurve_traits_2* geom_traits = m_traits.subcurve_traits_2();
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto equal = geom_traits->equal_2_object();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition((! equal(m_traits.
                                 construct_min_vertex_2_object()(xcv), p)));
      CGAL_precondition((! equal(m_traits.
                                 construct_max_vertex_2_object()(xcv), p)));

      CGAL_precondition_msg(xcv.number_of_subcurves() > 0,
                            "Cannot split a polycurve of length zero.");

      Comparison_result dir = cmp_seg_endpts(xcv[0]);

      // Locate the subcurve on the polycurve xcv that contains p.
      std::size_t i = m_traits.locate(xcv, p);
//      std::cerr << "SPLIT AT " << i << "(" << (dir == SMALLER ? "smaller" : "larger") << "): ";

      CGAL_precondition(i != Traits::INVALID_INDEX);

      /*
          i
      A B C D E F

      1 = A B C D (0 -> i+2)
      2 = D E F (i+1 -> end)

      1 = A B C (0 -> i+1)
      2 = C D E F (i-> end)
      */

      if (equal(max_vertex(xcv[i]), p)) {
//        std::cerr << "Split at " << p << " (" << i << "): " << xcv << std::endl;
        // The entire i'th subcurve belongs to xcv1:
        xcv1 = X_monotone_curve_2 (xcv.points_begin(), xcv[i+1]);
//        std::cerr << " -> " << xcv1 << std::endl;
        xcv2 = X_monotone_curve_2 (xcv[i], xcv.points_end());
//        std::cerr << " -> " << xcv2 << std::endl;
      }
      else if (equal(min_vertex(xcv[i]), p)) {
        // The entire i'th subcurves belongs to xcv2:
        xcv1 = X_monotone_curve_2 (xcv.points_begin(), xcv[i]);
        xcv2 = X_monotone_curve_2 (xcv[i-1], xcv.points_end());
      }
      else {
        // The i'th subcurve should be split: The left part(seg1)
        // goes to xcv1, and the right part(seg2) goes to xcv2.
        X_monotone_subcurve_2 seg1, seg2;
        auto p_ptr = std::make_shared<Point_2>(p);
        xcv1 = X_monotone_curve_2 (nullptr, xcv.points_begin(), xcv[i+1], p_ptr);
        xcv2 = X_monotone_curve_2 (p_ptr, xcv[i+1], xcv.points_end(), nullptr);
      }

      if (dir != SMALLER) std::swap(xcv1, xcv2);

      // std::cerr << "Splitted " << xcv << " at " << p << ", into "
      //           << xcv1 << " and " << xcv2 << std::endl;
    }

  };
  Split_2 split_2_object() const
  { return Split_2(*this); }

  class Intersect_2
  {
  protected:
    using Traits = Arr_lightweight_polyline_traits_2<Kernel, Point_iterator>;
    const Traits& m_traits;
    Intersect_2(const Traits& traits) : m_traits(traits) {}
    friend class Arr_lightweight_polyline_traits_2<Kernel, Point_iterator>;
  public:
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>        Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_subcurve_2>
                                                      Intersection_base_result;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                      Intersection_result;

      const Subcurve_traits_2* geom_traits = m_traits.subcurve_traits_2();
      auto cmp_y_at_x = m_traits.compare_y_at_x_2_object();
      auto equal = geom_traits->equal_2_object();
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto intersect = geom_traits->intersect_2_object();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();
      auto construct_opposite = geom_traits->construct_opposite_2_object();

      Comparison_result dir1 = cmp_seg_endpts(cv1[0]);
      Comparison_result dir2 = cmp_seg_endpts(cv2[0]);

      std::vector<X_monotone_subcurve_2> ocv; // Used to represent overlaps.
      const bool invert_ocv = ((dir1 == LARGER) && (dir2 == LARGER));
      const bool consistent = (dir1 == dir2);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      CGAL_assertion(consistent);
#endif

      const std::size_t n1 = cv1.number_of_subcurves();
      const std::size_t n2 = cv2.number_of_subcurves();

      std::size_t i1 = (dir1 == SMALLER) ? 0 : n1-1;
      std::size_t i2 = (dir2 == SMALLER) ? 0 : n2-1;

      auto compare_xy = m_traits.compare_xy_2_object();
      Comparison_result left_res =
        compare_xy(cv1[i1], ARR_MIN_END, cv2[i2], ARR_MIN_END);

      if (left_res == SMALLER) {
        // cv1's left endpoint is to the left of cv2's left endpoint:
        // Locate the index i1 of the subcurve in cv1 which contains cv2's
        // left endpoint.
        i1 = m_traits.locate_impl(cv1, cv2[i2], ARR_MIN_END,
                                       Are_all_sides_oblivious_tag());
        if (i1 == Traits::INVALID_INDEX) return oi;

        if (equal(max_vertex(cv1[i1]), min_vertex(cv2[i2]))) {
          if (((dir1 == SMALLER) && (i1 == n1-1)) ||
              ((dir1 == LARGER) && (i1 == 0))){
            // cv1's right endpoint equals cv2's left endpoint
            // Thus we can return this single(!) intersection point
            Intersection_point p(max_vertex(cv1[i1]), 0);
            *oi++ = Intersection_result(p);
            return oi;
          }
          dir1 == SMALLER ?
            ++i1 :
            (i1 != 0) ? --i1 : (std::size_t) Traits::INVALID_INDEX;
          left_res = EQUAL;
        }
      }
      else if (left_res == LARGER) {
        // cv1's left endpoint is to the right of cv2's left endpoint:
        // Locate the index i2 of the subcurve in cv2 which contains cv1's
        // left endpoint.
        i2 = m_traits.locate_impl(cv2, cv1[i1], ARR_MIN_END,
                                       Are_all_sides_oblivious_tag());
        if (i2 == Traits::INVALID_INDEX) return oi;

        if (equal(max_vertex(cv2[i2]), min_vertex(cv1[i1]))) {
          if (((dir2 == SMALLER) && (i2 == n2-1)) ||
              ((dir2 == LARGER) && (i2 == 0))){
            // cv2's right endpoint equals cv1's left endpoint
            // Thus we can return this single(!) intersection point
            Intersection_point p(max_vertex(cv2[i2]), 0);
            *oi++ = Intersection_result(p);
            return oi;
          }

          dir2 == SMALLER ?
            ++i2 :
            (i2 != 0) ? --i2 : (std::size_t) Traits::INVALID_INDEX;
          left_res = EQUAL;
        }
      }

      // Check if the the left endpoint lies on the other polycurve.
      bool left_coincides = (left_res == EQUAL);
      bool left_overlap = false;

      if (left_res == SMALLER)
        left_coincides = (cmp_y_at_x(cv2[i2], ARR_MIN_END, cv1[i1]) == EQUAL);
      else if (left_res == LARGER)
        left_coincides = (cmp_y_at_x(cv1[i1], ARR_MIN_END, cv2[i2]) == EQUAL);

      // The main loop: Go simultaneously over both polycurves.
      Comparison_result right_res = left_res;
      bool right_coincides = left_coincides;
      bool right_overlap = false;

      while (((dir1 == SMALLER) && (dir2 == SMALLER) &&
              (i1 < n1) && (i2 < n2)) ||
             ((dir1 != SMALLER) && (dir2 == SMALLER) &&
              (i1 != Traits::INVALID_INDEX) && (i2 < n2)) ||
             ((dir1 == SMALLER) && (dir2 != SMALLER) && (i1 < n1) &&
              (i2 != Traits::INVALID_INDEX)) ||
             ((dir1 != SMALLER) && (dir2 != SMALLER) &&
              (i1 != Traits::INVALID_INDEX) &&
              (i2 != Traits::INVALID_INDEX)))
      {
        right_res = compare_xy(cv1[i1], ARR_MAX_END, cv2[i2], ARR_MAX_END);

        right_coincides = (right_res == EQUAL);
        if (right_res == SMALLER)
          right_coincides =
            (cmp_y_at_x(cv1[i1], ARR_MAX_END, cv2[i2]) == EQUAL);
        else if (right_res == LARGER)
          right_coincides =
            (cmp_y_at_x(cv2[i2], ARR_MAX_END, cv1[i1]) == EQUAL);

        right_overlap = false;

        //! EF: the following code is abit suspicious. It may erroneously
        //      assume that the subcurves cannot overlap more than once.
        if (! right_coincides && ! left_coincides) {
          // Non of the endpoints of the current subcurve of one polycurve
          // coincides with the curent subcurve of the other polycurve:
          // Output the intersection if exists.
          CGAL_error_msg("Overlapping polyline not handled yet.");
          std::vector<Intersection_base_result> xections;
          intersect(cv1[i1], cv2[i2], std::back_inserter(xections));
          for (const auto& xection : xections) {
            const X_monotone_subcurve_2* subcv_p =
              boost::get<X_monotone_subcurve_2>(&xection);
            if (subcv_p != nullptr) {
              ocv.push_back(*subcv_p);
              oi = output_ocv (ocv, invert_ocv, oi);
              continue;
            }

            const Intersection_point* p_p =
              boost::get<Intersection_point>(&xection);
            if (p_p != nullptr) *oi++ = Intersection_result(*p_p);
          }
        }
        else if (right_coincides && left_coincides) {
          // An overlap exists between the current subcurves of the
          // polycurves: Output the overlapping subcurve.
          right_overlap = true;

          CGAL_error_msg("Overlapping polyline not handled yet.");
          std::vector<Intersection_base_result> sub_xections;
          intersect(cv1[i1], cv2[i2], std::back_inserter(sub_xections));

          for (const auto& item : sub_xections) {
            const X_monotone_subcurve_2* x_seg =
              boost::get<X_monotone_subcurve_2>(&item);
            if (x_seg != nullptr) {
              X_monotone_subcurve_2 seg = *x_seg;
              // We maintain the variant that if the input curves have opposite
              // directions (! consistent), the overalpping curves are directed
              // left=>right. This, however, is not guaranteed for the
              // subcurves. Therefore, we need to enforce it. That is, we make
              // sure the subcurves are also directed left=>right in this case.
              if (! consistent && (cmp_seg_endpts(seg) == LARGER))
                seg = construct_opposite(seg);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
              CGAL_assertion(cmp_seg_endpts(seg) == SMALLER);
#endif
              ocv.push_back(seg);
            }

            const Intersection_point* p_ptr =
              boost::get<Intersection_point>(&item);
            if (p_ptr != nullptr) {
              // Any point that is not equal to the max_vertex of the
              // subcurve should be inserted into oi.
              // The max_vertex of the current subcurve (if intersecting)
              // will be taken care of as the min_vertex of in the next
              // iteration.
              if (! equal(p_ptr->first, max_vertex(cv1[i1])))
                *oi++ = Intersection_result(*p_ptr);
            }
          }
        }

        else if (left_coincides && ! right_coincides) {
          // std::cout << "Left is coinciding but right is not." << std::endl;
          // The left point of the current subcurve of one polycurve
          // coincides with the current subcurve of the other polycurve.
          if (left_overlap) {
            // An overlap occurred at the previous iteration:
            // Output the overlapping polycurve.
            CGAL_assertion(ocv.size() > 0);
            oi = output_ocv (ocv, invert_ocv, oi);
          }
          else {
            // The left point of the current subcurve of one
            // polycurve coincides with the current subcurve of the
            // other polycurve, and no overlap occurred at the
            // previous iteration: Output the intersection
            // point. The derivative of at least one of the
            // polycurves is not defined at this point, so we give
            // it multiplicity 0.
            if (left_res == SMALLER) {
              Intersection_point p(min_vertex(cv2[i2]), 0);
              *oi++ = Intersection_result(p);
            }
            else {
              Intersection_point p(min_vertex(cv1[i1]), 0);
              *oi++ = Intersection_result(p);
            }
          }
        }

        // Proceed forward.
        if (right_res != SMALLER) {
          if (dir2 == SMALLER) ++i2;
          else {
            if (i2 == 0) i2 = Traits::INVALID_INDEX;
            else --i2;
          }
        }
        if (right_res != LARGER) {
          if (dir1 == SMALLER)
            ++i1;
          else {
            if (i1 == 0) i1 = Traits::INVALID_INDEX;
            else --i1;
          }
        }
        left_res = (right_res == SMALLER) ? LARGER :
          (right_res == LARGER) ? SMALLER : EQUAL;

        left_coincides = right_coincides;
        left_overlap = right_overlap;
      } // END of while loop

        // Output the remaining overlapping polycurve, if necessary.
      if (ocv.size() > 0) {
        oi = output_ocv (ocv, invert_ocv, oi);
      }
      else if (right_coincides) {
        typedef std::pair<Point_2,Multiplicity> return_point;
        return_point ip;
        if (right_res == SMALLER) {
          ip = (dir1 == SMALLER) ?
            return_point(max_vertex(cv1[i1-1]), 0) :
            (i1 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv1[i1+1]), 0) :
            return_point(max_vertex(cv1[0]), 0);
          *oi++ = Intersection_result(ip);
        }
        else if (right_res == LARGER) {
          ip = (dir2 == SMALLER) ?
            return_point(max_vertex(cv2[i2-1]), 0) :
            (i2 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv2[i2+1]), 0) :
            return_point(max_vertex(cv2[0]), 0);
          *oi++ = Intersection_result(ip);
        }
        else if (((i1 > 0) && (dir1 == SMALLER)) ||
                 ((i1 < n1) && (dir1 != SMALLER)) ||
                 ((i1 == Traits::INVALID_INDEX) &&
                  (dir1 != SMALLER)))
        {
          ip = (dir1 == SMALLER) ?
            return_point(max_vertex(cv1[i1-1]), 0) :
            (i1 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv1[i1+1]), 0) :
            return_point(max_vertex(cv1[0]), 0);
          *oi++ = Intersection_result(ip);
        }
        else {
          CGAL_assertion_msg((dir2 == SMALLER && i2 > 0) ||
                             (dir2 != SMALLER && i2 < n2) ||
                             (dir2 != SMALLER &&
                              ((i1 == Traits::INVALID_INDEX) ||
                               (i2 == Traits::INVALID_INDEX))),
                             "Wrong index for xcv2 in Intersect_2 of "
                             "polycurves.");
          ip = (dir2 == SMALLER) ?
            return_point(max_vertex(cv2[i2-1]), 0) :
            (i2 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv2[i2+1]), 0) :
            return_point(max_vertex(cv2[0]), 0);
          *oi++ = Intersection_result(ip);
        }
      }

      return oi;
    }

  private:
    template <typename OutputIterator>
    inline OutputIterator output_ocv
    (std::vector<X_monotone_subcurve_2>& ocv, bool invert_ocv, OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>        Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                      Intersection_result;
      X_monotone_curve_2 curve;
      if (invert_ocv)
        std::reverse (ocv.begin(), ocv.end());
      for (X_monotone_subcurve_2& sc : ocv)
        curve.push_back (sc);
      *(oi ++) = Intersection_result(curve);

      ocv.clear();

      return oi;
    }

  };
  Intersect_2 intersect_2_object() const
  { return Intersect_2(*this); }

  class Are_mergeable_2
  {
  public:
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      CGAL_assertion_msg(false, "Are_mergeable_2 not implemented");
      return true;
    }
  };
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(); }

  class Merge_2
  {
  public:
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_assertion_msg(false, "Merge_2 not implemented");
    }
  };
  Merge_2 merge_2_object() const
  { return Merge_2(); }
};


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
