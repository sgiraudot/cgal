#ifndef CGAL_LIGHTWEIGHT_POLYLINE_2_H
#define CGAL_LIGHTWEIGHT_POLYLINE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/iterator.h>

#include <boost/iterator/iterator_facade.hpp>

namespace CGAL {

namespace internal {

#if 0 // Draft idea if isolated points on the sides are not enough
template <typename Iterator>
class Augmented_range
{
public:

  using Base_iterator = Iterator;
  using value_type = typename std::iterator_traits<Base_iterator>::value_type;

private:

  std::vector<value_type> m_front;
  Base_iterator m_begin;
  Base_iterator m_end;
  std::vector<value_type> m_back;

public:

};
#endif



template <typename Kernel_, typename PointIterator>
class Lightweight_polyline_2_iterator;

template <typename Kernel_, typename PointIterator>
class Lightweight_polyline_2
{
public:

  using Kernel = Kernel_;
  using Base_iterator = PointIterator;
  using Self = Lightweight_polyline_2<Kernel, Base_iterator>;

  using Point_type_2 = typename Kernel::Point_2;
  using Point_ptr = std::shared_ptr<Point_type_2>;
  using Size = std::size_t;
  using size_type = std::size_t;

  using iterator = Lightweight_polyline_2_iterator<Kernel_, Base_iterator>;
  friend iterator;

  using Subcurve_type_2 = iterator;
  using Subcurve_iterator = Prevent_deref<iterator>;
  using Subcurve_const_iterator = Prevent_deref<iterator>;

protected:

  Base_iterator m_begin;
  Base_iterator m_end;
  Point_ptr m_first;
  Point_ptr m_last;

public:

  Lightweight_polyline_2() { }

#if 0
  Lightweight_polyline_2(const Subcurve_type_2&) { }
#endif

  Lightweight_polyline_2 (Base_iterator begin, Base_iterator end)
    : m_begin(begin), m_end(end)
  {
    CGAL_assertion (std::distance (begin, end) >= 2);
  }

  Lightweight_polyline_2 (iterator begin, iterator end)
  {
    const Self& support = begin.support();
    CGAL_assertion (&support == &end.support());

//    CGAL_assertion (std::distance(begin, end) >= 2);
    if (std::distance(begin, end) < 2) // Polyline with less than 2 points is empty
      return;

    if (begin.base() == support.m_begin - 1)
    {
      CGAL_assertion (support.m_first != nullptr);
      m_first = support.m_first;
      m_begin = support.m_begin;
    }
    else
      m_begin = begin.base();

    if (end.base() == support.m_end + 1)
    {
      CGAL_assertion (support.m_last != nullptr);
      m_last = support.m_last;
      m_end = support.m_end;
    }
    else
      m_end = end.base();

    // This constructor should only create x-monotone polylines
    CGAL_assertion (is_x_monotone());

//    std::cerr << "New polyline from 2 iterators " << *this << std::endl;
  }

  Lightweight_polyline_2 (Point_ptr first, iterator begin, iterator end, Point_ptr last)
  {
#if 0
    if (first) std::cerr << *first;
    else std::cerr << "nullptr";

    for (iterator it = begin; it != end; ++ it)
      std::cerr << " " << *it;

    if (last) std::cerr << " " << *last << std::endl;
    else std::cerr << " nullptr" << std::endl;
#endif

    const Self& support = begin.support();
    CGAL_assertion (&support == &end.support());

    if (begin.base() == support.m_begin - 1)
    {
      CGAL_assertion (first == nullptr);
      m_first = support.m_first;
      m_begin = support.m_begin;
    }
    else
    {
      m_first = first;
      m_begin = begin.base();
    }

    if (end.base() == support.m_end + 1)
    {
      CGAL_assertion (last == nullptr);
      m_last = support.m_last;
      m_end = support.m_end;
    }
    else
    {
      m_last = last;
      m_end = end.base();
    }

    // This constructor should only create x-monotone polylines
    CGAL_assertion (is_x_monotone());
//    std::cerr << "New polyline with hanging points " << *this << std::endl;
  }

  bool is_x_monotone() const
  {
    auto compare_x_2 = Kernel().compare_x_2_object();
    iterator b = points_begin(), e = points_end() - 1;
    Comparison_result comp = Kernel().compare_x_2_object()(*b, *(b+1));
    for (iterator it = b + 1; it < e; ++ it)
      if (comp != Kernel().compare_x_2_object()(*it, *(it+1)))
        return false;
    return true;
  }

  inline void push_back(const Subcurve_type_2& seg)
  {
    CGAL_error_msg("push_back not available");
  }

  inline void push_front(const Subcurve_type_2& seg)
  {
    CGAL_error_msg("push_front not available");
  }

  Bbox_2 bbox() const
  {
    return bbox_2 (points_begin(), points_end());
  }

  iterator points_begin() const { return iterator (this, Tag_true()); }
  iterator points_end() const { return iterator (this, Tag_false()); }

  Subcurve_const_iterator subcurves_begin() const { return Subcurve_const_iterator(points_begin()); }
  Subcurve_const_iterator subcurves_end() const { return Subcurve_const_iterator(points_end() - 1); }

  size_type number_of_subcurves() const { return std::distance (points_begin(), points_end()) - 1; }

  inline Subcurve_type_2 operator[](const std::size_t i) const
  {
    return points_begin() + i;
  }

  inline void clear()
  {
    m_begin = Base_iterator();
    m_end = Base_iterator();
    m_first = nullptr;
    m_last = nullptr;
  }

  friend std::ostream& operator<< (std::ostream& os, const Self& p)
  {
    os << p.number_of_subcurves();
    for (iterator it = p.points_begin(); it != p.points_end(); ++ it)
      os << " " << *it;
    return os;
  }


};

template <typename Kernel_, typename PointIterator>
class Lightweight_polyline_2_iterator
  : public boost::iterator_facade<Lightweight_polyline_2_iterator<Kernel_, PointIterator>,
                                  typename Kernel_::Point_2,
                                  std::random_access_iterator_tag>
{
public:

  using Kernel = Kernel_;
  using Base_iterator = PointIterator;
  using Self = Lightweight_polyline_2_iterator<Kernel, Base_iterator>;
  using Polyline = Lightweight_polyline_2<Kernel, Base_iterator>;
  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;
  friend Polyline;

protected:
  const Polyline* m_support;
  Base_iterator m_base;
  bool m_reverse;

public:

  Lightweight_polyline_2_iterator () : m_reverse(false) { }
  Lightweight_polyline_2_iterator (const Polyline* support, const Tag_true&) // begin
    : m_support (support), m_reverse(false)
  {
    if (m_support->m_first)
      m_base = m_support->m_begin - 1;
    else
      m_base = m_support->m_begin;
  }

  Lightweight_polyline_2_iterator (const Polyline* support, const Tag_false&) // end
    : m_support (support), m_reverse(false)
  {
    if (m_support->m_last)
      m_base = m_support->m_end + 1;
    else
      m_base = m_support->m_end;
  }

  const Polyline& support() const { return *m_support; }
  Base_iterator base() const { return m_base; }

  // The iterator is also used as a wrapper for a segment (using it as
  // source and using the immediate following iterator as target)

  const Point_2& source() const { return const_dereference(); }
  const Point_2& target() const { Self copy(*this); ++ copy; return copy.const_dereference(); }

  bool is_vertical() const
  {
    const Point_2& ps = source();
    const Point_2& pt = target();

    return (Kernel().compare_x_2_object()(ps, pt) == EQUAL &&
            Kernel().compare_y_2_object()(ps, pt) != EQUAL);
  }

  bool is_directed_right() const
  { return Kernel().compare_xy_2_object()(source(), target()) == SMALLER; }

  const Point_2& left() const
  {
    if (is_directed_right())
      return source();
    // else
    return target();
  }

  const Point_2& right() const
  {
    if (is_directed_right())
      return target();
    // else
    return source();
  }

  Self opposite() const
  {
    Self copy(*this);
    ++ copy;
    copy.m_reverse = true;
  }

  Line_2 line() const { return Line_2 (source(), target()); }

private:
  friend class boost::iterator_core_access;

  void increment()
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (m_base != Base_iterator());
    if (m_reverse)
      -- m_base;
    else
      ++ m_base;
  }

  void decrement()
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (m_base != Base_iterator());
    if (m_reverse)
      ++ m_base;
    else
      -- m_base;
  }

  void advance (std::ptrdiff_t n)
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (m_base != Base_iterator());
    if (m_reverse)
      m_base -= n;
    else
      m_base += n;
  }

  std::ptrdiff_t distance_to (const Self& other) const
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (other.m_support != nullptr);
    if (m_reverse)
      return std::distance (other.m_base, m_base);
    // else
      return std::distance (m_base, other.m_base);
  }

  bool equal (const Self& other) const
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (other.m_support != nullptr);
    return m_base == other.m_base;
  }

  Point_2& dereference() const
  {
    return const_cast<Point_2&>(const_dereference());
  }

  const Point_2& const_dereference() const
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (m_base != Base_iterator());
    if (m_base == m_support->m_begin - 1)
    {
      CGAL_assertion (m_support->m_first != nullptr);
      return *(m_support->m_first);
    }
    else if (m_base == m_support->m_end)
    {
      CGAL_assertion (m_support->m_last != nullptr);
      return *(m_support->m_last);
    }

    // else
    CGAL_assertion (std::distance (m_support->m_begin, m_base) >= 0);
    CGAL_assertion (std::distance (m_base, m_support->m_end) > 0);
    return *m_base;
  }
};

} // namespace internal

} //namespace CGAL

#endif
