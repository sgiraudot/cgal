// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// $Date$
//
// Author(s): Simon Giraudot <simon.giraudot@geometryfactory.com>

#ifndef CGAL_ARR_DO_INTERSECT_INEXACT_KERNEL_EXCEPTION_H
#define CGAL_ARR_DO_INTERSECT_INEXACT_KERNEL_EXCEPTION_H

#include <CGAL/license/Arrangement_on_surface_2.h>

namespace CGAL
{

struct Do_intersect_inexact_kernel_exception : public std::logic_error
{
  Do_intersect_inexact_kernel_exception()
    : std::logic_error("Segment/Segment intersection computed using inexact constructions.")
  { }
};

template <typename FT, typename Point>
bool throw_exception_if_new_point_in_inexact_kernel_impl
(std::pair<Point, unsigned int>& ip,
 const Point& p0l, const Point& p0r,
 const Point& p1l, const Point& p1r,
 const Tag_true&) // Is inexact
{
  // Check if segment/segment intersection its a simple
  // vertex-to-vertex case, in which case we replace the content of ip
  // by the common point to avoid using inexact constructions.
  //
  // If the point is a new point, then we throw an exception: by
  // default, users will get the exception and be advised to use an
  // exact kernel, except for `do_intersect` where the exception is
  // caught and the function returns `true` (as an intersection was
  // detected).
  bool p0l_inter = (p0l == p1l || p0l == p1r);
  bool p0r_inter = (p0r == p1l || p0r == p1r);

  if (p0l_inter && !p0r_inter)
  {
    ip = std::make_pair (p0l, 1);
    return false;
  }
  else if (!p0l_inter && p0r_inter)
  {
    ip = std::make_pair (p0r, 1);
    return false;
  }
  else if (p0l_inter && p0r_inter)
    return true;

  throw Do_intersect_inexact_kernel_exception();
  return true; // never reached
}

template <typename FT, typename Point>
bool throw_exception_if_new_point_in_inexact_kernel_impl
(std::pair<Point, unsigned int>&,
 const Point&, const Point&,
 const Point&, const Point&,
 const Tag_false&) // Is exact
{
  // Nothing to do in exact
  return true;
}

template <typename FT, typename Point>
bool throw_exception_if_new_point_in_inexact_kernel
(std::pair<Point, unsigned int>& ip,
 const Point& p0l, const Point& p0r,
 const Point& p1l, const Point& p1r)
{
  return throw_exception_if_new_point_in_inexact_kernel_impl<FT>
    (ip, p0l, p0r, p1l, p1r, Boolean_tag<std::is_floating_point<FT>::value>());
}

} // namespace CGAL


#endif // CGAL_ARR_DO_INTERSECT_INEXACT_KERNEL_EXCEPTION_H
