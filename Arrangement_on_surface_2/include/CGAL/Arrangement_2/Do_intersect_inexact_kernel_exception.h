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

struct Do_intersect_inexact_kernel_exception { };

template <typename FT, typename Point>
void throw_exception_if_new_point_in_inexact_kernel_impl
(std::pair<Point, unsigned int>& ip,
 const Point& p0l, const Point& p0r,
 const Point& p1l, const Point& p1r,
 const Tag_true&) // Is inexact
{
  if (p0l == p1l || p0l == p1r)
    ip = std::make_pair (p0l, 1);
  else if (p0r == p1l || p0r == p1r)
    ip = std::make_pair (p0r, 1);
  else
    throw Do_intersect_inexact_kernel_exception();
}

template <typename FT, typename Point>
void throw_exception_if_new_point_in_inexact_kernel_impl
(std::pair<Point, unsigned int>&,
 const Point&, const Point&,
 const Point&, const Point&,
 const Tag_false&) // Is exact
{ }

template <typename FT, typename Point>
void throw_exception_if_new_point_in_inexact_kernel
(std::pair<Point, unsigned int>& ip,
 const Point& p0l, const Point& p0r,
 const Point& p1l, const Point& p1r)
{
  throw_exception_if_new_point_in_inexact_kernel_impl<FT>
    (ip, p0l, p0r, p1l, p1r, Boolean_tag<std::is_floating_point<FT>::value>());
}

} // namespace CGAL


#endif // CGAL_ARR_DO_INTERSECT_INEXACT_KERNEL_EXCEPTION_H
