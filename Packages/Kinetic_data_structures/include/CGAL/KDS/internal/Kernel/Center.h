// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_INTERNAL_CENTER_H
#define CGAL_KDS_INTERNAL_CENTER_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Center
{
    typedef typename K::Point_3 argument_type;
    typedef typename K::Point_3 result_type;

    const result_type& operator()(const argument_type &p) const
    {
        return p;
    }
    const result_type& operator()(const typename K::Weighted_point_3 &wp) const
    {
        return wp.point();
    }
};

CGAL_KDS_END_INTERNAL_NAMESPACE
#endif
