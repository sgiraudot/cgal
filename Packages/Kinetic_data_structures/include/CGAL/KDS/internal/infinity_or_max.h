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

#ifndef CGAL_KDS_INFINITY_OR_MAX_H
#define CGAL_KDS_INFINITY_OR_MAX_H

#include <limits>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
template <class T>
T infinity_or_max()
{
    if (std::numeric_limits<T>::has_infinity) return std::numeric_limits<T>::infinity();
    else return std::numeric_limits<T>::max();
}


template <class T>
T infinity_or_max(T)
{
    if (std::numeric_limits<T>::has_infinity) return std::numeric_limits<T>::infinity();
    else return std::numeric_limits<T>::max();
}


CGAL_KDS_END_INTERNAL_NAMESPACE
#endif
