// Copyright (c) 2018 GeometryFactory Sarl (France).
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

#ifndef CGAL_CLASSIFICATION_COMPRESSED_FLOAT_H
#define CGAL_CLASSIFICATION_COMPRESSED_FLOAT_H

#include <CGAL/license/Classification.h>

namespace CGAL {
namespace Classification {

/// \cond SKIP_IN_MANUAL

#if defined(CGAL_CLASSIFICATION_DO_NOT_USE_FLOATS_INTERNALLY)
typedef double internal_float;
#else
typedef float internal_float;
#endif

#if defined(CGAL_CLASSIFICATION_DO_NOT_COMPRESS_FLOATS)
typedef internal_float compressed_float;

inline internal_float compress_float (const internal_float& f, const internal_float& = 0.f, const internal_float& = 1.f)
{
  return f;
}

inline internal_float decompress_float (const internal_float& t, const internal_float& = 0.f, const internal_float& = 1.f)
{
  return t;
}

#else

#  if defined(CGAL_CLASSIFICATION_COMPRESS_FLOATS_WITH_USHORT)
typedef unsigned short compressed_float;
#  else // Default = compress with unsigned char
typedef unsigned char compressed_float;
#  endif

inline compressed_float compress_float (const internal_float& f, const internal_float& min = 0.f, const internal_float& max = 1.f)
{
  return static_cast<compressed_float>
    (internal_float(std::numeric_limits<compressed_float>::max()) * (f - min) / (max - min));
}

inline internal_float decompress_float (const compressed_float& t, const internal_float& min = 0.f, const internal_float& max = 1.f)
{
  return ((max - min) * (t / internal_float(std::numeric_limits<compressed_float>::max())) + min);
}
  
#endif


  /// \endcond

} // namespace Classification
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_COMPRESSED_FLOAT_H
