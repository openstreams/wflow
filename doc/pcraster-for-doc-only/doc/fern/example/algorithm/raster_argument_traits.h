// -----------------------------------------------------------------------------
// Fern © Geoneric
//
// This file is part of Geoneric Fern which is available under the terms of
// the GNU General Public License (GPL), version 2. If you do not want to
// be bound by the terms of the GPL, you may purchase a proprietary license
// from Geoneric (http://www.geoneric.eu/contact).
// -----------------------------------------------------------------------------
#pragma once
#include "fern/algorithm/core/argument_traits.h"
#include "fern/algorithm/policy/skip_no_data.h"
#include "fern/algorithm/policy/dont_mark_no_data.h"
#include "fern/example/algorithm/raster.h"


namespace fern {
namespace algorithm {

template<
    typename T>
struct ArgumentTraits<
    example::Raster<T>>
{

    using Mask = example::Raster<T>;

    using InputNoDataPolicy = algorithm::SkipNoData;

    using OutputNoDataPolicy = algorithm::DontMarkNoData;

};

} // namespace algorithm
} // namespace fern
