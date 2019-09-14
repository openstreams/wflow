// -----------------------------------------------------------------------------
// Fern © Geoneric
//
// This file is part of Geoneric Fern which is available under the terms of
// the GNU General Public License (GPL), version 2. If you do not want to
// be bound by the terms of the GPL, you may purchase a proprietary license
// from Geoneric (http://www.geoneric.eu/contact).
// -----------------------------------------------------------------------------
#pragma once
#include <utility>
#include "fern/example/algorithm/operation.h"


namespace example {

template<
    typename Value1,
    typename Value2>
fern::algorithm::add::result_type<Value1, Value2> operator+(
    Value1 const& lhs,
    Value2 const& rhs)
{
    return add(lhs, rhs);
}

} // namespace example
