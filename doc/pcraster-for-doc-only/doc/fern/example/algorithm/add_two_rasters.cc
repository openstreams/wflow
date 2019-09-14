// -----------------------------------------------------------------------------
// Fern © Geoneric
//
// This file is part of Geoneric Fern which is available under the terms of
// the GNU General Public License (GPL), version 2. If you do not want to
// be bound by the terms of the GPL, you may purchase a proprietary license
// from Geoneric (http://www.geoneric.eu/contact).
// -----------------------------------------------------------------------------
// Example showing how to make it possible to pass instances of user-defined
// classes, unknown to fern, to algorithms.
// - User-defined raster class: example::Raster
// - Glue code to allow passing rasters to algorithm: example::DataTypeTraits
// - Operators to support nice syntax.
#include <cstdlib>
#include <numeric>
#include "fern/example/algorithm/operation.h"
#include "fern/example/algorithm/operator.h"


namespace example {

// Execution policy to use.
fern::algorithm::ExecutionPolicy execution_policy =
    fern::algorithm::ParallelExecutionPolicy{};

}  // namespace example


int main(
    int /* argc */,
    char** /* argv */)
{
    using namespace example;

    double const cell_size{5.0};
    size_t const nr_rows{6000};
    size_t const nr_cols{4000};

    // [0, 1, 2, 3, ...]
    Raster<int32_t> raster1(cell_size, nr_rows, nr_cols);
    std::iota(begin(raster1), end(raster1), 0);

    // [5, 5, 5, ...]
    Raster<int32_t> raster2(cell_size, nr_rows, nr_cols);
    std::fill(begin(raster2), end(raster2), 5);

    // [5, 6, 7, 8, ...]
    // Operator syntax.
    auto raster3 = raster1 + raster2;
    assert(raster3.get(100) == 100 + 5);

    // Function call syntax.
    auto raster4 = slope(cast<double>(raster3));

    return EXIT_SUCCESS;
}
