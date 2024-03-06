#include <cmath>
#include "DRutils.h"

namespace DDDRCaloTubes
{
    int fast_floor(double x)
    {
        return (int) x - (x < (int) x);
    }

    int fast_ceil(double x)
    {
        return (int) x + (x > (int) x);
    }

    bool check_for_integer(double x)
    {
        return (std::abs(x - std::round(x)) < 1e-6);
    }
} // namespace DDDRCaloTubes
