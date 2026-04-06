#include "models/ntl.h"
#include "models/wpd.h"
#include "models/tj.h"
#include "optimisation/NTL_opt.h"
#include "optimisation/WPD_opt.h"
#include "optimisation/TJ_opt.h"

#include <limits.h>
#include <algorithm>
#include <vector>
#include <utility>

inline std::pair<double, double> find_Z_range(const NTL::NTL& ntl)
{
    auto Z = ntl.get_Z_vec(1e-5);

    if (Z.empty())
    {
        return { 0.0, 0.0 };
    }

    auto min_max = std::minmax_element(Z.begin(), Z.end());

    return ( *min_max.first, *min_max.second );
}

inline std::pair<double, double> find_Z_range(const WPD::WPD& wpd)
{
    auto range2 = find_Z_range(wpd.get_ntl2());
    auto range3 = find_Z_range(wpd.get_ntl3());

    return { (range2.first < range3.first) ? range2.first : range3.first,
             (range2.second > range3.second) ? range2.second : range3.second };
}

inline std::pair<double, double> find_Z_range(const TJ::TJ& tj)
{
    auto range2 = find_Z_range(tj.get_ntl2());
    auto range3 = find_Z_range(tj.get_ntl3());

    return { (range2.first < range3.first) ? range2.first : range3.first,
             (range2.second > range3.second) ? range2.second : range3.second };
}

inline bool is_Z_range_bound(const std::pair<double, double>& Z_min_max, double Z_min, double Z_max)
{
    if (Z_min_max.first == 0 && Z_min_max.second == 0)
        return false;

    return (Z_min_max.first >= Z_min && Z_min_max.second <= Z_max);
}

inline bool validate(const NTL::NTL& ntl)
{

    if (ntl.get_Cn().empty())
    {
        std::cout << "NTL not solved yet\n";
        return false;
    }

    if (ntl.get_Z0() <= 1e-12)
    {
        std::cout << "Z0 not yet set\n";
        return false;
    }

    if (ntl.get_er() <= 1e-12)
    {
        std::cout << "Substrate not yet set\n";
        return false;
    }

    if (ntl.get_d() <= 1e-12)
    {
        std::cout << "NTL size equal to zero\n";
        return false;
    }

    return true;
}

inline bool validate(const WPD::WPD& wpd)
{
    if (!validate(wpd.get_ntl2()))
    {
        std::cout << "Invalid NTL2\n";
        return false;
    }

    if (!validate(wpd.get_ntl3()))
    {
        std::cout << "Invalid NTL3\n";
        return false;
    }

    if (wpd.get_R() <= 1e-12)
    {
        std::cout << "Invalid R\n";
        return false;
    }

    return true;
}

inline bool validate(const TJ::TJ& tj)
{
    if (!validate(tj.get_ntl2()))
    {
        std::cout << "Invalid NTL2\n";
        return false;
    }

    if (!validate(tj.get_ntl3()))
    {
        std::cout << "Invalid NTL3\n";
        return false;
    }

    return true;
}

inline bool validate_result(const NTL::opt_setup& setup, const NTL::opt_result result)
{
    if (!validate(result.ntl))
    {
        std::cout << "Invalid NTL\n";
        return false;
    }

    if (!is_Z_range_bound(find_Z_range(result.ntl), setup.Z_min, setup.Z_max))
    {
        std::cout << "Z(z) not bounded by constraints properly\n";
        return false;
    }

    return true;
}

inline bool validate_result(const WPD::opt_setup& setup, const WPD::opt_result result)
{
    if (!validate(result.wpd))
    {
        std::cout << "Invalid WPD\n";
        return false;
    }

    if (!validate(result.output2))
    {
        std::cout << "Invalid Out2\n";
        return false;
    }

    if (!validate(result.output3))
    {
        std::cout << "Invalid Out3\n";
        return false;
    }

    if (!is_Z_range_bound(find_Z_range(result.wpd), setup.Z_min, setup.Z_max))
    {
        std::cout << "Z(z) not bounded by constraints properly\n";
        return false;
    }

    if (!is_Z_range_bound(find_Z_range(result.output2), setup.Z_min, setup.Z_max))
    {
        std::cout << "Z(z) not bounded by constraints properly\n";
        return false;
    }

    if (!is_Z_range_bound(find_Z_range(result.output3), setup.Z_min, setup.Z_max))
    {
        std::cout << "Z(z) not bounded by constraints properly\n";
        return false;
    }

    return true;
}

inline bool validate_result(const TJ::opt_setup& setup, const TJ::opt_result result)
{
    if (!validate(result.tj))
    {
        std::cout << "Invalid TJ\n";
        return false;
    }

    if (!is_Z_range_bound(find_Z_range(result.tj), setup.Z_min, setup.Z_max))
    {
        std::cout << "Z(z) not bounded by constraints properly\n";
        return false;
    }

    return true;
}