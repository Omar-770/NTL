#pragma once

#include <vector>
#include <nlopt.hpp>
#include <random>


using complex = std::complex<double>;

inline std::vector<double> rand_start(const nlopt::opt& opt)
{
    std::vector<double> start(opt.get_dimension());
    const std::vector<double>& lb = opt.get_lower_bounds();
    const std::vector<double>& ub = opt.get_upper_bounds();

    static thread_local std::mt19937_64 gen(std::random_device{}());

    for (int i = 0; i < opt.get_dimension(); i++)
    {
        std::uniform_real_distribution<double> dist(lb[i], ub[i]);
        start[i] = dist(gen);
    }

    return start;
}

inline double rand_start(double min, double max)
{
    static thread_local std::mt19937_64 gen(std::random_device{}());
    double start;

    std::uniform_real_distribution<double> dist(min, max);
    start = dist(gen);


    return start;
}


inline std::ostream& operator<<(std::ostream& stream, const std::vector<double>& vec)
{
    stream << '{';
    for (int i = 0; i < vec.size(); i++)
    {
        stream << vec[i];
        if (i < vec.size() - 1)
            stream << ", ";
    }
    stream << '}';

    return stream;
}