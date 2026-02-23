#pragma once

#include <vector>
#include <nlopt.hpp>
#include <iostream>
#include <utility>
#include <string>
#include <random>
#include "models/wpd.h"
#include "models/ntl.h"
#include "optimisation/WPD_opt.h"
#include "common/file_handler.h"


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

inline std::ostream& operator<<(std::ostream& stream, const std::vector<std::complex<double>>& vec)
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

inline std::pair<bool, std::string> save_results(const WPD::opt_setup& setup, const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3)
{
    bool saved = false;
    std::cout << "\n\n\n";
    std::string save;
    do
    {
        std::cout << "Save results? [Y/N]: ";
        std::cin >> save;
    } while (save != "y" && save != "Y" && save != "N" && save != "n");
    std::string folder;

    if (save == "Y" || save == "y")
    {
        do
        {
            try
            {
                std::cout << "Enter folder name: ";
                std::cin >> folder;
                if (folder == "0") break;
                double H = 4e-3;
                ::NTL::fh::wpd_to_file(wpd, folder + "/wpd");
                ::NTL::fh::ntl_to_file(out2, folder + "/out2");
                ::NTL::fh::ntl_to_file(out3, folder + "/out3");
                ::NTL::fh::setup_to_file<WPD::opt_setup>(setup, folder + "/setup");
                ::NTL::fh::export_geometry_scr(out2, H, folder + "/out2");
                ::NTL::fh::export_geometry_scr(out3, H, folder + "/out3");
                ::NTL::fh::export_geometry_scr(wpd.get_ntl2(), H, folder + "/ntl2");
                ::NTL::fh::export_geometry_scr(wpd.get_ntl3(), H, folder + "/ntl3");
                ::NTL::fh::export_geometry_csv(out2, H, folder + "/out2");
                ::NTL::fh::export_geometry_csv(out3, H, folder + "/out3");
                ::NTL::fh::export_geometry_csv(wpd.get_ntl2(), H, folder + "/ntl2");
                ::NTL::fh::export_geometry_csv(wpd.get_ntl3(), H, folder + "/ntl3");
                saved = true;
            }
            catch (std::exception& e)
            {
                std::cerr << "Saving ERROR: " << e.what() << ", check if file exists...\n";
            }
        } while (saved != true);
    }

    return { saved, folder };
}