#include "latex.h"
#include <fstream>

bool ntl_profiles(const NTL::NTL& ntl, const std::string& folder, const std::string& file)
{
    std::filesystem::path dir = std::filesystem::path(folder);
    std::filesystem::create_directory(dir); 

    std::ofstream profile(dir / (file + "_profile.csv"));

    if (profile.is_open())
    {
        double step = 1e-4;
        const auto& z_profile = ntl.get_Z_vec(step);
        const auto& w_h_profile = ntl.get_w_h_vec(step);
        
        profile << "Z_mm,Z_Ohm,W_H,W_H_plus,W_H_minus" << std::endl;
        
        for (int i = 0; i < z_profile.size(); i++)
        {
            profile << z_profile[i].first * 1e3 << ',';
            profile << z_profile[i].second << ',';
            profile << w_h_profile[i].second << ',';
            profile << w_h_profile[i].second / 2.0 << ',';
            profile << w_h_profile[i].second / -2.0 << '\n';
        }
    }
    else
    {
        return false;
    }

    return true;
}

bool wpd_profiles(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3, const std::string& folder, const std::string& file)
{
    if (!ntl_profiles(wpd.get_ntl2(), folder, file + "_arm2"))   return false;
    if (!ntl_profiles(wpd.get_ntl3(), folder, file + "_arm3"))   return false;
    if (!ntl_profiles(out2, folder, file + "_out2"))   return false;
    if (!ntl_profiles(out3, folder, file + "_out3"))   return false;


    return true;
}   

bool wpd_sparams(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3,
    double Zref, double f_min, double f_max, double f_step, const std::vector<double>& freqs,
    const std::string& folder, const std::string& file)
{
    std::filesystem::path dir = std::filesystem::path(folder);
    std::filesystem::create_directory(dir);

    std::ofstream sp_file(dir / (file + "_sparams_code.csv"));
    std::ofstream f_targets(dir / (file + "_targets_code.csv"));

    if (!sp_file.is_open()) return false;

    // Header: Freq in GHz, S-parameters in dB
    sp_file << "Freq_GHz,S11,S22,S33,S21,S31,S23\n";

    auto to_dB = [](std::complex<double> val) {
        return 20.0 * std::log10(std::max(std::abs(val), 1e-15)); // Prevent log(0)
        };

    for (double f = f_min; f <= f_max + (f_step * 0.1); f += f_step)
    {
        // 50 integration steps is standard in your sim engine
        auto S = wpd.system_S_matrix(Zref, out2, out3, f, 50);

        sp_file << (f / 1e9) << ','
            << to_dB(S(0, 0)) << ','
            << to_dB(S(1, 1)) << ','
            << to_dB(S(2, 2)) << ','
            << to_dB(S(1, 0)) << ','  // S21
            << to_dB(S(2, 0)) << ','  // S31
            << to_dB(S(1, 2)) << '\n'; // S23
    }
    
    if (f_targets.is_open()) 
    {
        f_targets << "Freq_GHz,S11,S22,S33,S12,S13,S23\n";
        for (size_t i = 0; i < freqs.size(); ++i) 
        {
            double f = freqs[i];
            auto S = wpd.system_S_matrix(Zref, out2, out3, f, 50);
            f_targets << (f / 1e9) << ","
                << to_dB(S(0, 0)) << ","
                << to_dB(S(1, 1)) << ","
                << to_dB(S(2, 2)) << ","
                << to_dB(S(0, 1)) << ","
                << to_dB(S(0, 2)) << ","
                << to_dB(S(1, 2)) << "\n";
        }
        f_targets.close();
    }

    return true;
}

bool wpd_table_data(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3,
    const std::string& folder, const std::string& file)
{
    std::filesystem::path dir = std::filesystem::path(folder);
    std::filesystem::create_directory(dir);

    // 1. Export the Resistor and basic scalars
    std::ofstream f_scal(dir / (file + "_scalars.csv"));
    if (f_scal.is_open()) {
        f_scal << "Parameter,Value\n";
        f_scal << "R," << wpd.get_R() << "\n";
        f_scal << "Z0," << wpd.get_Z0() << "\n";
        f_scal.close();
    }
    else return false;

    // 2. Export the Coefficients
    std::ofstream f_coef(dir / (file + "_coeffs.csv"));
    if (f_coef.is_open()) {
        f_coef << "Term,Arm2,Arm3,Out2,Out3\n";

        auto ntl2 = wpd.get_ntl2();
        auto ntl3 = wpd.get_ntl3();

        // Helper to extract a specific term (C or S) from an NTL safely
        auto get_coeff = [](const NTL::NTL& ntl, char type, int index) -> std::string {
            const auto& cn = ntl.get_Cn();
            if (cn.empty()) return "-";

            int M = ntl.get_M();
            int N_cos = cn.size() - M;

            if (type == 'C') {
                if (index < N_cos) return std::to_string(cn[index]);
            }
            else if (type == 'S') {
                if (index > 0 && index <= M) return std::to_string(cn[N_cos + index - 1]);
            }
            return "-";
            };

        // Find the maximum number of Cosine and Sine terms across all 4 lines
        int max_C = std::max({
            (int)ntl2.get_Cn().size() - ntl2.get_M(),
            (int)ntl3.get_Cn().size() - ntl3.get_M(),
            (int)out2.get_Cn().size() - out2.get_M(),
            (int)out3.get_Cn().size() - out3.get_M()
            }) - 1; // -1 because C is 0-indexed

        int max_S = std::max({
            ntl2.get_M(), ntl3.get_M(), out2.get_M(), out3.get_M()
            });

        // Write Cosine terms (C_0, C_1, ...)
        for (int i = 0; i <= max_C; ++i) {
            f_coef << "C_{" << i << "},";
            f_coef << get_coeff(ntl2, 'C', i) << ",";
            f_coef << get_coeff(ntl3, 'C', i) << ",";
            f_coef << get_coeff(out2, 'C', i) << ",";
            f_coef << get_coeff(out3, 'C', i) << "\n";
        }

        // Write Sine terms (S_1, S_2, ...)
        for (int i = 1; i <= max_S; ++i) {
            f_coef << "S_{" << i << "},";
            f_coef << get_coeff(ntl2, 'S', i) << ",";
            f_coef << get_coeff(ntl3, 'S', i) << ",";
            f_coef << get_coeff(out2, 'S', i) << ",";
            f_coef << get_coeff(out3, 'S', i) << "\n";
        }

        f_coef << "d (mm),"
            << ntl2.get_d() * 1000.0 << ","
            << ntl3.get_d() * 1000.0 << ","
            << out2.get_d() * 1000.0 << ","
            << out3.get_d() * 1000.0 << "\n";

        f_coef.close();
    }
    else return false;

    return true;
}