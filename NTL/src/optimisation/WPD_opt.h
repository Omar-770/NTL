#pragma once

#include "optimiser.h"
#include "models/wpd.h"
#include "models/ntl.h"
#include <omp.h>
#include <iostream>
#include <nlohmann/json.hpp>
#include <chrono>
#include <common/file_handler.h>

namespace WPD
{
	class WPD_opt;
	struct WPD_opt_setup;
	struct WPD_opt_result;

	using NTL::opt;
	using NTL::opt_setup;
	using NTL::opt_result;
	using NTL::console;

	struct WPD_opt_setup : public opt_setup
	{
		WPD_opt_setup() {};
		WPD_opt_setup(const WPD_opt_setup& setup);
		WPD_opt_setup(const nlohmann::json& j);

		double Z0{0};
		double er{0};
		double d{ 0 };
		double Zref{ 0 };
		NTL::NTL ntl2;
		NTL::NTL ntl3; 
		std::vector<double> freqs{};
		int K{ 0 };
		double Z_min{ 0 };
		double Z_max{ 0 };
		double R_min{ 0 };
		double R_max{ 0 };
		double matching_dB{ 0 };
		double isolation_dB{0};
		std::vector<double> split{};
		
		nlohmann::json get_json() const override;
	};

	struct WPD_opt_result : public opt_result
	{
		WPD wpd;
	};

	class WPD_opt : public opt
	{
	public:
		WPD_opt(const WPD_opt_setup& setup);

		void optimise(console mode = console::active) {};


	private:
		double m_Z0;
		double m_er;
		double m_d;
		double m_Zref; //port1 impedance and load impedance for output transformers
		NTL::NTL m_ntl2; //port2 output transformer
		NTL::NTL m_ntl3; //port3 output transformer
		std::vector<double> m_freqs;
		int m_K;
		double m_Z_min;
		double m_Z_max;
		double m_R_min;
		double m_R_max;
		double m_matching_dB;
		double m_isolation_dB;
		std::vector<double> m_split; //P3/P2
		
	private:
		std::vector<std::array<std::complex<double>, 3>> m_Zl;
		int m_N2;
		int m_N3;
		double m_matching_abs;
		double m_isolation_abs;
		std::vector <std::array<double, 2>> m_split_abs;

	private:
		//nlopt functions
		double min_objective(const std::vector<double>& Cn) const override;
		void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const override;
		double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const override;
	};
}