#pragma once

#include <vector>
#include <nlopt.hpp>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "models/ntl.h"
#include "common/helpers.h"

namespace NTL
{
	class optimiser;
	struct optimiser_setup;
	struct optimiser_result;
	enum class console;


	struct optimiser_setup
	{
		optimiser_setup() {};
		optimiser_setup(const nlohmann::json& j);
		virtual ~optimiser_setup() = default;

		int N{ 0 };
		std::vector<double> lb{};
		std::vector<double> ub{};
		std::vector<double> toll_bounds{};
		std::vector<double> toll_z{};
		double GBL_MAX{ 0 };
		double LCL_MAX{ 0 };
		double accepted_error{ 0 };
		int max_attempts{ 10'000 };

		virtual nlohmann::json get_json() const;
	};

	struct optimiser_result
	{
		std::vector<double> optimised_cn;
		double final_error;
		int number_of_attempts;
	};

	enum class console
	{
		active = true, inactive = false
	};

	class optimiser
	{
	public:
		optimiser(const optimiser_setup& setup);
			

	protected:
		optimiser_result run_optimiser(console mode = console::inactive);
		double m_accepted_error;
		int m_max_attempts;
		int m_N;
		std::vector<double> m_lb;
		std::vector<double> m_ub;
		std::vector<double> m_toll_bounds;
		std::vector<double> m_toll_z;
		double m_GBL_MAX;
		double m_LCL_MAX;	
		

	private:
		//nlopt functions
		virtual double min_objective(const std::vector<double>& Cn) const = 0;
		virtual void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const = 0;

	};
}