#include "NTL_opt.h"

namespace NTL
{
	NTL_opt_setup::NTL_opt_setup(const nlohmann::json& j) : opt_setup(j)
	{
		if (j.at("setup_type") != "NTL_opt")
			throw(std::logic_error("Attempted to read a setup from a different json object"));

		Z0 = j.at("Z0").get<double>();
		er = j.at("er").get<double>();
		d = j.at("d").get<double>();
		Zs = j.at("Zs").get<double>();
		Zl = j.at("Zl").get<std::vector<double>>();
		freqs = j.at("freqs").get<std::vector<double>>();
		K = j.at("K").get<double>();
		Z_min = j.at("Z_min").get<double>();
		Z_max = j.at("Z_max").get<double>();
	}

	nlohmann::json NTL_opt_setup::get_json() const
	{
		return {
			{ "json_type", "setup" },
			{ "setup_type", "NTL_opt"},
			{ "N", N },
			{ "lb", lb },
			{ "ub", ub },
			{ "toll_bounds", toll_bounds },
			{ "toll_z", toll_z },
			{ "GBL_MAX", GBL_MAX },
			{ "LCL_MAX", LCL_MAX },
			{ "accepted_error", accepted_error },
			{ "max_attempts", max_attempts },
			{ "Z0", Z0 },
			{ "er", er },
			{ "d", d },
			{ "Zs", Zs },
			{ "Zl", Zl },
			{ "freqs", freqs },
			{ "K", K },
			{ "Z_min", Z_min },
			{ "Z_max", Z_max }
		};
	}

	NTL_opt::NTL_opt(const NTL_opt_setup& setup) : opt(setup),
		m_Z0(setup.Z0), m_er(setup.er), m_d(setup.d), m_Zs(setup.Zs), m_Zl(setup.Zl),
		m_freqs(setup.freqs), m_K(setup.K), m_Z_min(setup.Z_min), m_Z_max(setup.Z_max)
	{
		if (m_Zl.size() < m_freqs.size() && m_Zl.size() != 1)
			throw(std::invalid_argument("Number of load impedances & frequency points mismatch"));

		if (m_Zl.size() == 1)
		{
			m_Zl.resize(m_freqs.size());
			for (int i = 1; i < m_Zl.size(); i++)
				m_Zl[i] = m_Zl[0];
		}

		if (m_Z0 < 1e-6 || m_er < 1e-6 || m_d < 1e-6)
			throw(std::invalid_argument("Invalid NTL physical charachteristics"));
		if (m_Zs < 1e-6)
			throw(std::invalid_argument("Invalid source impedance"));
		for (auto& z : m_Zl)
			if (z < 1e-6)
				throw(std::invalid_argument("Invalid load impedance(s)"));

		if (m_freqs.empty())
			throw(std::invalid_argument("Empty frequency vector"));
		if (m_Z_min < 1e-6 || m_Z_max < 1e-6 || m_Z_max < m_Z_min)
			throw(std::invalid_argument("Invalid min/max impedance(s)"));
	}

	NTL_opt_result NTL_opt::optimise(console mode)
	{
		bool out = (mode == console::active) ? true : false;
		NTL ntl(m_Z0, m_er, m_d);
		auto start_time = std::chrono::high_resolution_clock::now();
		opt_result opt_result = optimiser(mode);
		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end_time - start_time;
		ntl.set_Cn(opt_result.optimised_cn);

		if (out)
			std::cout << "*** Optimisation finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		return NTL_opt_result(opt_result, ntl);
	}

	NTL_opt_result NTL_opt::optimise_d(double resolution, console mode)
	{
		bool out = (mode == console::active) ? true : false;

		NTL ntl(m_Z0, m_er, m_d);

		double init_d = m_d;
		int init_max_attempts = m_max_attempts;

		auto start_time = std::chrono::high_resolution_clock::now();
		if (out)
			std::cout << "\n\t===>>> Starting with d = " << m_d * 1000 << "mm" << std::endl;
		opt_result init_result = optimiser(mode);
		opt_result new_result = init_result;

		if (init_result.final_error < m_accepted_error)
		{
			double error_this_attempt = init_result.final_error;

			while (true)
			{
				m_d -= resolution;

				if (out)
					std::cout << "\n\t===>>> Trimming NTL by " << resolution * 1000 << "mm to:\t"
					<< m_d * 1000 << "mm" << std::endl;

				opt_result result_this_attempt = optimiser(mode);

				if (result_this_attempt.final_error > m_accepted_error)
				{
					m_d += resolution;
					break;
				}

				new_result = result_this_attempt;
			}

		}

		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end_time - start_time;

		ntl.set_d(m_d);
		ntl.set_Cn(new_result.optimised_cn);

		if (out)
		{
			std::cout << "Final NTL length:\t" << ntl.get_d() * 1000 << "mm" << std::endl;
			std::cout << "*** Optimisation finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		}

		m_d = init_d;
		m_max_attempts = init_max_attempts;

		return NTL_opt_result(new_result, ntl);
	}

	NTL_opt_result NTL_opt::optimise_d(console mode)
	{
		return optimise_d(1.0e-3, mode);
	}


	double NTL_opt::min_objective(const std::vector<double>& Cn) const
	{

		double sum_squares{};

		#pragma omp parallel
		{
			double thread_sum_squares{};
			#pragma omp for nowait
			for (int i = 0; i < m_freqs.size(); i++)
			{
				double f = m_freqs[i];
				double Zl = m_Zl[i];

				matrix2x2cd T = calculate_T_matrix(m_Z0, m_er, m_d, Cn, f, m_K);

				std::complex<double> a = T(0, 0), b = T(0, 1), c = T(1, 0), d = T(1, 1);
				std::complex<double> Zin = (Zl * a + b) / (Zl * c + d);
				std::complex<double> gamma = (Zin - m_Zs) / (Zin + m_Zs);
				thread_sum_squares += std::pow(std::norm(gamma), 2);
			}

			#pragma omp critical
			{
				sum_squares += thread_sum_squares;
			}
		}

		return std::sqrt(sum_squares);
	}

	void NTL_opt::equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		res[0] = calculate_Z(Cn, n, m_Z0, m_d, 0) - m_Z0;
		res[1] = calculate_Z(Cn, n, m_Z0, m_d, m_d) - m_Z0;
	}

	void NTL_opt::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = calculate_Z(Cn, n, m_Z0, m_d, (i + 0.5) * dz) - m_Z_max;
	}

	void NTL_opt::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = m_Z_min - calculate_Z(Cn, n, m_Z0, m_d, (i + 0.5) * dz);
	}

	double NTL_opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
	{
		double total_sum_squares{};
		int F = m_freqs.size();
		int N = Cn.size();

		#pragma omp parallel
		{
			double thread_sum_squares{};
			std::vector<double> thread_grad(N, 0);

			#pragma omp for nowait
			for (int i = 0; i < F; i++)
			{
				double f = m_freqs[i];
				double Zl = m_Zl[i];

				auto [T, dT] = calculate_T_matrix_with_grad(m_Z0, m_er, m_d, Cn, f, m_K);

				std::complex<double> A = T(0, 0);
				std::complex<double> B = T(0, 1);
				std::complex<double> C = T(1, 0);
				std::complex<double> D = T(1, 1);

				std::complex<double> num = A * Zl + B;
				std::complex<double> den = C * Zl + D;

				std::complex<double> Zin = num / den;
				std::complex<double> Gamma = (Zin - m_Zs) / (Zin + m_Zs);

				thread_sum_squares += std::pow(std::norm(Gamma), 2);
				std::complex<double> dGamma;

				for (int n = 0; n < N; n++)
				{
					std::complex<double> dA = dT[n](0, 0);
					std::complex<double> dB = dT[n](0, 1);
					std::complex<double> dC = dT[n](1, 0);
					std::complex<double> dD = dT[n](1, 1);

					dGamma += Zl * dA + dB - Zin * (Zl * dC + dD);
					dGamma *= 2 * m_Zs / (den * std::pow(Zin + m_Zs, 2));

					thread_grad[n] += 4.0 * std::norm(Gamma) * std::real(Gamma * std::conj(dGamma));

				}
			}

			#pragma omp critical
			{
				total_sum_squares += thread_sum_squares;
				if(!grad.empty())
					for (int n = 0; n < N; n++)
						grad[n] += thread_grad[n];
			}

		}

		return total_sum_squares;
	}

	//double NTL_opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
	//{
	//	// 1. Prepare global accumulators
	//	double total_sum_squares = 0.0;
	//	bool calc_grad = !grad.empty();

	//	if (calc_grad)
	//		std::fill(grad.begin(), grad.end(), 0.0);

	//	int num_freqs = m_freqs.size();
	//	int num_coeffs = Cn.size();

	//	// 2. Parallel Region
	//	#pragma omp parallel
	//	{
	//		// A. Thread-Local Accumulators (Stack allocated for speed)
	//		double thread_sum_squares = 0.0;
	//		std::vector<double> thread_grad;
	//		if (calc_grad)
	//			thread_grad.assign(num_coeffs, 0.0); // Initialize with zeros

	//		// B. Distribute Frequency Loop
	//		#pragma omp for nowait 
	//		for (int i = 0; i < num_freqs; i++)
	//		{
	//			double f = m_freqs[i];
	//			double Zl = m_Zl[i];

	//			// Call Model (Serial execution per thread)
	//			auto [T, T_grads] = GMN_calculate_T_matrix_with_grad(m_Z0, m_er, m_d, Cn, f, m_K);

	//			// --- Physics Calculation ---
	//			std::complex<double> A = T(0, 0), B = T(0, 1), C = T(1, 0), D = T(1, 1);
	//			std::complex<double> Zin_num = A * Zl + B;
	//			std::complex<double> Zin_den = C * Zl + D;
	//			std::complex<double> Zin = Zin_num / Zin_den;
	//			std::complex<double> Gamma = (Zin - m_Zs) / (Zin + m_Zs);

	//			// Accumulate Error (Local)
	//			thread_sum_squares += std::pow(std::norm(Gamma), 2);

	//			// Accumulate Gradient (Local)
	//			if (calc_grad)
	//			{
	//				for (int n = 0; n < num_coeffs; n++)
	//				{
	//					std::complex<double> dA = T_grads[n](0, 0);
	//					std::complex<double> dB = T_grads[n](0, 1);
	//					std::complex<double> dC = T_grads[n](1, 0);
	//					std::complex<double> dD = T_grads[n](1, 1);

	//					// Quotient Rule for Zin
	//					std::complex<double> dNum = dA * Zl + dB;
	//					std::complex<double> dDen = dC * Zl + dD;
	//					std::complex<double> dZin = (dNum * Zin_den - Zin_num * dDen) / (Zin_den * Zin_den);

	//					// Quotient Rule for Gamma
	//					std::complex<double> dGamma = (2.0 * m_Zs * dZin) / std::pow(Zin + m_Zs, 2);

	//					// Chain Rule for |G|^4
	//					thread_grad[n] += 4.0 * std::norm(Gamma) * std::real(Gamma * std::conj(dGamma));
	//				}
	//			}
	//		}

	//		// C. Critical Section: Merge Local Results to Global
	//		#pragma omp critical
	//		{
	//			total_sum_squares += thread_sum_squares;
	//			if (calc_grad)
	//			{
	//				for (int n = 0; n < num_coeffs; n++)
	//					grad[n] += thread_grad[n];
	//			}
	//		}
	//	}

	//	// 3. Final Scaling
	//	double objective_val = std::sqrt(total_sum_squares);

	//	if (calc_grad)
	//	{
	//		if (objective_val > 1e-12)
	//		{
	//			double scale = 0.5 / objective_val;
	//			for (auto& g : grad) g *= scale;
	//		}
	//		else
	//		{
	//			std::fill(grad.begin(), grad.end(), 0.0);
	//		}
	//	}

	//	return objective_val;
	//}

	//double NTL_opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
	//{

	//	const double f_base = min_objective(Cn);

	//	if (!grad.empty())
	//	{
	//		const double h = 1e-8;
	//		int N = Cn.size();

	//		#pragma omp parallel for
	//		for (int i = 0; i < N; ++i)
	//		{
	//			std::vector<double> Cn_nudged = Cn;
	//			Cn_nudged[i] += h;

	//			double f_nudged = min_objective(Cn_nudged);

	//			grad[i] = (f_nudged - f_base) / h;
	//		}
	//	}

	//	return f_base;
	//}
}