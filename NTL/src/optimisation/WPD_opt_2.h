#pragma once
#include "optimisation/optimiser.h"
#include "optimisation/WPD_opt.h"
#include <memory>

namespace WPD
{
	class opt_2;

	class opt_2 : public optimiser
	{
	public:
		opt_2(const opt_setup& setup);

		opt_result optimise(console mode = console::active);

	private:
		double m_Z0;
		double m_er;
		double m_d;
		int m_M;
		double m_Zref; //port1 impedance and load impedance for output transformers
		std::vector<double> m_freqs;
		int m_K;
		double m_Z_min;
		double m_Z_max;
		double m_R_min;
		double m_R_max;
		std::vector<double> m_split; //P3/P2

		int m_F;

	private:
		NTL::opt_setup m_ntl_setup;
		std::array<std::unique_ptr<NTL::opt>, 4> m_ntl_opt;

	private:
		NTL::NTL m_out2, m_out3;
		NTL::NTL m_ntl2, m_ntl3;
	};
}