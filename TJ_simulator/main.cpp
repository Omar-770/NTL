#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <complex>
#include <algorithm>
#include "qt_plot.h"
#include "models/ntl.h"
#include "models/tj.h"
#include "simulation/NTL_sim.h"
#include "simulation/TJ_sim.h"
#include "optimisation/TJ_opt.h"
#include "common/file_handler.h"
#include "common/enums.h"
#include "common/helpers.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);

		std::string folder_name;
		std::cout << "Folder Name: ";
		std::cin >> folder_name;
		std::cout << "\n\n\n";

		//build

		NTL::NTL ntl2 = fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl2");
		NTL::NTL ntl3 = fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl3");

		auto setup = fh::file_to_setup<TJ::opt_setup>(folder_name + "/" + folder_name + "_setup");
		int K = 2 * setup.K;

		TJ::TJ tj(ntl2, ntl3);

		double f_min = 1e7;
		double f_max = *std::max_element(setup.freqs.cbegin(), setup.freqs.cend());
		double f_step = 1e6;

		//simulate 	

		std::cout << "================================\n";
		std::cout << "     TRANSMISSION MATRICES\n";
		std::cout << "================================\n";
		

		for (auto& f : setup.freqs)
		{
			auto T2 = ntl2.T_matrix(f, K);
			auto T3 = ntl3.T_matrix(f, K);

			std::cout << "\n[Frequency: " << f << "]\n";
			std::cout << "\nT2:\n";
			std::cout << T2 << std::endl;
			std::cout << "\nT3:\n";
			std::cout << T3 << std::endl;
		}

		std::cout << "\n================================\n";
		std::cout << "       Electrical Lengths\n";
		std::cout << "================================\n";

		for (auto& f : setup.freqs)
		{
			double theta2 = ntl2.electrical_length(f, K);
			double theta3 = ntl3.electrical_length(f, K);

			std::cout << "\n[Frequency: " << f << "]\n";
			std::cout << "\ntheta2:\n";
			std::cout << theta2 << std::endl;
			std::cout << "\ntheta3:\n";
			std::cout << theta3 << std::endl;
		}	

		std::cout << "\n================================\n";
		std::cout << "   Input Impedance/Admittance\n";
		std::cout << "================================\n";

		for (auto& f : setup.freqs)
		{
			std::complex<double> zin2 = ntl2.Zin(setup.Zref, f, K);
			std::complex<double> zin3 = ntl3.Zin(setup.Zref, f, K);
			
			std::complex<double> yin2 = 1.0 / zin2;
			std::complex<double> yin3 = 1.0 / zin3;

			std::complex<double> yin = yin2 + yin3;
			std::complex<double> zin = 1.0 / yin;

			std::cout << "\n[Frequency: " << f << "]\n";
			std::cout << "\nZin2\tYin2:\n";
			std::cout << zin2 << '\t' << yin2 << std::endl;
			std::cout << "\nZin3\tYin3:\n";
			std::cout << zin3 << '\t' << yin3 << std::endl;
			std::cout << "\nZin\tYin:\n";
			std::cout << zin << '\t' << yin << std::endl;

		}		

		std::cout << "\n================================\n";
		std::cout << "         S-Parameters\n";
		std::cout << "================================\n\n";

		NTL::sim nsim(f_min, 1.2 * f_max, f_step);
		TJ::sim tsim(f_min, 1.2 * f_max, f_step, setup.freqs);


		tsim.add_window(nsim.w_h_profile(ntl2, "NTL2"));
		tsim.add_window(nsim.z_profile(ntl2, "NTL2"));
		tsim.add_window(nsim.zin(ntl2, setup.Zref).magnitude());
		tsim.add_window(nsim.zin(ntl2, setup.Zref).phase());

		tsim.add_window(nsim.w_h_profile(ntl3));
		tsim.add_window(nsim.z_profile(ntl3));
		tsim.add_window(nsim.zin(ntl3, setup.Zref).magnitude());
		tsim.add_window(nsim.zin(ntl3, setup.Zref).phase());


		std::array<std::complex<double>, 3> Zref;
		Zref[0] = { setup.Zref, 0 };
		Zref[1] = { setup.Zref, 0 };
		Zref[2] = { setup.Zref, 0 };

		tsim.sparams(tj, Zref, K);

		tsim.merge("TJ");

		app.exec();
	}
	catch (const std::exception& e)
	{
		std::cerr << "\nCRITICAL ERROR: " << e.what() << std::endl;
		std::cerr << "Press Enter to exit..." << std::endl;
		std::cin.get();
	}

	return 0;
}