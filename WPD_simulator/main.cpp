#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <complex>
#include <algorithm>
#include "qt_plot.h"
#include "models/ntl.h"
#include "models/wpd.h"
#include "simulation/NTL_sim.h"
#include "simulation/WPD_sim.h"
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

		NTL::NTL ntl2 = fh::file_to_ntl(folder_name + "/" + folder_name +"_ntl2");
		NTL::NTL ntl3 = fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl3");
		NTL::NTL out2 = fh::file_to_ntl(folder_name + "/" + folder_name + "_out2");
		NTL::NTL out3 = fh::file_to_ntl(folder_name + "/" + folder_name + "_out3");
		auto R = fh::file_to_json(folder_name + "/" + folder_name + "_wpd").at("R");
		auto setup = fh::file_to_setup<WPD::opt_setup>(folder_name + "/" + folder_name + "_setup");

		int K = 2 * setup.K;

		WPD::WPD wpd(ntl2, ntl3, R);

		double f_min = 1e7;
		double f_max = *std::max_element(setup.freqs.cbegin(), setup.freqs.cend());
		double f_step = 1e6;

		//simulate 	

		std::cout << "\n\tR = " << wpd.get_R() << std::endl;
		
		std::cout << "================================\n";
		std::cout << "     TRANSMISSION MATRICES\n";
		std::cout << "================================\n";

		std::cout << "\n\t-> Main Arms\n";

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

		std::cout << "\n\t-> Output Arms\n";

		for (auto& f : setup.freqs)
		{
			auto T2 = out2.T_matrix(f, K);
			auto T3 = out3.T_matrix(f, K);

			std::cout << "\n[Frequency: " << f << "]\n";
			std::cout << "\nT2:\n";
			std::cout << T2 << std::endl;
			std::cout << "\nT3:\n";
			std::cout << T3 << std::endl;
		}

		std::cout << "\n================================\n";
		std::cout << "         Odd Mode Gamma\n";
		std::cout << "================================\n";
		
		auto split = setup.split.cbegin();
		for (auto& f : setup.freqs)
		{
			auto T2 = ntl2.T_matrix(f, K);
			auto T3 = ntl3.T_matrix(f, K);
			
			double R2 = std::sqrt(*split) * setup.Zref;
			std::complex<double> ratio2 = T2(0, 1) / T2(0, 0);
			std::complex<double> Zeq2 = ratio2 * R2 / (ratio2 + R2);
			std::complex<double> GammaOdd2 = (Zeq2 - R2) / (Zeq2 + R2);

			double R3 = setup.Zref / std::sqrt(*split) ;
			std::complex<double> ratio3 = T3(0, 1) / T3(0, 0);
			std::complex<double> Zeq3 = ratio3 * R3 / (ratio3 + R3);
			std::complex<double> GammaOdd3 = (Zeq3 - R3) / (Zeq3 + R3);
			++split;

			std::cout << "\n[Frequency: " << f << "]\n";
			std::cout << "\nGamma2Odd:\n";
			std::cout << 20 * std::log10(std::max<double>(std::abs(GammaOdd2), 1e-12)) << std::endl;
			std::cout << "\nGamma3Odd:\n";
			std::cout << 20 * std::log10(std::max<double>(std::abs(GammaOdd3), 1e-12)) << std::endl;
			
		}

		split = setup.split.cbegin();

		std::cout << "\n================================\n";
		std::cout << "       Electrical Lengths\n";
		std::cout << "================================\n";

		std::cout << "\n\t-> Main Arms\n";

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

		std::cout << "\n\t-> Output Arms\n";

		for (auto& f : setup.freqs)
		{
			double theta2 = out2.electrical_length(f, K);
			double theta3 = out3.electrical_length(f, K);

			std::cout << "\n[Frequency: " << f << "]\n";
			std::cout << "\ntheta2:\n";
			std::cout << theta2 << std::endl;
			std::cout << "\ntheta3:\n";
			std::cout << theta3 << std::endl;
		}

		std::cout << "\n================================\n";
		std::cout << "         S-Parameters\n";
		std::cout << "================================\n\n";

		NTL::sim nsim;
		WPD::sim wsim(f_min, 1.2 * f_max, f_step, setup.freqs);

		wsim.add_window(nsim.w_h_profile(ntl2, "NTL2"));
		wsim.add_window(nsim.z_profile(ntl2, "NTL2"));

		wsim.add_window(nsim.w_h_profile(ntl3, "NTL3"));
		wsim.add_window(nsim.z_profile(ntl3, "NTL3"));

		wsim.add_window(nsim.w_h_profile(out2, "Out2"));
		wsim.add_window(nsim.z_profile(out2, "Out2"));

		wsim.add_window(nsim.w_h_profile(out3, "Out3"));
		wsim.add_window(nsim.z_profile(out3, "Out3"));

		wsim.sparams(wpd, setup.Zref, out2, out3, K);

		wsim.merge("WPD");		
		
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