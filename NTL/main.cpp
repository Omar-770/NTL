#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <omp.h>
#include <complex>
#include <algorithm>
#include "qt_plot.h"
#include "models/ntl.h"
#include "models/wpd.h"
#include "optimisation/NTL_opt.h"
#include "optimisation/WPD_opt.h"
#include "optimisation/WPD_opt_2.h"
#include "simulation/NTL_sim.h"
#include "simulation/WPD_sim.h"
#include "common/file_handler.h"
#include "common/enums.h"
#include "common/helpers.h"
#include "common/latex.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);

		
		//setup and run optimisation
		auto setup = fh::file_to_setup<WPD::opt_setup>("wpd_setup1");
		WPD::opt_2 opt(setup);
		WPD::opt_result result = opt.optimise();

		//build the result
		NTL::NTL out2 = result.output2;
		NTL::NTL out3 = result.output3;
		WPD::WPD wpd = result.wpd;

		//simulate the result
		double f_max = *std::max_element(setup.freqs.cbegin(), setup.freqs.cend());
		double f_min = 1e7;
		WPD::sim wsim(f_min, f_max * 1.2, 1e6, setup.freqs);
		NTL::sim nsim;

		wsim.w_h_profile(wpd);
		wsim.add_window(nsim.w_h_profile(out2, "Out2"));
		wsim.add_window(nsim.w_h_profile(out3, "Out3"));
		
		wsim.z_profile(wpd);
		wsim.add_window(nsim.z_profile(out2, "Out2-Z"));
		wsim.add_window(nsim.z_profile(out3, "Out3-Z"));
		
		wsim.sparams(wpd, setup.Zref, out2, out3);
		wsim.merge("WPD");

		app.exec();

		//Saving
		
		std::cout << "\n\n\nSave Results? [y/n]: ";
		std::string save;
		std::cin >> save;
		if(save == "y" || save == "Y")
		{	
			std::string folder;
			std::cout << "Folder Name: ";
			std::cin >> folder;
			//Latex

			std::cout << "Saving to CSVs...\n";

			if (NTL::latex::wpd(wpd, out2, out3, folder, folder) &&
				NTL::latex::sparams(wpd, out2, out3, setup.Zref, f_min, f_max, 1e6, setup.freqs, folder, folder) &&
				NTL::latex::table_data(wpd, out2, out3, folder, folder))
				std::cout << "CSVs save success\n";
			else
				std::cerr << "CSVs save failure\n";

			std::cout << "Saving objects...\n";
			double H = 4e-3;
			//json
			fh::wpd_to_file(wpd, folder + "/" + folder + "_wpd");
			fh::ntl_to_file(out2, folder + "/" + folder + "_out2");
			fh::ntl_to_file(out3, folder + "/" + folder + "_out3");
			fh::setup_to_file<WPD::opt_setup>(setup, folder + "/" + folder + "_setup");
			//AutoCAD Script
			fh::export_geometry_scr(out2, H, folder + "/" + folder + "_out2");
			fh::export_geometry_scr(out3, H, folder + "/" + folder + "_out3");
			fh::export_geometry_scr(wpd.get_ntl2(), H, folder + "/" + folder + "_ntl2");
			fh::export_geometry_scr(wpd.get_ntl3(), H, folder + "/" + folder + "_ntl3");

			std::cout << "\n\n\n=== Save Complete ===\n\n\n";
		}
		
		/*
		
		//Extract
		double Z0 = 50;
		double er = 4.6;
		double R = 115;

		std::vector<double> f{ 1.2e9, 2e9 };

		std::vector<double> ntl2_cn =
		{
			0.0945231, 0.314151, -0.453895, -0.362074, -0.306795, 0.745765,
			-0.211532, 0.108787, 0.336469, -0.138876, 0.012232, -0.0656455,
			-0.0731207, 0.000854035, -0.00123451, 0.00112109, -0.00114144, 0.000569087
				};

		std::vector<double> ntl3_cn =
		{
			-0.245902, 0.267936, 0.301251, 0.0834044, -0.480617, 0.63242,
			-0.147128, -0.212577, -0.272693, -0.368662, 0.175423, 0.123331,
			0.143876, -0.0398623, -0.0222308, -0.016033, 0.00392588, 0.0235757
		};
		
		std::vector<double> out2_cn =
		{
			-0.379811, -0.180068, 0.303148, 0.0589602, 0.0843374, 0.0678859,
			0.117573, -0.065974, 0.0772449, 0.00305118, -0.0737445, 0.0778793,
			-0.0904853, 0.0087686, -0.0547364, 0.044911, -0.0225169, 0.0112633
		};

		std::vector<double> out3_cn =
		{
			0.39695, 0.220506, -0.361802, -0.103416, -0.113679, -0.0305861,
			-0.0353034, -0.0513112, -0.00592409, 0.0162741, 0.0252745, 0.0435272,
			-0.000485493, 0.017384, -0.0338536, -0.190641, 0.0848196, 0.0569552
		};

		NTL::NTL ntl2(Z0, er, 76e-3, ntl2_cn, 5);
		NTL::NTL ntl3(Z0, er, 76e-3, ntl3_cn, 5);
		NTL::NTL out2(Z0, er, 32e-3, out2_cn, 5);
		NTL::NTL out3(Z0, er, 32e-3, out3_cn, 5);

		WPD::WPD wpd(ntl2, ntl3, R);		

		//Simulate
		
		NTL::sim nsim;
		WPD::sim wsim(0.9e9, 2.1e9, 1e6, f);

		wsim.w_h_profile(wpd);
		wsim.add_window(nsim.w_h_profile(out2, "Out2"));
		wsim.add_window(nsim.w_h_profile(out3, "Out3"));
		wsim.z_profile(wpd);
		wsim.add_window(nsim.z_profile(out2, "Out2"));
		wsim.add_window(nsim.z_profile(out3, "Out3"));

		wsim.sparams(wpd, Z0, out2, out3);
		
		wsim.merge("WPD");

		
		*/
	}
	catch (const std::exception& e)
	{
		std::cerr << "\nCRITICAL ERROR: " << e.what() << std::endl;
		std::cerr << "Press Enter to exit..." << std::endl;
		std::cin.get();
	}

	return 0;
}