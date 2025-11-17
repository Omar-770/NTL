#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "qt_plot.h"
#include "models/ntl.h"
#include "optimisation/ntl_opt.h"


int main(int argc, char* argv[])
{
	NTL::NTL_opt_setup opt_setup =
	{
		50,4.6,80e-3,50,100,{0.5e9, 1.3e9, 2e9},50,10,0.1 * 50,
		2.6 * 50,1e-6,1e-6,350e3,
		50e3,6e-3
	};

	NTL::NTL_opt opt(opt_setup);
	NTL::NTL ntl = opt.optimise();


	/*try
	{
		QApplication app(argc, argv);


		double Z0 = 50;
		double Zl = 100;
		double Zs = 50;
		std::vector<double> f = { 0.5e9, 1.3e9, 2e9 };
		double e_r = 4.6;
		int K = 50;
		double d = 80e-3;
		int N = 10;

		NTL::NTL ntl(Z0, e_r, d);
		std::vector<double> Cn = { 0.27029971731871727, -0.1406883396667979, -0.19437319260269542,
			-0.6318277441297461, 0.2638584614123822, 0.11608394752540892, 0.06819894141855014,
			0.13368659184146162, 0.07296676694708622, 0.041794852549067875 };


		ntl.set_Cn(Cn);

		double f_max = *max_element(f.begin(), f.end());
		double f_min = 1;
		double f_step = 1e7;
		double points = 1.1 * f_max / f_step + 1;

		std::vector<std::pair<double, double>> Z_vec = ntl.get_Z_vec(0.1e-3);
		std::vector<std::pair<double, double>> w_h_vec = ntl.get_w_h_vec(0.1e-3);
		std::vector<std::pair<double, double>> A, B, C, D;
		std::vector<std::pair<double, double>> S11, S12, S21, S22;
		S11.reserve(points); S12.reserve(points); S21.reserve(points); S22.reserve(points);
		A.reserve(points); B.reserve(points); C.reserve(points); D.reserve(points);


		for (double freq = f_min; freq <= 1.1 * f_max; freq += f_step)
		{
			auto T = ntl.T_matrix(freq, K);
			auto S = ntl.S_matrix(freq, Zs, Zl, K);

			A.push_back({ freq, std::abs((T(0, 0))) });
			B.push_back({ freq, std::abs((T(0, 1))) });
			C.push_back({ freq, std::abs((T(1, 0))) });
			D.push_back({ freq, std::abs((T(1, 1))) });

			S11.push_back({ freq, 20 * std::log10(std::abs((S(0, 0)))) });
			S12.push_back({ freq, 20 * std::log10(std::abs((S(0, 1)))) });
			S21.push_back({ freq, 20 * std::log10(std::abs((S(1, 0)))) });
			S22.push_back({ freq, 20 * std::log10(std::abs((S(1, 1)))) });
		}

		qt_plot qt;

		qt.combine({ qt.plot(Z_vec, "Impedance Profile Z(z)"),
			qt.plot_mirror(w_h_vec, "W_H(z)"),
			qt.plot(S11 , "S11"), qt.plot(S12 , "S12"),
			qt.plot(S21 , "S21"), qt.plot(S22 , "S22"),
			qt.plot(A , "A"), qt.plot(B , "B"),
			qt.plot(C , "C"), qt.plot(D , "D")
			});

		app.exec();
	}
	catch (const std::runtime_error& e)
	{
		std::cerr << "NTL_sim, main.cpp: ";
		std::cerr << e.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "NTL_sim, main.cpp: ";
		std::cerr << "Ran into an unknown exception..." << std::endl;
	}*/

	return 0;

}