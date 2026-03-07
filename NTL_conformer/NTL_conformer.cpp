#include <iostream>
#include "models/ntl.h"
#include "models/wpd.h"
#include "common/file_handler.h"
#include "common/helpers.h"


namespace fh = NTL::fh;

int main()
{
	while (true)
	{
		std::string file_name;
		double angle{};
		std::cout << "Enter file name: ";
		std::cin >> file_name;
		std::string name = "../../dir/" + file_name;
		NTL::NTL ntl;

		try
		{
			ntl = fh::file_to_ntl(name);
		}
		catch (const std::exception& e)
		{
			std::cerr << e.what() << "\n\n";
			continue;
		}

		std::cout << "Enter rotation angle: ";
		std::cin >> angle;

		fh::bend_geometry_scr(ntl, 1.6e-3, name, angle);
	}
}