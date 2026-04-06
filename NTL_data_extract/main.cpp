#include <iostream>
#include "latex.h"
#include "common/file_handler.h"
#include "optimisation/wpd_opt.h"
#include "optimisation/ntl_opt.h"
#include "optimisation/tj_opt.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		std::string folder_name;
		std::cout << "Folder Name: ";
		std::cin >> folder_name;
		std::cout << "\n\n\n";

		//build

		NTL::NTL ntl2 = fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl2");
		NTL::NTL ntl3 = fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl3");
		NTL::NTL out2 = fh::file_to_ntl(folder_name + "/" + folder_name + "_out2");
		NTL::NTL out3 = fh::file_to_ntl(folder_name + "/" + folder_name + "_out3");
		auto R = fh::file_to_json(folder_name + "/" + folder_name + "_wpd").at("R");
		auto setup = fh::file_to_setup<WPD::opt_setup>(folder_name + "/" + folder_name + "_setup");
	

		WPD::WPD wpd(ntl2, ntl3, R);

		//NTL::NTL ntl2 = fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl2");
		//NTL::NTL ntl3 = fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl3");
		//auto setup = fh::file_to_setup<TJ::opt_setup>(folder_name + "/" + folder_name + "_setup");

		//TJ::TJ tj(ntl2, ntl3);

		//extract
		
		if (wpd_profiles(wpd, out2, out3, folder_name, folder_name))
			std::cout << "Extracted profiles...\n";
			
		if (wpd_table_data(wpd, out2, out3, folder_name, folder_name))
			std::cout << "Extracted table data...\n";

		if (wpd_sparams(wpd, out2, out3, setup.Zref, 0.1e9, 3e9, 1e6, setup.freqs, folder_name, folder_name))
			std::cout << "Extracted sparams...\n";
		 
		//if (tj_profiles(tj, folder_name, folder_name))
		//	std::cout << "Extracted profiles...\n";
		//	
		//if (tj_table_data(tj, folder_name, folder_name))
		//	std::cout << "Extracted table data...\n";

		//if (tj_sparams(tj, setup.Zref, 0.1e9, 3e9, 1e6, setup.freqs, folder_name, folder_name))
		//	std::cout << "Extracted sparams...\n";
		
	}
	catch (const std::exception& e)
	{
		std::cerr << "\nCRITICAL ERROR: " << e.what() << std::endl;
		std::cerr << "Press Enter to exit..." << std::endl;
		std::cin.get();
	}

	return 0;
}