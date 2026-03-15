#include "models/ntl.h"
#include "models/wpd.h"
#include "common/file_handler.h"
#include <iostream>
#include <sstream>
#include <string>

#include "AutoCAD.h"

inline void export_WPD_AutoCAD(const std::string& folder_name, double substrate_height, double step = 1e-4) 
{
	std::cout << "Checking validity...\n";

	NTL::NTL ntl2 = NTL::fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl2");
	NTL::NTL ntl3 = NTL::fh::file_to_ntl(folder_name + "/" + folder_name + "_ntl3");
	NTL::NTL out2 = NTL::fh::file_to_ntl(folder_name + "/" + folder_name + "_out2");
	NTL::NTL out3 = NTL::fh::file_to_ntl(folder_name + "/" + folder_name + "_out3");

	scr(ntl2, substrate_height, folder_name + "/" + folder_name + "_ntl2", step);
	scr(ntl3, substrate_height, folder_name + "/" + folder_name + "_ntl3", step);
	scr(out2, substrate_height, folder_name + "/" + folder_name + "_out2", step);
	scr(out3, substrate_height, folder_name + "/" + folder_name + "_out3", step);
}