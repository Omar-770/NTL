#include <iostream>
#include "models/ntl.h"
#include "common/file_handler.h"
#include "common/helpers.h"
#include "NTL_bender.h"
#include "NTL_freestyle.h"

namespace fh = NTL::fh;

int main()
{
	auto check_input = [](std::istream& cin)
	{
		if (cin.fail())
		{
			cin.clear();
			cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Invalid input. Please enter numbers only.\n";
		}
	};

	int choice{};

	while (true)
	{		
		std::cout << "Select Mode: ";
		std::cout << "1) Free Mode\n2) WPD Assembly\n\n";
		std::cout << "=> ";
		std::cin >> choice;
		check_input(std::cin);

		switch (choice)
		{
		case 1:
			NTL_free_style();
			break;
		case 2:
			//assemble_WPD();
			break;
		default:
			break;
		}

	}
	

	return 0;
}