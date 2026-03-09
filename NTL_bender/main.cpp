#include <iostream>
#include "models/ntl.h"
#include "common/file_handler.h"
#include "common/helpers.h"
#include "bender.h"

namespace fh = NTL::fh;

int main()
{

	std::string file_name{};
	double start_pos{};
	double end_pos{};
	double angle{};
	double substrate_height{ 1.6e-3 };
	double step{ 1e-4 };
	int choice{};
	NTL::NTL ntl;

	while (true)
	{
		try
		{
			std::cout << "File Name: ";
			std::cin >> file_name;
			ntl = fh::file_to_ntl(file_name);
		}
		catch (const std::exception& e)
		{
			std::cerr << e.what() << std::endl;
			continue;
		}

		while (true)
		{
			std::cout << "\nSelect an option\n";
			std::cout << "1) Segment\n2) Arc\n3) Back\n\n";
			std::cout << "=> ";
			std::cin >> choice;
			if (std::cin.fail()) 
			{
				std::cin.clear(); 
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				std::cout << "Invalid input. Please enter numbers only.\n";
				continue;
			}

			if (choice == 3)
			{
				std::cout << "\n==========\n\n";
				break;
			}

			switch (choice)
			{
			case 1:
				std::cout << "Start Position (mm): ";
				std::cin >> start_pos;
				if (std::cin.fail())
				{
					std::cin.clear();
					std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					std::cout << "Invalid input. Please enter numbers only.\n";
					continue;
				}
				std::cout << "End Position (mm): ";
				std::cin >> end_pos;
				if (std::cin.fail())
				{
					std::cin.clear();
					std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					std::cout << "Invalid input. Please enter numbers only.\n";
					continue;
				}
				try
				{
					seg(ntl, start_pos / 1000.0, end_pos / 1000.0, substrate_height, file_name, step);
				}
				catch (const std::exception& e)
				{
					std::cerr << e.what() << std::endl;
					continue;
				}
				catch (...)
				{
					std::cerr << "\n\nERROR: UNKOWN EXCEPTION!\n\n";
					continue;
				}
				break;

			case 2:
				std::cout << "Start Position (mm): ";
				std::cin >> start_pos;
				if (std::cin.fail())
				{
					std::cin.clear();
					std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					std::cout << "Invalid input. Please enter numbers only.\n";
					continue;
				}
				std::cout << "End Position (mm): ";
				std::cin >> end_pos;
				if (std::cin.fail())
				{
					std::cin.clear();
					std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					std::cout << "Invalid input. Please enter numbers only.\n";
					continue;
				}
				std::cout << "Angle(deg): ";
				std::cin >> angle;
				if (std::cin.fail())
				{
					std::cin.clear();
					std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					std::cout << "Invalid input. Please enter numbers only.\n";
					continue;
				}
				try
				{
					arc(ntl, start_pos / 1000.0, end_pos / 1000.0, angle, substrate_height, file_name, step);
				}
				catch (const std::exception& e)
				{
					std::cerr << e.what() << std::endl;
					continue;
				}
				catch (...)
				{
					std::cerr << "\n\nERROR: UNKOWN EXCEPTION!\n\n";
					continue;
				}
				break;

			default:
				std::cout << "Invalid Input" << std::endl;
				continue;
			}
		}
	}

	return 0;
}