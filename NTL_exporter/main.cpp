#include <iostream>
#include "export_WPD.h"

int main()
{
	std::string folder_name{};
	int choice{};

	while (true)
	{
		std::cout << "Folder Name: ";
		std::cin >> folder_name;

		std::cout << "\nSelect a device\n";
		std::cout << "1) Back\n2) WPD\n\n";
		std::cout << "=> ";
		std::cin >> choice;
		if (std::cin.fail())
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Invalid input. Please enter numbers only.\n";
			continue;
		}

		if (choice == 1)
		{
			std::cout << "\n==========\n\n";
			continue;
		}

		double substrate_height{};
		double step{};
		switch (choice)
		{
		case 2:		

			std::cout << "Substrate Height (mm): ";
			std::cin >> substrate_height;
			if (std::cin.fail())
			{
				std::cin.clear();
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				std::cout << "Invalid input. Please enter numbers only.\n";
				continue;
			}
			std::cout << "Step Size (mm): ";
			std::cin >> step;
			if (std::cin.fail())
			{
				std::cin.clear();
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				std::cout << "Invalid input. Please enter numbers only.\n";
				continue;
			}
			substrate_height /= 1000.0;
			step /= 1000.0;
			try
			{
				export_WPD_AutoCAD(folder_name, substrate_height, step);
			}
			catch (const std::exception& e)
			{
				std::cerr << e.what() << std::endl;
				continue;
			}
			break;

		default:
			std::cout << "Invalid Input" << std::endl;
			continue;
		}
	}

	return 0;
}