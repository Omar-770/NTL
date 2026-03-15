#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <utility>
#include <math.h>
#include "models/ntl.h"


std::vector<std::pair<double, double>> inline get_w_vector(const NTL::NTL& ntl,
	double start_pos, double end_pos, double substrate_height, double step = 1e-4)
{
	int n_steps = std::round((end_pos - start_pos) / step);

	std::vector<std::pair<double, double>> w_h;
	w_h.reserve(n_steps + 1);

	for (int i = 0; i < n_steps + 1; i++)
	{
		double pos = start_pos + double(i) * step;
		if (pos >= end_pos)
		{
			pos = end_pos;
		}

		w_h.emplace_back(pos, substrate_height * ntl.W_H(pos));
	}

	return w_h;
}

void inline seg(const NTL::NTL& ntl, double start_pos, double end_pos,
	double substrate_height, const std::string& filename, double step = 1e-4)
{
	if (start_pos < 0 || start_pos >= end_pos || end_pos >(ntl.get_d() + 1e-6))
		throw std::runtime_error("Invalid start or end position");

	std::fstream stream;
	std::stringstream ss;
	ss << std::fixed << std::setprecision(1);
	ss << filename << "_seg_" << (start_pos * 1000.0) << "_" << (end_pos * 1000.0) << ".scr";
	std::string name = ss.str();
	stream.open(name, std::ios_base::out);

	if (!stream.is_open())
		throw std::runtime_error("Could not open file for AutoCAD export");

	stream << "_OSMODE 0\n";
	stream << "_UCS _World\n";
	stream << "_PLINE\n";

	auto w = get_w_vector(ntl, start_pos, end_pos, substrate_height, step);

	for (int i = 0; i < w.size(); i++)
	{
		stream << 1000 * w[i].first << "," << 0.5 * 1000 * w[i].second << "\n";
	}

	for (int i = w.size() - 1; i >= 0; i--)
	{
		stream << 1000 * w[i].first << "," << -0.5 * 1000 * w[i].second << "\n";
	}

	stream << "_CLOSE\n";
	stream << "_REGION\n_L\n\n";
	stream << "_ZOOM _E\n";

	stream.close();
	std::cout << "AutoCAD script exported to " << name << std::endl;
	std::cout << std::endl;
}


void inline arc(const NTL::NTL& ntl, double start_pos, double end_pos, double angle,
	double substrate_height, const std::string& filename, double step = 1e-4)
{
	if (start_pos < 0 || start_pos >= end_pos || end_pos >(ntl.get_d() + 1e-6))
		throw std::runtime_error("Invalid start or end position");

	if (angle <= 0 || angle > 360)
		throw std::invalid_argument("Angle must be between 0 and 360");

	std::fstream stream;
	std::stringstream ss;
	ss << std::fixed << std::setprecision(1);
	ss << filename << "_arc_" << (start_pos * 1000.0) << "_" << (end_pos * 1000.0)
		<< "_ang_" << angle << ".scr";
	std::string name = ss.str();
	stream.open(name, std::ios_base::out);

	if (!stream.is_open())
		throw std::runtime_error("Could not open file for AutoCAD export");

	stream << "_OSMODE 0\n";
	stream << "_UCS _World\n";
	stream << "_PLINE\n";

	auto w = get_w_vector(ntl, start_pos, end_pos, substrate_height, step);

	double Theta = angle * NTL::M_PI / 180.0;
	double R = (end_pos - start_pos) / Theta;

	double step_theta = Theta / (w.size() - 1);

	for (int i = 0; i < w.size(); i++)
	{
		double theta = double(i) * step_theta;
		if (theta >= Theta)
		{
			theta = Theta;
		}

		double x = 1000 * (R - 0.5 * w[i].second) * std::cos(theta);
		double y = 1000 * (R - 0.5 * w[i].second) * std::sin(theta);

		stream << x << "," << y << "\n";
	}

	for (int i = w.size() - 1; i >= 0; i--)
	{
		double theta = double(i) * step_theta;
		if (theta <= 0)
		{
			theta = 0;
		}

		double x = 1000 * (R + 0.5 * w[i].second) * std::cos(theta);
		double y = 1000 * (R + 0.5 * w[i].second) * std::sin(theta);

		stream << x << "," << y << "\n";
	}

	stream << "_CLOSE\n";
	stream << "_REGION\n_L\n\n";
	stream << "_ZOOM _E\n";

	stream.close();
	std::cout << "AutoCAD script exported to " << name << std::endl;
	std::cout << std::endl;
}