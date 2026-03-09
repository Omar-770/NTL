#pragma once

#include <string>
#include <vector>
#include <utility>
#include <math.h>
#include "models/ntl.h"


std::vector<std::pair<double, double>> get_w_vector(const NTL::NTL& ntl,
	double start_pos, double end_pos, double substrate_height, double step = 1e-4);

void seg(const NTL::NTL& ntl, double start_pos, double end_pos,
	double substrate_height, const std::string& filename, double step = 1e-4);

void arc(const NTL::NTL& ntl, double start_pos, double end_pos, double angle,
	double substrate_height, const std::string& filename, double step = 1e-4);