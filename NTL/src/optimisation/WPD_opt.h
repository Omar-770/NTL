#pragma once

#include "optimiser.h"
#include <omp.h>
#include <iostream>
#include <nlohmann/json.hpp>
#include <chrono>

namespace NTL
{
	class WPD_opt;
	struct WPD_opt_setup;
	struct WPD_opt_result;
}