#pragma once
#include "qt_plot.h"

void plot1(qt_plot& plotter)
{
	plotter.combine({ plotter.plot([](double x)-> double {return x * x; }, 0, 10, 0, 100, 0.5),
		plotter.plot([](double x)-> double {return x * x * x; }, 0, 10, 0, 1000, 0.5) });
}

void plot2(qt_plot& plotter)
{
	plotter.plot([](double x)-> double {return x; }, 0, 10, 0, 10, 0.5);
}