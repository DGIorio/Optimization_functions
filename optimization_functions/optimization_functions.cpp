// optmization_functions.cpp : Define as funções exportadas para o aplicativo DLL.
//

#include "stdafx.h"

extern "C" __declspec(dllexport) double rosen(unsigned int dim, double* x);

extern "C" __declspec(dllexport) double ackley(unsigned int dim, double* x);

extern "C" __declspec(dllexport) double rastrigin(unsigned int dim, double* x);

extern "C" __declspec(dllexport) double griewank(unsigned int dim, double* x);

extern "C" __declspec(dllexport) double schwefel(unsigned int dim, double* x);

int main(int argc, char** argv);

extern "C"
{
	__declspec(dllexport) double rosen(unsigned int dim, double* x)
	{
		unsigned int i = 0u;
		double ret = 0.0;
		for (i = 0u; i < dim - 1u; i++)
		{
			ret += 100. * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) + (x[i] - 1) * (x[i] - 1);
		}
		return ret;
	}
}

extern "C"
{
	__declspec(dllexport) double ackley(unsigned int dim, double* x)
	{
		unsigned int i = 0u;
		double ret = 0.0;
		double omega = 2.0 * 3.141592653589793238463;
		double s1 = 0.0;
		double s2 = 0.0;
		double nepero = std::exp(1.0);

		for (i = 0u; i < dim; i++)
		{
			s1 += x[i] * x[i];
			s2 += std::cos(omega * x[i]);
		}
		ret = -20. * std::exp(-0.2 * std::sqrt(1.0 / static_cast<double>(dim) * s1))
			- std::exp(1.0 / static_cast<double>(dim) * s2) + 20. + nepero;
		return ret;
	}
}

extern "C"
{
	__declspec(dllexport) double rastrigin(unsigned int dim, double* x)
	{
		unsigned int i = 0u;
		double ret = 0.0;
		const auto omega = 2.0 * 3.141592653589793238463;
		for (i = 0u; i < dim; ++i)
		{
			ret += x[i] * x[i] - 10. * std::cos(omega * x[i]);
		}
		ret += 10. * static_cast<double>(dim);
		return ret;
	}
}

extern "C"
{
	__declspec(dllexport) double griewank(unsigned int dim, double* x)
	{
		unsigned int i = 0u;
		double ret = 0.0;
		double fr = 4000.;
		double retval = 0.0;
		double p = 1.0;

		for (i = 0u; i < dim; i++) {
			retval += x[i] * x[i];
		}
		for (i = 0u; i < dim; i++) {
			p *= std::cos(x[i] / std::sqrt(static_cast<double>(i) + 1.0));
		}
		ret = (retval / fr - p + 1.);
		return ret;
	}
}

extern "C"
{
	__declspec(dllexport) double schwefel(unsigned int dim, double* x)
	{
		unsigned int i = 0u;
		double ret = 0.0;
		for (i = 0u; i < dim; i++) {
			ret += x[i] * std::sin(std::sqrt(std::abs(x[i])));
		}
		ret = 418.9828872724338 * static_cast<double>(dim) - ret;
		return ret;
	}
}

int main(int argc, char** argv)
{
	/*
		argc is the argument count
		argv contains the arguments

		argv[0] = path
		argv[1] = function
		argv[2] = dim
		argv[3] = x1
		argv[4] = x2
		argv[N] = xN
	*/

	if (argc < 4)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << "argv[0] = " << " function" << std::endl;
		std::cerr << "argv[1] = " << " dimension" << std::endl;
		std::cerr << "argv[2] = " << " x1" << std::endl;
		std::cerr << "argv[2] = " << " x1" << std::endl;
		std::cerr << "argv[N] = " << " xN" << std::endl;
		return -1;
	}

	int i;
	//int dim;
	unsigned int dim;
	double* x;
	double retval;

	//dim = atoi(argv[2]);
	dim = (unsigned int)strtoul(argv[2], (char **)NULL, 10);
	x = new double[dim]();

	for (i = 3; i < argc; ++i)
	{
		x[i - 3] = atof(argv[i]);
	}

	if (strcmp(argv[1], "rosen") == 0)
	{
		retval = rosen(dim, x);
	}

	if (strcmp(argv[1], "ackley") == 0)
	{
		retval = ackley(dim, x);
	}

	if (strcmp(argv[1], "rastrigin") == 0)
	{
		retval = rastrigin(dim, x);
	}
	
	if (strcmp(argv[1], "griewank") == 0)
	{
		retval = griewank(dim, x);
	}

	if (strcmp(argv[1], "schwefel") == 0)
	{
		retval = schwefel(dim, x);
	}

	std::cout << std::setw(24) << std::setprecision(16) << retval << std::endl;
}