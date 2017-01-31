/*
* UtilityFunctions.cpp
*
*  Created on: 20 May 2014
*      Author: Humayun
*/

#include "UtilityFunctions.h"

double string_to_double(const string& s)
{
	istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

double ComputePDFByFormula(double mean, double variance, double value)
{
	double e = exp(-(value - mean) * (value - mean) / 2 * variance);
	double sD = sqrt(variance * 2 * 3.14);
	double pDF = e / sD;
	return pDF;
}

void ComputeVectorMeanVariance(vector<double> data, double &mean, double &variance)
{
	double sum = 0.0;
	for (int i = 0; i<data.size(); i++)
		sum += data[i];
	mean = sum / data.size();

	double temp = 0;
	for (int i = 0; i<data.size(); i++)
		temp += (mean - data[i])*(mean - data[i]);
	variance = temp / data.size();
	double standardDeviation = sqrt(variance);
}

