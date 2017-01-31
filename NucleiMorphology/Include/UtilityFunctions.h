/*
 * UtilityFunctions.h
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */
#ifndef __UTILITYFUNCTIONS_H__
#define __UTILITYFUNCTIONS_H__

#include <string>
#include <vector>
#include <sstream>
#include <math.h>

using namespace std;

double string_to_double(const string& s);
double ComputePDFByFormula(double mean, double variance, double value);
void ComputeVectorMeanVariance(vector<double> data, double &mean, double &variance);

#endif // __UTILITYFUNCTIONS_H__
