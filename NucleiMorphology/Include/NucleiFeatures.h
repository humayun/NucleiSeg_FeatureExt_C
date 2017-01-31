/*
 * NucleiFeatures.h
 *
 *  Created on: 21 May 2014
 *      Author: Humayun
 */

#ifndef __NUCLEIFEATURES_H__
#define __NUCLEIFEATURES_H__

#include <algorithm> // std::min_element
#include <vector>

using namespace std;

class NucleiFeatures {

public:
	NucleiFeatures() { };
	~NucleiFeatures() { };

	void normaliseFeatures(
					vector< vector< double > >& inData2D,
					vector< vector< double > >& outData2D,
					int normaliseRange = 1000 );

	void normaliseFeatures(
					vector< vector< double > >& inData2D,
					vector< vector< long > >& outData2D,
					int normaliseRange = 1000 );

};

#endif /* __NUCLEIFEATURES_H_ */
