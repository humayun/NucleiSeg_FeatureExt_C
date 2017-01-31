/*
 * NucleiFeatures.cpp
 *
 *  Created on: 21 May 2014
 *      Author: Humayun
 */

#include "NucleiFeatures.h"

void NucleiFeatures::normaliseFeatures(
					vector< vector< double > >& inData2D,
					vector< vector< double > >& outData2D,
					int normaliseRange /* = 1000 */ ) {

	for( int r = 0; r < inData2D.size(); r++ ) {
		vector< double > tmp;
		for( int c = 0; c < inData2D[r].size(); c++)
			tmp.push_back( inData2D[r][c] );
		outData2D.push_back( tmp );
	}

	// Start from First features then N-1 Features because last column is class label
	for( int c = 0; c < inData2D[0].size()-1; c++) {
		vector< double > oneFeature;
		// Extract one feature values for all instances
		for( int r = 0; r < inData2D.size(); r++ )
			oneFeature.push_back( inData2D[r][c] );

		double min = *std::min_element(oneFeature.begin(), oneFeature.end());
		double max = *std::max_element(oneFeature.begin(), oneFeature.end());

		// Normalize one feature value for all instances
		for( int r = 0; r < inData2D.size(); r++ )
		if( 0.0 != (max - min) )
				outData2D[r][c] = normaliseRange * ( oneFeature[r] - min ) / (max - min);
			else
				outData2D[r][c] = 0;
	}
}

void NucleiFeatures::normaliseFeatures(
					vector< vector< double > >& inData2D,
					vector< vector< long > >& outData2D,
					int normaliseRange /* = 1000 */ ) {

	for( int r = 0; r < inData2D.size(); r++ ) {
		vector< long > tmp;
		for( int c = 0; c < inData2D[r].size(); c++)
			tmp.push_back( inData2D[r][c] );
		outData2D.push_back( tmp );
	}

	// Start from First features then N-1 Features because last column is class label
	for( int c = 0; c < inData2D[0].size()-1; c++) {
		vector< double > oneFeature;
		// Extract one feature values for all instances
		for( int r = 0; r < inData2D.size(); r++ )
			oneFeature.push_back( inData2D[r][c] );

		double min = *min_element(oneFeature.begin(), oneFeature.end());
		double max = *max_element(oneFeature.begin(), oneFeature.end());

		// Normalize one feature value for all instances
		for( int r = 0; r < inData2D.size(); r++ )
			if( 0.0 != (max - min) )
				outData2D[r][c] = normaliseRange * ( oneFeature[r] - min ) / (max - min);
			else
				outData2D[r][c] = 0;
	}
}
