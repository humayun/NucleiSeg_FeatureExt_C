/*
 * RoutineFunctions.h
 *
 *  Created on: 20 May 2012
 *      Author: humayun
 */
#ifndef __ROUTINEFUNCTIONS_H__
#define __ROUTINEFUNCTIONS_H__

#include "ITKDeclarations.h"

CharImagePointer SegmentationPipeLine( CharImagePointer inImage , CharImageIndexType seedPoint );

CharImagePointer SegmentationPipeLine( RGBImagePointer inImage, CharImageIndexType seedPoint );

void ComputeFeatures( CharImagePointer inImage, CharImagePointer maskImage, std::vector<double> &patchFeatures );

#endif // __ROUTINEFUNCTIONS_H__
