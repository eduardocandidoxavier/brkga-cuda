#ifndef DECODER_H
#define DECODER_H

#include <stdio.h>
#include <iostream>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/device_ptr.h>
#include <algorithm>

#include "BRKGA.h"
#include "CommonStructs.h"
#include "SetCoveringDecoder.h"


__device__ float device_decode(float *chromosome, int n, void *d_instance_info);
float host_decode(float *chromosome, int n, void *d_instance_info);
__device__ float device_decode_chromosome_sorted(ChromosomeGeneIdxPair *chromosome, int n, void *d_instance_info);

#endif
