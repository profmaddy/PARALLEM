/*
 * bypass.h
 *
 *  Created on: 27 Feb 2020
 *      Author: ndm12
 */

#ifndef BYPASS_H_
#define BYPASS_H_

#include "Data.h"
#include "production.h"
#include "lem.h"

void transferFAtoFAgrid(Data* data);
void bypassrouting(Data* data, Data* device, int iter);

#endif


