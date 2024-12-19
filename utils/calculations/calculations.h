#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "inmost.h"
using namespace INMOST;


void gradient(Cell c, BlockEntry& UVW, vMatrix& G);

void gradient_to_strain(const vMatrix& G, vMatrix& strain);

void stress_to_traction(const vMatrix& stress, const rMatrix& n, vMatrix& trac);


#endif //CALCULATIONS_H
