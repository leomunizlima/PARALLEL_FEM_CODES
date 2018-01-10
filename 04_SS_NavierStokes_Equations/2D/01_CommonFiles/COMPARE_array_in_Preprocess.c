#include "SSNavierStokesEquations.h"

int COMPARE_array_in_Preprocess (const void * a, const void * b)
{
	if (((ARRAYType*)a)->array1 <  ((ARRAYType*)b)->array1) return -1;
	if (((ARRAYType*)a)->array1 >  ((ARRAYType*)b)->array1) return  1;
	
	return 0;
}


