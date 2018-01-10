#include "obliquo.h"

double OBLIQUO_v2presc(double x, double y){
	double v2;
	
	//v2 = -sin 10°
	if (fabs(y)<=1e-15 && fabs(x)>0){
		v2 = 0.0;
	}else{
		v2 = -0.17365;//-0.173648177;
	}
	
	return v2;
} 


