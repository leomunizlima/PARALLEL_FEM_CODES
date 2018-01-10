#include "cylinder.h"

void CYLINDER_BC_no_penetrability(int J1, int J2, int J3, NodeType *Node, double alpha_dt,  
                                    double delta, double deltaNMV, double Area, double y23, double y12, 
				    double y31, double x32, double x21, double x13, 
				    double Ax[4][4], double Ay[4][4], double M1[12][12], double N1[12][4], double M2[4][12], 
				    double R1[12], double R2[4], double Ue[12], double dUe[12], double UeB[4], double dUeB[4])
{
	int tag = 0;
	double theta1, theta2, theta3, sin_theta1, sin_theta2, sin_theta3, cos_theta1, cos_theta2, cos_theta3;

	if ((Node[J1].Type[1] == -1 || Node[J1].Type[2] == -1)&&(fabs(Node[J1].x)>1e-5)){
		tag = 1;
		theta1 = CYLINDER_theta(Node[J1].x, Node[J1].y); 
		sin_theta1 = sin(theta1);
		cos_theta1 = cos(theta1);
	}
	else{
		sin_theta1 = 0.0;
		cos_theta1 = 1.0; 
	}

	if ((Node[J2].Type[1] == -1 || Node[J2].Type[2] == -1)&&(fabs(Node[J2].x)>1e-5)){
		tag = 1;
		theta2 = CYLINDER_theta(Node[J2].x, Node[J2].y); 
		sin_theta2 = sin(theta2);
		cos_theta2 = cos(theta2);
	}
	else{
		sin_theta2 = 0.0;
		cos_theta2 = 1.0; 
	}

	if ((Node[J3].Type[1] == -1 || Node[J3].Type[2] == -1)&&(fabs(Node[J3].x)>1e-5)){
		tag = 1;
		theta3 = CYLINDER_theta(Node[J3].x, Node[J3].y); 
		sin_theta3 = sin(theta3);
		cos_theta3 = cos(theta3);
	}
	else{
		sin_theta3 = 0.0;
		cos_theta3 = 1.0; 
	}

	//This function applies rotation in matrix M1, N1, and M2 when at least one node in the element has boundary conditions of no penetrability
	if (tag){ 
		no_penetrability(sin_theta1, sin_theta2, sin_theta3, cos_theta1, cos_theta2, cos_theta3, alpha_dt, delta, deltaNMV, 
				Area, y23, y12, y31, x32, x21, x13, Ax, Ay, M1, N1, M2, R1, R2, Ue, dUe, UeB, dUeB); 
	}
}


