#include "EulerEquations.h"

void BC_theta_OK(int J1, int J2, int J3, NodeType *Node, double theta[3], double (*BC_theta)(double, double))
{
		if (Node[J1].Type[1]<0 || Node[J1].Type[2]<0){
			theta[0] = BC_theta(Node[J1].x, Node[J1].y);
		}
		
		if (Node[J2].Type[1]<0 || Node[J2].Type[2]<0){
			theta[1] = BC_theta(Node[J2].x, Node[J2].y);
		}

		if (Node[J3].Type[1]<0 || Node[J3].Type[2]<0){
			theta[2] = BC_theta(Node[J3].x, Node[J3].y);
		}
}

void BC_theta_NO(int J1, int J2, int J3, NodeType *Node, double theta[3], double (*BC_theta)(double, double))
{
	theta[0] = 0.;
	theta[1] = 0.;
	theta[2] = 0.;
}


