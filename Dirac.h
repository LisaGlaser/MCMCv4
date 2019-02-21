#ifndef __DIRAC_H__
#define __DIRAC_H__

#include "matrices.h"
#include "progParams.h"
//#include "Randomc/randomc.h"  
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <cmath>

#include "Randomc/stocc.h"                     // define random library classes




int moveA_raw(Eigen::MatrixXcd &M, int her, double p, int s);
int moveT_raw(Eigen::MatrixXcd &M, double p, int s);


class Dirac
{
	public:
//		Dirac();
		Dirac(programParams iniV);
		int size;
		int truesize;
		int dimgamma;
		int type;
		double gD2,gD22,gD4;
		Eigen::MatrixXcd D;
		Eigen::MatrixXcd L23,L13,L12,H123,H0,L1,L2,L3;

		double getS();
		void moveA(double pA);
		void moveT(double pA);
		void makeD();
		void printD(FILE *file);
		void printAll(char *file);
		
	private:
		double S;
		int moved;
		double action();
		double action01();
		double action10();
		double action02();
		double action11();
		double action20();
		double action03();
		double action031();
		double action13();
		double action130();
		int getdimgamma();
		
		void makeD01();
		void makeD10();
		void makeD02();
		void makeD11();
		void makeD20();
		void makeD03();
		void makeD031();
		void makeD13();
		void makeD130();
		
		void initial(int i,char *inifile);
		
};

#endif
