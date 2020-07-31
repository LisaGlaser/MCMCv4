/*
 * MCMC.cxx
 *
 * Copyright 2016 Lisa Glaser <glaser@nbi.ku.dk>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 *
 * I know I should have commented more, but it's too late now.
 *
 */

/// new improved version, now with an attempt at making classes.
#include "progParams.h"
#include "Dirac.h"
#include "matrices.h"
#include "Randomc/stocc.h"                     // define random library classes

#include <typeinfo>
#include <stdio.h>
#include <iostream>
#include <time.h> // guess what, this measures time!
#include <sys/time.h>
//#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include "Randomc/randomc.h"

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files,
// If compiled as a project then compile and link in these cpp files.
   // Include code for the chosen random number generator:
   #include "Randomc/rancombi.cpp"
   #include "Randomc/mersenne.cpp"
   #include "Randomc/mother.cpp"
  #include "Randomc/stoc1.cpp"                // random library source code
   // define system specific user interface:
   #include "Randomc/userintf.cpp"
#endif

#define DEBUGMC 0
#define SUPERBUGMC 0
#define DEBUG	0
#define STEST 0
#define DDBUG 0


#define GAUSS 0 // uncommment this to get back to ordinary moves // this variable also exists in Dirac.h
//#define MOREEV 0

using Eigen::MatrixXcd;


using namespace Eigen;
using namespace std;

// I only define the random number generator once, so that successive calls give successive random numbers

TRandomCombined<CRandomMersenne,CRandomMother> RanGen(time(NULL));

SelfAdjointEigenSolver<MatrixXcd> es;


int measure( FILE *outfile, programParams iniV, Dirac D)
{


	// Now what do I want to measure?

	// First ideas S, Tr(D) and the traces of the m

	// traces don't seem to be particularly usefull as observables.
	// especially not of the Dirac, since that is trace-less
	// the determinant isn't interesting either, as soon as there is one zero eigenvalue it's screwed up
	// at least for the Dirac

	complex<double> temp;
	double tID,tD;

if( iniV.measure==2)
{
	//SelfAdjointEigenSolver<MatrixXcd> es(D); is called above and globally available
	VectorXd eivals;

	es.compute(D.D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
	eivals =es.eigenvalues();



	for(int i=0;i<D.truesize;i++)
		{
		fprintf (outfile,	" %10g ", eivals(i));
		}

		fprintf (outfile,"\n" );

}
else if( iniV.measure==3)
{

	D.printD(outfile);


}
else
{
	cout<< "Sorry we don't have that many measurement options!"<<endl;
	return 1;
}

	return 0;

}


int MChain(programParams iniV)
{
	Dirac Dtemp(iniV);
	Dirac D(iniV);
	int sweep,k,i=0,tempering,ttemp,pretemp=1;
	double weightA,weightM,p,ar=0.,ptemp,n1=0.;

	Dtemp=D;


	// and everything I need for a simple recording function
	char *filename;
	FILE *actionF;
	FILE *measureF;

	filename= (char*) calloc (30,sizeof(char));
	strcpy(filename,"action_monitor.txt");

	actionF = fopen(filename, "w");
	measureF = fopen(iniV.outfile, "w");


	sweep=iniV.matrixsize*4; // only measure the observables after 1 sweep, which I define to be the size of the matrix
	//	while(k*sweep<=500){k+=1.;}
	k=1000;

	weightA=iniV.wmoveA;


	if(weightA==0.)
	{
		tempering=1; // at first we want to temper


			weightA=1./pow(iniV.matrixsize,1.5);
	}
	else
	{
		tempering=0;
	}



	// this function implements the monte carlo chain

	while(i<iniV.stepnumber*sweep+1 || tempering==1)
	{
	i++;
	Dtemp=D;
	// then do a move

		Dtemp.moveA(weightA);
		if(iniV.moveT==1){ Dtemp.moveT(weightA);}

	// count moves
		n1+=1.;




	ptemp=exp(D.getS()-Dtemp.getS());
	p=RanGen.Random();
	if(DEBUGMC)	printf("The old action is %g the new action is %g\n ",D.getS(),Dtemp.getS());

	if(SUPERBUGMC) printf("Move happens if %g > %g \n",ptemp, p);
	if(ptemp>p)
	{
		D=Dtemp;
		ar+=1.;

	if(DEBUGMC) printf("We got a move! \n");

	}


	// action only get's recorded in a special file


	fprintf(actionF," %10f %10f %10f \n",D.getS(),ar/n1,weightA);

	// then measure the observables
	if(i%sweep== 0 && tempering==0)
	{
			ttemp=measure(measureF,iniV,D);

			if(ttemp==1) // check if the measurement failed, if it did throw me out.
			{
				return 1;
			}

	}
	else if(i%k== 0 && tempering==1 )//&& pretemp > k*D.size ) // or adjust the weight if we are still tempering /// I should probably make independent tempering processes for my two moves
	{
				printf("%d -th try \n",i/k);


				if(fabs(ar/(n1)-0.5)>0.01)
				{
					weightA+=weightA*(-0.5+ar/(n1));

				printf("The acceptance rate is %f10 and the new move weightA is %g \n",ar/(n1),weightA);

				ar=0.;
				n1=0.;
				ttemp=0;
				}
				else
				{
					ttemp+=1;
					printf("The acceptance rate stays %f10 and the move weightA stays %g  \n",ar/(n1),weightA);
					if(ttemp>5)
						{
							tempering=0;
							printf("We are done tempering. %g is the final weightA  \n",weightA);
							i=0;
							}

				}



			}
/*		else if(pretemp<= k*D.size)
		{
		pretemp++;
		i=0;
		//printf(" %d out of %d",pretemp,k*D.size);
		}*/

	}






	//printf("The acceptance rate over the whole run is %g \n",ar/nsteps);

	D.printAll(iniV.finfile);

	fclose(actionF);

	fclose(measureF);
	free(filename);
//	*/
	return 0;
}

int main(int argc, char *argv[])
{
	// Some variables I will be using during the main loop
	programParams parameters;

  // The first thing I do is initialize my program. For that I read a bunch of things from an input file
	// the programm expects the argument to be an input file
		if(argc==1)
		{
			cout<<"This programm needs an input file as an argument."<<endl;
			return 1;
		}
		else
		{

/// A function to initialize everything
			parameters.initialize(argv[1]);

			parameters.announce();
			srand (time(NULL));


			MChain(parameters);

			return 0;
		}

}
