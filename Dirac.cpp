#include "Dirac.h"


/*Dirac::Dirac()
{
	// these are shitty defaults should have called the nicer constructor
	size=1;
	truesize=1;
	dimgamma=1;
	type=1;
	D<<1;
	S=0; //action(0,0,0); can't call action  bc. it's virtual

} */


Dirac::Dirac(programParams iniV)
{
	type=iniV.Type;
	dimgamma=getdimgamma();
	size=iniV.matrixsize;
	truesize=dimgamma*size*size;
	gD2=iniV.gD2;
	gD22=iniV.gD22;
	gD4=iniV.gD4;

	initial(iniV.initialconfig,iniV.inifile);
	makeD();
	S=action();
	moved=0;
}

// Defines a function that sets the dimension of the Clifford module for each of the Clifford Types currently available.
int Dirac::getdimgamma()
{
	int dg;

		if(type==01)
		{
			dg=1;
		}
		else if(type==10)
		{
			dg=1;
		}
		else if(type==20)
		{
			dg=2;
		}
		else if(type==11)
		{
			dg=2;
		}
		else if(type==02)
		{
			dg=2;
		}
		else if(type==03)
		{
			dg=2;
		}
		else if(type==031)
		{
			dg=2;
		}
		else  if(type==13)
		{
			dg=4;
		}
		else if(type==130)
		{
			dg=4;
		}
		else
		{
			dg=1;
		}


	return dg;
}


// A function that initialises vital variables
void Dirac::initial(int ini,char *inifile)
{


	Eigen::MatrixXcd mtemp(size,size);

	if(ini==1) // the fuzzy 2 sphere
	{
		mtemp.setIdentity();

		// Get the necessary representations of su(2) generators
		L23=Lxy(size,1);
		L13=Lxy(size,2);
		L12=Lxy(size,3);
		H0=mtemp;
		H123=mtemp;

		mtemp=Eigen::MatrixXcd::Zero(size,size);

		L1=mtemp;
		L2=mtemp;
		L3=mtemp;

	}
	else if(ini==2) // an identity Dirac operator
	{
		mtemp=Eigen::MatrixXcd::Zero(size,size);

		L23=mtemp;
		L13=mtemp;
		L12=mtemp;
		H123=mtemp;
		L1=mtemp;
		L2=mtemp;
		L3=mtemp;
		mtemp.setIdentity();
		H0=mtemp;

	}
	else if(ini==3) // identity matrices ( or I * identity for the anti hermitian matrices)
	{
		mtemp.setIdentity();
		std::complex <double> I(0,1);

		L23=I*mtemp;
		L13=I*mtemp;
		L12=I*mtemp;
		H123=mtemp;
		L1=I*mtemp;
		L2=I*mtemp;
		L3=I*mtemp;
		H0=mtemp;

	}
/*
Paul: I've commented this out, because it causes errors to do with the file scan, something about not being able to
take the address of a rvalue. It might be something to do with .real and .imag, as they aren't assignable, so might not have an address? Unsure.
Code now compiles
*/

// 	else if(ini==4) // initial config from a file
// 	{
// 		L23.resize(size,size);
// 		L13.resize(size,size);
// 		L12.resize(size,size);
// 		H123.resize(size,size);
// 		H0.resize(size,size);
// 		L1.resize(size,size);
// 		L2.resize(size,size);
// 		L3.resize(size,size);
//
// 		FILE *eF;
// 		eF = fopen (inifile, "r");
//
// 		int j=0;
// //		j=fscanf(inifile, "%d \n", matrixsize);
// 				// I am using that k, k2, i, j are integers!
// //				mp(k,k2)=a(k/j,k2/j)*b(k%j,k2%j);
// 		for( int m=0;m<size*size;m++)
// 		{
//
// 	//	cout << m/n <<" =i" << m%n << "=j" <<endl;
// 		j=fscanf(eF, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n", 	&L23(m/size,m%size).real(),
// 																							&L23(m/size,m%size).imag(),
// 																							&L13(m/size,m%size).real(),
// 																							&L13(m/size,m%size).imag(),
// 																							&L12(m/size,m%size).real(),
// 																							&L12(m/size,m%size).imag(),
// 																							&H123(m/size,m%size).real(),
// 																							&H123(m/size,m%size).imag(),
// 																							&H0(m/size,m%size).real(),
// 																							&H0(m/size,m%size).imag(),
// 																							&L1(m/size,m%size).real(),
// 																							&L1(m/size,m%size).imag(),
// 																							&L2(m/size,m%size).real(),
// 																							&L2(m/size,m%size).imag(),
// 																							&L3(m/size,m%size).real(),
// 																							&L3(m/size,m%size).imag());
//
// 		if(j==-1)
// 		{
// 			printf("There was a problem with the file in line %d",m);
// 		}
//
// 		}
	//
	//
	// fclose(eF);
	// }
	else if(ini==5) // random Dirac operator
	{
		mtemp=Eigen::MatrixXcd::Random(size,size);

		L23=mtemp-mtemp.adjoint();

		mtemp=Eigen::MatrixXcd::Random(size,size);
		L13=mtemp-mtemp.adjoint();

		mtemp=Eigen::MatrixXcd::Random(size,size);
		L12=mtemp-mtemp.adjoint();

		mtemp=Eigen::MatrixXcd::Random(size,size);
		H123=mtemp+mtemp.adjoint();

		mtemp=Eigen::MatrixXcd::Random(size,size);
		L1=mtemp-mtemp.adjoint();

		mtemp=Eigen::MatrixXcd::Random(size,size);
		L2=mtemp-mtemp.adjoint();

		mtemp=Eigen::MatrixXcd::Random(size,size);
		L3=mtemp-mtemp.adjoint();

		mtemp=Eigen::MatrixXcd::Random(size,size);
		H0=mtemp+mtemp.adjoint();
	}
	else
	{
		std::cout<< " Sorry we don't support that many initial configurations yet"<<std::endl;
		}


}

void Dirac::printD(FILE *file)
{

	for(int i=0;i<truesize;i++)
		{
			for(int j=0;j<truesize;j++)
			{

			fprintf (file, "%10f+I%10f ", D(i,j).real(), D(i,j).imag() );

			}

		}

		fprintf (file,"\n" );
}

void Dirac::printAll(char *file)
{

	FILE *eF;
	eF = fopen (file, "w");


	for(int i=0; i<size; i++)
	{
		for(int k=0; k<size; k++)
		{


			fprintf(eF, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n",L23(i,k).real(),
																							L23(i,k).imag(),
																							L13(i,k).real(),
																							L13(i,k).imag(),
																							L12(i,k).real(),
																							L12(i,k).imag(),
																							H123(i,k).real(),
																							H123(i,k).imag(),
																							H0(i,k).real(),
																							H0(i,k).imag(),
																							L1(i,k).real(),
																							L1(i,k).imag(),
																							L2(i,k).real(),
																							L2(i,k).imag(),
																							L3(i,k).real(),
																							L3(i,k).imag());

				if(DEBUG) printf(" at element %d %d \n", i,k);
		}

	}
//	fprintf(thefile, L23.real());

	fclose(eF);

}

double Dirac::getS()
{
	if(moved==1)
	{
		if(type==01)
		{
			S=action01();
		}
		else if(type==10)
		{
			S=action10();
		}
		else if(type==20)
		{
			S=action20();
		}
		else if(type==11)
		{
			S=action11();
		}
		else if(type==02)
		{
			S=action02();
		}
		else if(type==03)
		{
			S=action03();
		}
		else if(type==031)
		{
			S=action031();
		}
		else  if(type==13)
		{
			S=action13();
		}
		else if(type==130)
		{
			S=action130();
		}
		else
		{
			S=action();
		}
	moved=0;
	}


	return S;

}

double Dirac::action()
{
	std::complex<double> S2,S4;
	double S;
	// New simpler action, we just want D^2

	Eigen::MatrixXcd temp1=D*D;
	Eigen::MatrixXcd temp2=temp1*temp1;
	// in GR the action should be additive, and that's why we add the extrinsic curvature. But what about here?

	S2=(temp1).trace();
	S4=0;
	if(gD4!=0)
	{
		S4= (temp2).trace();
	}
	if(S2.imag()>10e-10 || S4.imag()>10e-10)
	{
		std::cout<<"This sucks, your actionG isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}

	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();


	return S;
}


void Dirac::moveA(double pA)
{
int hermitian;
Eigen::MatrixXcd temp;

	hermitian=-1;
	moveA_raw(L23,hermitian,pA,size);
	if(!L23.isApprox(-1.*L23.adjoint())) std::cout<<"L23 warning"<<std::endl;
	temp=L23-L23.adjoint();
	L23=0.5*temp;

	moveA_raw(L13,hermitian,pA,size);
	if(!L13.isApprox(-1.*L13.adjoint())) std::cout<<"L13 warning"<<std::endl;
	temp=L13-L13.adjoint();
	L13=0.5*temp;

	moveA_raw(L12,hermitian,pA,size);
	if(!L12.isApprox(-1.*L12.adjoint())) std::cout<<"L12 warning"<<std::endl;
	temp=L12-L12.adjoint();
	L12=0.5*temp;

	moveA_raw(L1,hermitian,pA,size);
	if(!L1.isApprox(-1.*L1.adjoint())) std::cout<<"L1 warning"<<std::endl;
	temp=L1-L1.adjoint();
	L1=0.5*temp;

	moveA_raw(L2,hermitian,pA,size);
	if(!L2.isApprox(-1.*L2.adjoint())) std::cout<<"L2 warning"<<std::endl;
	temp=L2-L2.adjoint();
	L2=0.5*temp;

	moveA_raw(L3,hermitian,pA,size);
	if(!L3.isApprox(-1.*L3.adjoint())) std::cout<<"L3 warning"<<std::endl;
	temp=L3-L3.adjoint();
	L3=0.5*temp;


	hermitian=1;
	moveA_raw(H123,hermitian,pA,size);
	if(!H123.isApprox(H123.adjoint())) std::cout<<"H123 warning"<<std::endl;
	temp=H123+H123.adjoint();
	H123=0.5*temp;

	moveA_raw(H0,hermitian,pA,size);
	if(!H0.isApprox(H0.adjoint())) std::cout<<"H0 warning"<<std::endl;
	temp=H0+H0.adjoint();
	H0=0.5*temp;

	makeD();
	moved=1;
	if(!D.isApprox(D.adjoint()))
	{
		std::cout<<"That Dirac is not selfadjoint after moveM, something is wrong here!"<<std::endl;
		std::cout<<D-D.adjoint()<<std::endl;
		}

}


void Dirac::moveT(double pA)
{

Eigen::MatrixXcd temp;

	moveT_raw(H123,pA,size);
	moveT_raw(H0,pA,size);

	makeD();
	moved=1;
	if(!D.isApprox(D.adjoint()))
	{
		std::cout<<"That Dirac is not selfadjoint after moveM, something is wrong here!"<<std::endl;
		std::cout<<D-D.adjoint()<<std::endl;
		}

}



int seed = (int)time(0);    // random seed
StochasticLib1 sto(seed);           // make instance of random library



int moveA_raw(Eigen::MatrixXcd &M, int her, double p,int s)
{
	Eigen::MatrixXcd diff(s,s);
	diff=Eigen::MatrixXcd::Random(s,s);

	M+=sto.Normal(p,p)*(diff+ her*diff.adjoint()); // here comes the anti hermitian


	return 0;

}




int moveT_raw(Eigen::MatrixXcd &M, double p,int s)
{
	Eigen::MatrixXcd diff(s,s);
	diff.setIdentity();

	M+=sto.Normal(0.,p)*M.trace()*diff.setIdentity()/s; // here comes the anti hermitian


	return 0;

}




void Dirac::makeD10()
{

	D= acom(H0);


}


void Dirac::makeD01()
{

	D= com(H0);


}

// A small routine to caluclate the Dirac operator
// implicit assumptions are that L23,L13,L12,H0 have the same size and that they are squared
void Dirac::makeD130()
{

	Eigen::MatrixXcd g0,g1,g2,g3;

	g0=gamma13(0);
	g1=gamma13(1);
	g2=gamma13(2);
	g3=gamma13(3);



	D= tens(g0*g1*g2,com(L12)) + tens(g0*g1*g3,com(L13)) + tens(g0*g2*g3,com(L23)) + tens(g0,acom(H0));



}


// A small routine to calculate the general Dirac operator
// implicit assumptions are that L23,L13,L12,H0 have the same size and that they are squared
void Dirac::makeD13()
{

	Eigen::MatrixXcd g0,g1,g2,g3;

	g0=gamma13(0);
	g1=gamma13(1);
	g2=gamma13(2);
	g3=gamma13(3);




	D= tens(g0*g1*g2,com(L12)) + tens(g0*g1*g3,com(L13)) + tens(g0*g2*g3,com(L23))+ tens(g1*g2*g3,acom(H123))
			+ tens(g0,acom(H0))+ tens(g1,com(L1))+ tens(g2,com(L2))+ tens(g3,com(L3));



}


void Dirac::makeD20()
{

	Eigen::MatrixXcd g1,g2;
	g1=gamma20(1);
	g2=gamma20(2);

	D=tens(g1,acom(H0))+tens(g2,acom(H123));

}

void Dirac::makeD11()
{

	Eigen::MatrixXcd g1,g2;
	g1=gamma11(1);
	g2=gamma11(2);

	D=tens(g1,acom(H0))+tens(g2,com(L1));



}

void Dirac::makeD02()
{

	Eigen::MatrixXcd g1,g2;
	g1=gamma02(1);
	g2=gamma02(2);

	D=tens(g1,com(L1))+tens(g2,com(L2));


}

void Dirac::makeD03()
{

	Eigen::MatrixXcd g1,g2,g3;
	g1=gamma03(1);
	g2=gamma03(2);
	g3=gamma03(3);
	Eigen::MatrixXcd ident(g1.cols(),g1.cols());
	Eigen::MatrixXcd temp;
	ident.setIdentity();

	D=tens(ident,acom(H0))+tens(g1,com(L1))+tens(g2,com(L2))+tens(g3,com(L3));



}


void Dirac::makeD031()
{

	Eigen::MatrixXcd g1,g2,g3;
	g1=gamma03(1);
	g2=gamma03(2);
	g3=gamma03(3);

	D=tens(g1,com(L1))+tens(g2,com(L2))+tens(g3,com(L3));



}


void Dirac::makeD()
{


	if(type==13)
	{
		makeD13();
	}
	else if(type==130)
	{
		makeD130();
	}
	else if(type==01)
	{
		makeD01();
	}
	else if(type==10)
	{
		makeD10();
	}
	else if(type==20)
	{
		makeD20();
	}
	else if(type==11)
	{
		makeD11();
	}
	else if(type==02)
	{
		makeD02();
	}
	else if(type==31)
	{
		makeD031();
	}
	else if(type==03)
	{
		makeD03();
	}
	else
	{
		printf("Sorry, we haven't implemented Type %d",type);
	}


}


double Dirac::action130()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd temp1=H0*H0;
	Eigen::MatrixXcd temp2=L23*L23+L13*L13+L12*L12;
	//Eigen::MatrixXcd temp2=temp1*temp1;
	// in GR the action should be additive, and that's why we add the extrinsic curvature. But what about here?

	S2= 8.*size*temp1.trace()+8.*H0.trace()*H0.trace()-8.*size*temp2.trace()+8.*(L23.trace()*L23.trace()+L13.trace()*L13.trace()+L12.trace()*L12.trace());
//	S4= (temp2).trace();

	if(S2.imag()>10e-10 ) //|| S4.imag()>10e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
	}

//	std::cout<<"temp1"<< temp1.trace()<<std::endl;
//	std::cout<<"temp2"<< temp2.trace()<<std::endl;

	S= gD2*S2.real()+gD22*S2.real()*S2.real(); //+0.5*S4.real();


	return S;
}


double Dirac::action13()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd temp1=H0*H0-L1*L1-L2*L2-L3*L3;
	Eigen::MatrixXcd temp2=L23*L23+L13*L13+L12*L12-H123*H123;

	S2= 8.*size*temp1.trace()+8.*H0.trace()*H0.trace()+8.*L1.trace()*L1.trace()+8.*L2.trace()*L2.trace()+8.*L3.trace()*L3.trace()-8.*size*temp2.trace()+8.*(L23.trace()*L23.trace()+L13.trace()*L13.trace()+L12.trace()*L12.trace()+H123.trace()*H123.trace());
	S4=0;
	if(gD4!=0)
	{
// Mathematica generated code. I'm not gonna type all this shit in by hand :P
Eigen::MatrixXcd H0H0=H0*H0; std::complex<double> H0H0t=H0H0.trace();
//Eigen::MatrixXcd H0=H0;
std::complex<double> H0t=H0.trace();
Eigen::MatrixXcd H0H0H0=H0*H0*H0; std::complex<double> H0H0H0t=H0H0H0.trace();
Eigen::MatrixXcd H0H0H0H0=H0*H0*H0*H0; std::complex<double> H0H0H0H0t=H0H0H0H0.trace();
Eigen::MatrixXcd L1L1=L1*L1; std::complex<double> L1L1t=L1L1.trace();
//Eigen::MatrixXcd L1=L1;
std::complex<double> L1t=L1.trace();
Eigen::MatrixXcd L1L1L1=L1*L1*L1; std::complex<double> L1L1L1t=L1L1L1.trace();
Eigen::MatrixXcd L1L1L1L1=L1*L1*L1*L1; std::complex<double> L1L1L1L1t=L1L1L1L1.trace();
Eigen::MatrixXcd L12L12=L12*L12; std::complex<double> L12L12t=L12L12.trace();
//Eigen::MatrixXcd L12=L12;
std::complex<double> L12t=L12.trace();
Eigen::MatrixXcd L12L12L12=L12*L12*L12; std::complex<double> L12L12L12t=L12L12L12.trace();
Eigen::MatrixXcd L12L12L12L12=L12*L12*L12*L12; std::complex<double> L12L12L12L12t=L12L12L12L12.trace();
Eigen::MatrixXcd H123H123=H123*H123; std::complex<double> H123H123t=H123H123.trace();
//Eigen::MatrixXcd H123=H123;
std::complex<double> H123t=H123.trace();
Eigen::MatrixXcd H123H123H123=H123*H123*H123; std::complex<double> H123H123H123t=H123H123H123.trace();
Eigen::MatrixXcd H123H123H123H123=H123*H123*H123*H123; std::complex<double> H123H123H123H123t=H123H123H123H123.trace();
Eigen::MatrixXcd L13L13=L13*L13; std::complex<double> L13L13t=L13L13.trace();
//Eigen::MatrixXcd L13=L13;
std::complex<double> L13t=L13.trace();
Eigen::MatrixXcd L13L13L13=L13*L13*L13; std::complex<double> L13L13L13t=L13L13L13.trace();
Eigen::MatrixXcd L13L13L13L13=L13*L13*L13*L13; std::complex<double> L13L13L13L13t=L13L13L13L13.trace();
Eigen::MatrixXcd L2L2=L2*L2; std::complex<double> L2L2t=L2L2.trace();
//Eigen::MatrixXcd L2=L2;
std::complex<double> L2t=L2.trace();
Eigen::MatrixXcd L2L2L2=L2*L2*L2; std::complex<double> L2L2L2t=L2L2L2.trace();
Eigen::MatrixXcd L2L2L2L2=L2*L2*L2*L2; std::complex<double> L2L2L2L2t=L2L2L2L2.trace();
Eigen::MatrixXcd L23L23=L23*L23; std::complex<double> L23L23t=L23L23.trace();
//Eigen::MatrixXcd L23=L23;
std::complex<double> L23t=L23.trace();
Eigen::MatrixXcd L23L23L23=L23*L23*L23; std::complex<double> L23L23L23t=L23L23L23.trace();
Eigen::MatrixXcd L23L23L23L23=L23*L23*L23*L23; std::complex<double> L23L23L23L23t=L23L23L23L23.trace();
Eigen::MatrixXcd L3L3=L3*L3; std::complex<double> L3L3t=L3L3.trace();
//Eigen::MatrixXcd L3=L3;
std::complex<double> L3t=L3.trace();
Eigen::MatrixXcd L3L3L3=L3*L3*L3; std::complex<double> L3L3L3t=L3L3L3.trace();
Eigen::MatrixXcd L3L3L3L3=L3*L3*L3*L3; std::complex<double> L3L3L3L3t=L3L3L3L3.trace();
Eigen::MatrixXcd H0L1=H0*L1; std::complex<double> H0L1t=H0L1.trace();
Eigen::MatrixXcd H0L1H0L1=H0*L1*H0*L1; std::complex<double> H0L1H0L1t=H0L1H0L1.trace();
Eigen::MatrixXcd H0L1L1=H0*L1*L1; std::complex<double> H0L1L1t=H0L1L1.trace();
Eigen::MatrixXcd H0L12=H0*L12; std::complex<double> H0L12t=H0L12.trace();
Eigen::MatrixXcd H0L12H0L12=H0*L12*H0*L12; std::complex<double> H0L12H0L12t=H0L12H0L12.trace();
Eigen::MatrixXcd H0L12L12=H0*L12*L12; std::complex<double> H0L12L12t=H0L12L12.trace();
Eigen::MatrixXcd H0H123=H0*H123; std::complex<double> H0H123t=H0H123.trace();
Eigen::MatrixXcd H0H123H0H123=H0*H123*H0*H123; std::complex<double> H0H123H0H123t=H0H123H0H123.trace();
Eigen::MatrixXcd H0H123H123=H0*H123*H123; std::complex<double> H0H123H123t=H0H123H123.trace();
Eigen::MatrixXcd H0L13=H0*L13; std::complex<double> H0L13t=H0L13.trace();
Eigen::MatrixXcd H0L13H0L13=H0*L13*H0*L13; std::complex<double> H0L13H0L13t=H0L13H0L13.trace();
Eigen::MatrixXcd H0L13L13=H0*L13*L13; std::complex<double> H0L13L13t=H0L13L13.trace();
Eigen::MatrixXcd H0L2=H0*L2; std::complex<double> H0L2t=H0L2.trace();
Eigen::MatrixXcd H0L2H0L2=H0*L2*H0*L2; std::complex<double> H0L2H0L2t=H0L2H0L2.trace();
Eigen::MatrixXcd H0L2L2=H0*L2*L2; std::complex<double> H0L2L2t=H0L2L2.trace();
Eigen::MatrixXcd H0L23=H0*L23; std::complex<double> H0L23t=H0L23.trace();
Eigen::MatrixXcd H0L23H0L23=H0*L23*H0*L23; std::complex<double> H0L23H0L23t=H0L23H0L23.trace();
Eigen::MatrixXcd H0L23L23=H0*L23*L23; std::complex<double> H0L23L23t=H0L23L23.trace();
Eigen::MatrixXcd H0L3=H0*L3; std::complex<double> H0L3t=H0L3.trace();
Eigen::MatrixXcd H0L3H0L3=H0*L3*H0*L3; std::complex<double> H0L3H0L3t=H0L3H0L3.trace();
Eigen::MatrixXcd H0L3L3=H0*L3*L3; std::complex<double> H0L3L3t=H0L3L3.trace();
Eigen::MatrixXcd H0H0L1=H0*H0*L1; std::complex<double> H0H0L1t=H0H0L1.trace();
Eigen::MatrixXcd H0H0L1L1=H0*H0*L1*L1; std::complex<double> H0H0L1L1t=H0H0L1L1.trace();
Eigen::MatrixXcd H0H0L12=H0*H0*L12; std::complex<double> H0H0L12t=H0H0L12.trace();
Eigen::MatrixXcd H0H0L12L12=H0*H0*L12*L12; std::complex<double> H0H0L12L12t=H0H0L12L12.trace();
Eigen::MatrixXcd H0H0H123=H0*H0*H123; std::complex<double> H0H0H123t=H0H0H123.trace();
Eigen::MatrixXcd H0H0H123H123=H0*H0*H123*H123; std::complex<double> H0H0H123H123t=H0H0H123H123.trace();
Eigen::MatrixXcd H0H0L13=H0*H0*L13; std::complex<double> H0H0L13t=H0H0L13.trace();
Eigen::MatrixXcd H0H0L13L13=H0*H0*L13*L13; std::complex<double> H0H0L13L13t=H0H0L13L13.trace();
Eigen::MatrixXcd H0H0L2=H0*H0*L2; std::complex<double> H0H0L2t=H0H0L2.trace();
Eigen::MatrixXcd H0H0L2L2=H0*H0*L2*L2; std::complex<double> H0H0L2L2t=H0H0L2L2.trace();
Eigen::MatrixXcd H0H0L23=H0*H0*L23; std::complex<double> H0H0L23t=H0H0L23.trace();
Eigen::MatrixXcd H0H0L23L23=H0*H0*L23*L23; std::complex<double> H0H0L23L23t=H0H0L23L23.trace();
Eigen::MatrixXcd H0H0L3=H0*H0*L3; std::complex<double> H0H0L3t=H0H0L3.trace();
Eigen::MatrixXcd H0H0L3L3=H0*H0*L3*L3; std::complex<double> H0H0L3L3t=H0H0L3L3.trace();
Eigen::MatrixXcd L1L12=L1*L12; std::complex<double> L1L12t=L1L12.trace();
Eigen::MatrixXcd L1L12L1L12=L1*L12*L1*L12; std::complex<double> L1L12L1L12t=L1L12L1L12.trace();
Eigen::MatrixXcd L1L12L12=L1*L12*L12; std::complex<double> L1L12L12t=L1L12L12.trace();
Eigen::MatrixXcd L1H123=L1*H123; std::complex<double> L1H123t=L1H123.trace();
Eigen::MatrixXcd L1H123L1H123=L1*H123*L1*H123; std::complex<double> L1H123L1H123t=L1H123L1H123.trace();
Eigen::MatrixXcd L1H123H123=L1*H123*H123; std::complex<double> L1H123H123t=L1H123H123.trace();
Eigen::MatrixXcd L1L13=L1*L13; std::complex<double> L1L13t=L1L13.trace();
Eigen::MatrixXcd L1L13L1L13=L1*L13*L1*L13; std::complex<double> L1L13L1L13t=L1L13L1L13.trace();
Eigen::MatrixXcd L1L13L13=L1*L13*L13; std::complex<double> L1L13L13t=L1L13L13.trace();
Eigen::MatrixXcd L1L2=L1*L2; std::complex<double> L1L2t=L1L2.trace();
Eigen::MatrixXcd L1L2L1L2=L1*L2*L1*L2; std::complex<double> L1L2L1L2t=L1L2L1L2.trace();
Eigen::MatrixXcd L1L2L2=L1*L2*L2; std::complex<double> L1L2L2t=L1L2L2.trace();
Eigen::MatrixXcd L1L23=L1*L23; std::complex<double> L1L23t=L1L23.trace();
Eigen::MatrixXcd L1L23L1L23=L1*L23*L1*L23; std::complex<double> L1L23L1L23t=L1L23L1L23.trace();
Eigen::MatrixXcd L1L23L23=L1*L23*L23; std::complex<double> L1L23L23t=L1L23L23.trace();
Eigen::MatrixXcd L1L3=L1*L3; std::complex<double> L1L3t=L1L3.trace();
Eigen::MatrixXcd L1L3L1L3=L1*L3*L1*L3; std::complex<double> L1L3L1L3t=L1L3L1L3.trace();
Eigen::MatrixXcd L1L3L3=L1*L3*L3; std::complex<double> L1L3L3t=L1L3L3.trace();
Eigen::MatrixXcd L1L1L12=L1*L1*L12; std::complex<double> L1L1L12t=L1L1L12.trace();
Eigen::MatrixXcd L1L1L12L12=L1*L1*L12*L12; std::complex<double> L1L1L12L12t=L1L1L12L12.trace();
Eigen::MatrixXcd L1L1H123=L1*L1*H123; std::complex<double> L1L1H123t=L1L1H123.trace();
Eigen::MatrixXcd L1L1H123H123=L1*L1*H123*H123; std::complex<double> L1L1H123H123t=L1L1H123H123.trace();
Eigen::MatrixXcd L1L1L13=L1*L1*L13; std::complex<double> L1L1L13t=L1L1L13.trace();
Eigen::MatrixXcd L1L1L13L13=L1*L1*L13*L13; std::complex<double> L1L1L13L13t=L1L1L13L13.trace();
Eigen::MatrixXcd L1L1L2=L1*L1*L2; std::complex<double> L1L1L2t=L1L1L2.trace();
Eigen::MatrixXcd L1L1L2L2=L1*L1*L2*L2; std::complex<double> L1L1L2L2t=L1L1L2L2.trace();
Eigen::MatrixXcd L1L1L23=L1*L1*L23; std::complex<double> L1L1L23t=L1L1L23.trace();
Eigen::MatrixXcd L1L1L23L23=L1*L1*L23*L23; std::complex<double> L1L1L23L23t=L1L1L23L23.trace();
Eigen::MatrixXcd L1L1L3=L1*L1*L3; std::complex<double> L1L1L3t=L1L1L3.trace();
Eigen::MatrixXcd L1L1L3L3=L1*L1*L3*L3; std::complex<double> L1L1L3L3t=L1L1L3L3.trace();
Eigen::MatrixXcd L12H123=L12*H123; std::complex<double> L12H123t=L12H123.trace();
Eigen::MatrixXcd L12H123L12H123=L12*H123*L12*H123; std::complex<double> L12H123L12H123t=L12H123L12H123.trace();
Eigen::MatrixXcd L12H123H123=L12*H123*H123; std::complex<double> L12H123H123t=L12H123H123.trace();
Eigen::MatrixXcd L12L13=L12*L13; std::complex<double> L12L13t=L12L13.trace();
Eigen::MatrixXcd L12L13L12L13=L12*L13*L12*L13; std::complex<double> L12L13L12L13t=L12L13L12L13.trace();
Eigen::MatrixXcd L12L13L13=L12*L13*L13; std::complex<double> L12L13L13t=L12L13L13.trace();
Eigen::MatrixXcd L12L2=L12*L2; std::complex<double> L12L2t=L12L2.trace();
Eigen::MatrixXcd L12L2L12L2=L12*L2*L12*L2; std::complex<double> L12L2L12L2t=L12L2L12L2.trace();
Eigen::MatrixXcd L12L2L2=L12*L2*L2; std::complex<double> L12L2L2t=L12L2L2.trace();
Eigen::MatrixXcd L12L23=L12*L23; std::complex<double> L12L23t=L12L23.trace();
Eigen::MatrixXcd L12L23L12L23=L12*L23*L12*L23; std::complex<double> L12L23L12L23t=L12L23L12L23.trace();
Eigen::MatrixXcd L12L23L23=L12*L23*L23; std::complex<double> L12L23L23t=L12L23L23.trace();
Eigen::MatrixXcd L12L3=L12*L3; std::complex<double> L12L3t=L12L3.trace();
Eigen::MatrixXcd L12L3L12L3=L12*L3*L12*L3; std::complex<double> L12L3L12L3t=L12L3L12L3.trace();
Eigen::MatrixXcd L12L3L3=L12*L3*L3; std::complex<double> L12L3L3t=L12L3L3.trace();
Eigen::MatrixXcd L12L12H123=L12*L12*H123; std::complex<double> L12L12H123t=L12L12H123.trace();
Eigen::MatrixXcd L12L12H123H123=L12*L12*H123*H123; std::complex<double> L12L12H123H123t=L12L12H123H123.trace();
Eigen::MatrixXcd L12L12L13=L12*L12*L13; std::complex<double> L12L12L13t=L12L12L13.trace();
Eigen::MatrixXcd L12L12L13L13=L12*L12*L13*L13; std::complex<double> L12L12L13L13t=L12L12L13L13.trace();
Eigen::MatrixXcd L12L12L2=L12*L12*L2; std::complex<double> L12L12L2t=L12L12L2.trace();
Eigen::MatrixXcd L12L12L2L2=L12*L12*L2*L2; std::complex<double> L12L12L2L2t=L12L12L2L2.trace();
Eigen::MatrixXcd L12L12L23=L12*L12*L23; std::complex<double> L12L12L23t=L12L12L23.trace();
Eigen::MatrixXcd L12L12L23L23=L12*L12*L23*L23; std::complex<double> L12L12L23L23t=L12L12L23L23.trace();
Eigen::MatrixXcd L12L12L3=L12*L12*L3; std::complex<double> L12L12L3t=L12L12L3.trace();
Eigen::MatrixXcd L12L12L3L3=L12*L12*L3*L3; std::complex<double> L12L12L3L3t=L12L12L3L3.trace();
Eigen::MatrixXcd H123L13=H123*L13; std::complex<double> H123L13t=H123L13.trace();
Eigen::MatrixXcd H123L13H123L13=H123*L13*H123*L13; std::complex<double> H123L13H123L13t=H123L13H123L13.trace();
Eigen::MatrixXcd H123L13L13=H123*L13*L13; std::complex<double> H123L13L13t=H123L13L13.trace();
Eigen::MatrixXcd H123L2=H123*L2; std::complex<double> H123L2t=H123L2.trace();
Eigen::MatrixXcd H123L2H123L2=H123*L2*H123*L2; std::complex<double> H123L2H123L2t=H123L2H123L2.trace();
Eigen::MatrixXcd H123L2L2=H123*L2*L2; std::complex<double> H123L2L2t=H123L2L2.trace();
Eigen::MatrixXcd H123L23=H123*L23; std::complex<double> H123L23t=H123L23.trace();
Eigen::MatrixXcd H123L23H123L23=H123*L23*H123*L23; std::complex<double> H123L23H123L23t=H123L23H123L23.trace();
Eigen::MatrixXcd H123L23L23=H123*L23*L23; std::complex<double> H123L23L23t=H123L23L23.trace();
Eigen::MatrixXcd H123L3=H123*L3; std::complex<double> H123L3t=H123L3.trace();
Eigen::MatrixXcd H123L3H123L3=H123*L3*H123*L3; std::complex<double> H123L3H123L3t=H123L3H123L3.trace();
Eigen::MatrixXcd H123L3L3=H123*L3*L3; std::complex<double> H123L3L3t=H123L3L3.trace();
Eigen::MatrixXcd H123H123L13=H123*H123*L13; std::complex<double> H123H123L13t=H123H123L13.trace();
Eigen::MatrixXcd H123H123L13L13=H123*H123*L13*L13; std::complex<double> H123H123L13L13t=H123H123L13L13.trace();
Eigen::MatrixXcd H123H123L2=H123*H123*L2; std::complex<double> H123H123L2t=H123H123L2.trace();
Eigen::MatrixXcd H123H123L2L2=H123*H123*L2*L2; std::complex<double> H123H123L2L2t=H123H123L2L2.trace();
Eigen::MatrixXcd H123H123L23=H123*H123*L23; std::complex<double> H123H123L23t=H123H123L23.trace();
Eigen::MatrixXcd H123H123L23L23=H123*H123*L23*L23; std::complex<double> H123H123L23L23t=H123H123L23L23.trace();
Eigen::MatrixXcd H123H123L3=H123*H123*L3; std::complex<double> H123H123L3t=H123H123L3.trace();
Eigen::MatrixXcd H123H123L3L3=H123*H123*L3*L3; std::complex<double> H123H123L3L3t=H123H123L3L3.trace();
Eigen::MatrixXcd L13L2=L13*L2; std::complex<double> L13L2t=L13L2.trace();
Eigen::MatrixXcd L13L2L13L2=L13*L2*L13*L2; std::complex<double> L13L2L13L2t=L13L2L13L2.trace();
Eigen::MatrixXcd L13L2L2=L13*L2*L2; std::complex<double> L13L2L2t=L13L2L2.trace();
Eigen::MatrixXcd L13L23=L13*L23; std::complex<double> L13L23t=L13L23.trace();
Eigen::MatrixXcd L13L23L13L23=L13*L23*L13*L23; std::complex<double> L13L23L13L23t=L13L23L13L23.trace();
Eigen::MatrixXcd L13L23L23=L13*L23*L23; std::complex<double> L13L23L23t=L13L23L23.trace();
Eigen::MatrixXcd L13L3=L13*L3; std::complex<double> L13L3t=L13L3.trace();
Eigen::MatrixXcd L13L3L13L3=L13*L3*L13*L3; std::complex<double> L13L3L13L3t=L13L3L13L3.trace();
Eigen::MatrixXcd L13L3L3=L13*L3*L3; std::complex<double> L13L3L3t=L13L3L3.trace();
Eigen::MatrixXcd L13L13L2=L13*L13*L2; std::complex<double> L13L13L2t=L13L13L2.trace();
Eigen::MatrixXcd L13L13L2L2=L13*L13*L2*L2; std::complex<double> L13L13L2L2t=L13L13L2L2.trace();
Eigen::MatrixXcd L13L13L23=L13*L13*L23; std::complex<double> L13L13L23t=L13L13L23.trace();
Eigen::MatrixXcd L13L13L23L23=L13*L13*L23*L23; std::complex<double> L13L13L23L23t=L13L13L23L23.trace();
Eigen::MatrixXcd L13L13L3=L13*L13*L3; std::complex<double> L13L13L3t=L13L13L3.trace();
Eigen::MatrixXcd L13L13L3L3=L13*L13*L3*L3; std::complex<double> L13L13L3L3t=L13L13L3L3.trace();
Eigen::MatrixXcd L2L23=L2*L23; std::complex<double> L2L23t=L2L23.trace();
Eigen::MatrixXcd L2L23L2L23=L2*L23*L2*L23; std::complex<double> L2L23L2L23t=L2L23L2L23.trace();
Eigen::MatrixXcd L2L23L23=L2*L23*L23; std::complex<double> L2L23L23t=L2L23L23.trace();
Eigen::MatrixXcd L2L3=L2*L3; std::complex<double> L2L3t=L2L3.trace();
Eigen::MatrixXcd L2L3L2L3=L2*L3*L2*L3; std::complex<double> L2L3L2L3t=L2L3L2L3.trace();
Eigen::MatrixXcd L2L3L3=L2*L3*L3; std::complex<double> L2L3L3t=L2L3L3.trace();
Eigen::MatrixXcd L2L2L23=L2*L2*L23; std::complex<double> L2L2L23t=L2L2L23.trace();
Eigen::MatrixXcd L2L2L23L23=L2*L2*L23*L23; std::complex<double> L2L2L23L23t=L2L2L23L23.trace();
Eigen::MatrixXcd L2L2L3=L2*L2*L3; std::complex<double> L2L2L3t=L2L2L3.trace();
Eigen::MatrixXcd L2L2L3L3=L2*L2*L3*L3; std::complex<double> L2L2L3L3t=L2L2L3L3.trace();
Eigen::MatrixXcd L23L3=L23*L3; std::complex<double> L23L3t=L23L3.trace();
Eigen::MatrixXcd L23L3L23L3=L23*L3*L23*L3; std::complex<double> L23L3L23L3t=L23L3L23L3.trace();
Eigen::MatrixXcd L23L3L3=L23*L3*L3; std::complex<double> L23L3L3t=L23L3L3.trace();
Eigen::MatrixXcd L23L23L3=L23*L23*L3; std::complex<double> L23L23L3t=L23L23L3.trace();
Eigen::MatrixXcd L23L23L3L3=L23*L23*L3*L3; std::complex<double> L23L23L3L3t=L23L23L3L3.trace();
Eigen::MatrixXcd H0L1L12=H0*L1*L12; std::complex<double> H0L1L12t=H0L1L12.trace();
Eigen::MatrixXcd H0L1L13=H0*L1*L13; std::complex<double> H0L1L13t=H0L1L13.trace();
Eigen::MatrixXcd H0L1L2=H0*L1*L2; std::complex<double> H0L1L2t=H0L1L2.trace();
Eigen::MatrixXcd H0L1L23=H0*L1*L23; std::complex<double> H0L1L23t=H0L1L23.trace();
Eigen::MatrixXcd H0L1L3=H0*L1*L3; std::complex<double> H0L1L3t=H0L1L3.trace();
Eigen::MatrixXcd H0L12H123=H0*L12*H123; std::complex<double> H0L12H123t=H0L12H123.trace();
Eigen::MatrixXcd H0L12L2=H0*L12*L2; std::complex<double> H0L12L2t=H0L12L2.trace();
Eigen::MatrixXcd H0L12L3=H0*L12*L3; std::complex<double> H0L12L3t=H0L12L3.trace();
Eigen::MatrixXcd H0H123L13=H0*H123*L13; std::complex<double> H0H123L13t=H0H123L13.trace();
Eigen::MatrixXcd H0H123L23=H0*H123*L23; std::complex<double> H0H123L23t=H0H123L23.trace();
Eigen::MatrixXcd H0L13L2=H0*L13*L2; std::complex<double> H0L13L2t=H0L13L2.trace();
Eigen::MatrixXcd H0L13L3=H0*L13*L3; std::complex<double> H0L13L3t=H0L13L3.trace();
Eigen::MatrixXcd H0L2L23=H0*L2*L23; std::complex<double> H0L2L23t=H0L2L23.trace();
Eigen::MatrixXcd H0L2L3=H0*L2*L3; std::complex<double> H0L2L3t=H0L2L3.trace();
Eigen::MatrixXcd H0L23L3=H0*L23*L3; std::complex<double> H0L23L3t=H0L23L3.trace();
Eigen::MatrixXcd L1H0L12=L1*H0*L12; std::complex<double> L1H0L12t=L1H0L12.trace();
Eigen::MatrixXcd L1H0H123=L1*H0*H123; std::complex<double> L1H0H123t=L1H0H123.trace();
Eigen::MatrixXcd L1H0L13=L1*H0*L13; std::complex<double> L1H0L13t=L1H0L13.trace();
Eigen::MatrixXcd L1H0L2=L1*H0*L2; std::complex<double> L1H0L2t=L1H0L2.trace();
Eigen::MatrixXcd L1H0L23=L1*H0*L23; std::complex<double> L1H0L23t=L1H0L23.trace();
Eigen::MatrixXcd L1H0L3=L1*H0*L3; std::complex<double> L1H0L3t=L1H0L3.trace();
Eigen::MatrixXcd L1L12L2=L1*L12*L2; std::complex<double> L1L12L2t=L1L12L2.trace();
Eigen::MatrixXcd L1L12L23=L1*L12*L23; std::complex<double> L1L12L23t=L1L12L23.trace();
Eigen::MatrixXcd L1L12L3=L1*L12*L3; std::complex<double> L1L12L3t=L1L12L3.trace();
Eigen::MatrixXcd L1H123L2=L1*H123*L2; std::complex<double> L1H123L2t=L1H123L2.trace();
Eigen::MatrixXcd L1H123L23=L1*H123*L23; std::complex<double> L1H123L23t=L1H123L23.trace();
Eigen::MatrixXcd L1H123L3=L1*H123*L3; std::complex<double> L1H123L3t=L1H123L3.trace();
Eigen::MatrixXcd L1L13L2=L1*L13*L2; std::complex<double> L1L13L2t=L1L13L2.trace();
Eigen::MatrixXcd L1L13L23=L1*L13*L23; std::complex<double> L1L13L23t=L1L13L23.trace();
Eigen::MatrixXcd L1L13L3=L1*L13*L3; std::complex<double> L1L13L3t=L1L13L3.trace();
Eigen::MatrixXcd L1L2L3=L1*L2*L3; std::complex<double> L1L2L3t=L1L2L3.trace();
Eigen::MatrixXcd L12H0L1=L12*H0*L1; std::complex<double> L12H0L1t=L12H0L1.trace();
Eigen::MatrixXcd L12H0H123=L12*H0*H123; std::complex<double> L12H0H123t=L12H0H123.trace();
Eigen::MatrixXcd L12H0L13=L12*H0*L13; std::complex<double> L12H0L13t=L12H0L13.trace();
Eigen::MatrixXcd L12H0L2=L12*H0*L2; std::complex<double> L12H0L2t=L12H0L2.trace();
Eigen::MatrixXcd L12H0L23=L12*H0*L23; std::complex<double> L12H0L23t=L12H0L23.trace();
Eigen::MatrixXcd L12H0L3=L12*H0*L3; std::complex<double> L12H0L3t=L12H0L3.trace();
Eigen::MatrixXcd L12L1H0=L12*L1*H0; std::complex<double> L12L1H0t=L12L1H0.trace();
Eigen::MatrixXcd L12L1H123=L12*L1*H123; std::complex<double> L12L1H123t=L12L1H123.trace();
Eigen::MatrixXcd L12L1L13=L12*L1*L13; std::complex<double> L12L1L13t=L12L1L13.trace();
Eigen::MatrixXcd L12L1L2=L12*L1*L2; std::complex<double> L12L1L2t=L12L1L2.trace();
Eigen::MatrixXcd L12L1L23=L12*L1*L23; std::complex<double> L12L1L23t=L12L1L23.trace();
Eigen::MatrixXcd L12L1L3=L12*L1*L3; std::complex<double> L12L1L3t=L12L1L3.trace();
Eigen::MatrixXcd L12H123L13=L12*H123*L13; std::complex<double> L12H123L13t=L12H123L13.trace();
Eigen::MatrixXcd L12H123L2=L12*H123*L2; std::complex<double> L12H123L2t=L12H123L2.trace();
Eigen::MatrixXcd L12H123L23=L12*H123*L23; std::complex<double> L12H123L23t=L12H123L23.trace();
Eigen::MatrixXcd L12L13L2=L12*L13*L2; std::complex<double> L12L13L2t=L12L13L2.trace();
Eigen::MatrixXcd L12L13L23=L12*L13*L23; std::complex<double> L12L13L23t=L12L13L23.trace();
Eigen::MatrixXcd L12L2L23=L12*L2*L23; std::complex<double> L12L2L23t=L12L2L23.trace();
Eigen::MatrixXcd L12L2L3=L12*L2*L3; std::complex<double> L12L2L3t=L12L2L3.trace();
Eigen::MatrixXcd H123H0L1=H123*H0*L1; std::complex<double> H123H0L1t=H123H0L1.trace();
Eigen::MatrixXcd H123H0L12=H123*H0*L12; std::complex<double> H123H0L12t=H123H0L12.trace();
Eigen::MatrixXcd H123H0L13=H123*H0*L13; std::complex<double> H123H0L13t=H123H0L13.trace();
Eigen::MatrixXcd H123H0L2=H123*H0*L2; std::complex<double> H123H0L2t=H123H0L2.trace();
Eigen::MatrixXcd H123H0L23=H123*H0*L23; std::complex<double> H123H0L23t=H123H0L23.trace();
Eigen::MatrixXcd H123H0L3=H123*H0*L3; std::complex<double> H123H0L3t=H123H0L3.trace();
Eigen::MatrixXcd H123L1L12=H123*L1*L12; std::complex<double> H123L1L12t=H123L1L12.trace();
Eigen::MatrixXcd H123L1L13=H123*L1*L13; std::complex<double> H123L1L13t=H123L1L13.trace();
Eigen::MatrixXcd H123L1L2=H123*L1*L2; std::complex<double> H123L1L2t=H123L1L2.trace();
Eigen::MatrixXcd H123L1L23=H123*L1*L23; std::complex<double> H123L1L23t=H123L1L23.trace();
Eigen::MatrixXcd H123L1L3=H123*L1*L3; std::complex<double> H123L1L3t=H123L1L3.trace();
Eigen::MatrixXcd H123L12H0=H123*L12*H0; std::complex<double> H123L12H0t=H123L12H0.trace();
Eigen::MatrixXcd H123L12L13=H123*L12*L13; std::complex<double> H123L12L13t=H123L12L13.trace();
Eigen::MatrixXcd H123L12L2=H123*L12*L2; std::complex<double> H123L12L2t=H123L12L2.trace();
Eigen::MatrixXcd H123L12L23=H123*L12*L23; std::complex<double> H123L12L23t=H123L12L23.trace();
Eigen::MatrixXcd H123L12L3=H123*L12*L3; std::complex<double> H123L12L3t=H123L12L3.trace();
Eigen::MatrixXcd H123L13L2=H123*L13*L2; std::complex<double> H123L13L2t=H123L13L2.trace();
Eigen::MatrixXcd H123L13L23=H123*L13*L23; std::complex<double> H123L13L23t=H123L13L23.trace();
Eigen::MatrixXcd H123L13L3=H123*L13*L3; std::complex<double> H123L13L3t=H123L13L3.trace();
Eigen::MatrixXcd H123L2L23=H123*L2*L23; std::complex<double> H123L2L23t=H123L2L23.trace();
Eigen::MatrixXcd H123L23L3=H123*L23*L3; std::complex<double> H123L23L3t=H123L23L3.trace();
Eigen::MatrixXcd L13H0L1=L13*H0*L1; std::complex<double> L13H0L1t=L13H0L1.trace();
Eigen::MatrixXcd L13H0L12=L13*H0*L12; std::complex<double> L13H0L12t=L13H0L12.trace();
Eigen::MatrixXcd L13H0H123=L13*H0*H123; std::complex<double> L13H0H123t=L13H0H123.trace();
Eigen::MatrixXcd L13H0L2=L13*H0*L2; std::complex<double> L13H0L2t=L13H0L2.trace();
Eigen::MatrixXcd L13H0L23=L13*H0*L23; std::complex<double> L13H0L23t=L13H0L23.trace();
Eigen::MatrixXcd L13H0L3=L13*H0*L3; std::complex<double> L13H0L3t=L13H0L3.trace();
Eigen::MatrixXcd L13L1H0=L13*L1*H0; std::complex<double> L13L1H0t=L13L1H0.trace();
Eigen::MatrixXcd L13L1L12=L13*L1*L12; std::complex<double> L13L1L12t=L13L1L12.trace();
Eigen::MatrixXcd L13L1H123=L13*L1*H123; std::complex<double> L13L1H123t=L13L1H123.trace();
Eigen::MatrixXcd L13L1L2=L13*L1*L2; std::complex<double> L13L1L2t=L13L1L2.trace();
Eigen::MatrixXcd L13L1L23=L13*L1*L23; std::complex<double> L13L1L23t=L13L1L23.trace();
Eigen::MatrixXcd L13L1L3=L13*L1*L3; std::complex<double> L13L1L3t=L13L1L3.trace();
Eigen::MatrixXcd L13L12H123=L13*L12*H123; std::complex<double> L13L12H123t=L13L12H123.trace();
Eigen::MatrixXcd L13L12L2=L13*L12*L2; std::complex<double> L13L12L2t=L13L12L2.trace();
Eigen::MatrixXcd L13L12L23=L13*L12*L23; std::complex<double> L13L12L23t=L13L12L23.trace();
Eigen::MatrixXcd L13L12L3=L13*L12*L3; std::complex<double> L13L12L3t=L13L12L3.trace();
Eigen::MatrixXcd L13H123H0=L13*H123*H0; std::complex<double> L13H123H0t=L13H123H0.trace();
Eigen::MatrixXcd L13H123L12=L13*H123*L12; std::complex<double> L13H123L12t=L13H123L12.trace();
Eigen::MatrixXcd L13H123L2=L13*H123*L2; std::complex<double> L13H123L2t=L13H123L2.trace();
Eigen::MatrixXcd L13H123L23=L13*H123*L23; std::complex<double> L13H123L23t=L13H123L23.trace();
Eigen::MatrixXcd L13H123L3=L13*H123*L3; std::complex<double> L13H123L3t=L13H123L3.trace();
Eigen::MatrixXcd L13L2L3=L13*L2*L3; std::complex<double> L13L2L3t=L13L2L3.trace();
Eigen::MatrixXcd L13L23L3=L13*L23*L3; std::complex<double> L13L23L3t=L13L23L3.trace();
Eigen::MatrixXcd L2H0L1=L2*H0*L1; std::complex<double> L2H0L1t=L2H0L1.trace();
Eigen::MatrixXcd L2H0L12=L2*H0*L12; std::complex<double> L2H0L12t=L2H0L12.trace();
Eigen::MatrixXcd L2H0H123=L2*H0*H123; std::complex<double> L2H0H123t=L2H0H123.trace();
Eigen::MatrixXcd L2H0L13=L2*H0*L13; std::complex<double> L2H0L13t=L2H0L13.trace();
Eigen::MatrixXcd L2H0L23=L2*H0*L23; std::complex<double> L2H0L23t=L2H0L23.trace();
Eigen::MatrixXcd L2H0L3=L2*H0*L3; std::complex<double> L2H0L3t=L2H0L3.trace();
Eigen::MatrixXcd L2L1H0=L2*L1*H0; std::complex<double> L2L1H0t=L2L1H0.trace();
Eigen::MatrixXcd L2L1L12=L2*L1*L12; std::complex<double> L2L1L12t=L2L1L12.trace();
Eigen::MatrixXcd L2L1H123=L2*L1*H123; std::complex<double> L2L1H123t=L2L1H123.trace();
Eigen::MatrixXcd L2L1L13=L2*L1*L13; std::complex<double> L2L1L13t=L2L1L13.trace();
Eigen::MatrixXcd L2L1L23=L2*L1*L23; std::complex<double> L2L1L23t=L2L1L23.trace();
Eigen::MatrixXcd L2L1L3=L2*L1*L3; std::complex<double> L2L1L3t=L2L1L3.trace();
Eigen::MatrixXcd L2L12H0=L2*L12*H0; std::complex<double> L2L12H0t=L2L12H0.trace();
Eigen::MatrixXcd L2L12L1=L2*L12*L1; std::complex<double> L2L12L1t=L2L12L1.trace();
Eigen::MatrixXcd L2L12H123=L2*L12*H123; std::complex<double> L2L12H123t=L2L12H123.trace();
Eigen::MatrixXcd L2L12L13=L2*L12*L13; std::complex<double> L2L12L13t=L2L12L13.trace();
Eigen::MatrixXcd L2L12L23=L2*L12*L23; std::complex<double> L2L12L23t=L2L12L23.trace();
Eigen::MatrixXcd L2L12L3=L2*L12*L3; std::complex<double> L2L12L3t=L2L12L3.trace();
Eigen::MatrixXcd L2H123L1=L2*H123*L1; std::complex<double> L2H123L1t=L2H123L1.trace();
Eigen::MatrixXcd L2H123L12=L2*H123*L12; std::complex<double> L2H123L12t=L2H123L12.trace();
Eigen::MatrixXcd L2H123L13=L2*H123*L13; std::complex<double> L2H123L13t=L2H123L13.trace();
Eigen::MatrixXcd L2H123L23=L2*H123*L23; std::complex<double> L2H123L23t=L2H123L23.trace();
Eigen::MatrixXcd L2H123L3=L2*H123*L3; std::complex<double> L2H123L3t=L2H123L3.trace();
Eigen::MatrixXcd L2L13H0=L2*L13*H0; std::complex<double> L2L13H0t=L2L13H0.trace();
Eigen::MatrixXcd L2L13L1=L2*L13*L1; std::complex<double> L2L13L1t=L2L13L1.trace();
Eigen::MatrixXcd L2L13L12=L2*L13*L12; std::complex<double> L2L13L12t=L2L13L12.trace();
Eigen::MatrixXcd L2L13H123=L2*L13*H123; std::complex<double> L2L13H123t=L2L13H123.trace();
Eigen::MatrixXcd L2L13L23=L2*L13*L23; std::complex<double> L2L13L23t=L2L13L23.trace();
Eigen::MatrixXcd L2L13L3=L2*L13*L3; std::complex<double> L2L13L3t=L2L13L3.trace();
Eigen::MatrixXcd L2L23L3=L2*L23*L3; std::complex<double> L2L23L3t=L2L23L3.trace();
Eigen::MatrixXcd L23H0L1=L23*H0*L1; std::complex<double> L23H0L1t=L23H0L1.trace();
Eigen::MatrixXcd L23H0L12=L23*H0*L12; std::complex<double> L23H0L12t=L23H0L12.trace();
Eigen::MatrixXcd L23H0H123=L23*H0*H123; std::complex<double> L23H0H123t=L23H0H123.trace();
Eigen::MatrixXcd L23H0L13=L23*H0*L13; std::complex<double> L23H0L13t=L23H0L13.trace();
Eigen::MatrixXcd L23H0L2=L23*H0*L2; std::complex<double> L23H0L2t=L23H0L2.trace();
Eigen::MatrixXcd L23H0L3=L23*H0*L3; std::complex<double> L23H0L3t=L23H0L3.trace();
Eigen::MatrixXcd L23L1H0=L23*L1*H0; std::complex<double> L23L1H0t=L23L1H0.trace();
Eigen::MatrixXcd L23L1L12=L23*L1*L12; std::complex<double> L23L1L12t=L23L1L12.trace();
Eigen::MatrixXcd L23L1H123=L23*L1*H123; std::complex<double> L23L1H123t=L23L1H123.trace();
Eigen::MatrixXcd L23L1L13=L23*L1*L13; std::complex<double> L23L1L13t=L23L1L13.trace();
Eigen::MatrixXcd L23L1L2=L23*L1*L2; std::complex<double> L23L1L2t=L23L1L2.trace();
Eigen::MatrixXcd L23L1L3=L23*L1*L3; std::complex<double> L23L1L3t=L23L1L3.trace();
Eigen::MatrixXcd L23L12L1=L23*L12*L1; std::complex<double> L23L12L1t=L23L12L1.trace();
Eigen::MatrixXcd L23L12H123=L23*L12*H123; std::complex<double> L23L12H123t=L23L12H123.trace();
Eigen::MatrixXcd L23L12L13=L23*L12*L13; std::complex<double> L23L12L13t=L23L12L13.trace();
Eigen::MatrixXcd L23L12L2=L23*L12*L2; std::complex<double> L23L12L2t=L23L12L2.trace();
Eigen::MatrixXcd L23L12L3=L23*L12*L3; std::complex<double> L23L12L3t=L23L12L3.trace();
Eigen::MatrixXcd L23H123H0=L23*H123*H0; std::complex<double> L23H123H0t=L23H123H0.trace();
Eigen::MatrixXcd L23H123L1=L23*H123*L1; std::complex<double> L23H123L1t=L23H123L1.trace();
Eigen::MatrixXcd L23H123L12=L23*H123*L12; std::complex<double> L23H123L12t=L23H123L12.trace();
Eigen::MatrixXcd L23H123L13=L23*H123*L13; std::complex<double> L23H123L13t=L23H123L13.trace();
Eigen::MatrixXcd L23H123L2=L23*H123*L2; std::complex<double> L23H123L2t=L23H123L2.trace();
Eigen::MatrixXcd L23H123L3=L23*H123*L3; std::complex<double> L23H123L3t=L23H123L3.trace();
Eigen::MatrixXcd L23L13L1=L23*L13*L1; std::complex<double> L23L13L1t=L23L13L1.trace();
Eigen::MatrixXcd L23L13L12=L23*L13*L12; std::complex<double> L23L13L12t=L23L13L12.trace();
Eigen::MatrixXcd L23L13H123=L23*L13*H123; std::complex<double> L23L13H123t=L23L13H123.trace();
Eigen::MatrixXcd L23L13L2=L23*L13*L2; std::complex<double> L23L13L2t=L23L13L2.trace();
Eigen::MatrixXcd L23L13L3=L23*L13*L3; std::complex<double> L23L13L3t=L23L13L3.trace();
Eigen::MatrixXcd L23L2H0=L23*L2*H0; std::complex<double> L23L2H0t=L23L2H0.trace();
Eigen::MatrixXcd L23L2L12=L23*L2*L12; std::complex<double> L23L2L12t=L23L2L12.trace();
Eigen::MatrixXcd L23L2H123=L23*L2*H123; std::complex<double> L23L2H123t=L23L2H123.trace();
Eigen::MatrixXcd L23L2L3=L23*L2*L3; std::complex<double> L23L2L3t=L23L2L3.trace();
Eigen::MatrixXcd L3H0L1=L3*H0*L1; std::complex<double> L3H0L1t=L3H0L1.trace();
Eigen::MatrixXcd L3H0L12=L3*H0*L12; std::complex<double> L3H0L12t=L3H0L12.trace();
Eigen::MatrixXcd L3H0H123=L3*H0*H123; std::complex<double> L3H0H123t=L3H0H123.trace();
Eigen::MatrixXcd L3H0L13=L3*H0*L13; std::complex<double> L3H0L13t=L3H0L13.trace();
Eigen::MatrixXcd L3H0L2=L3*H0*L2; std::complex<double> L3H0L2t=L3H0L2.trace();
Eigen::MatrixXcd L3H0L23=L3*H0*L23; std::complex<double> L3H0L23t=L3H0L23.trace();
Eigen::MatrixXcd L3L1H0=L3*L1*H0; std::complex<double> L3L1H0t=L3L1H0.trace();
Eigen::MatrixXcd L3L1L12=L3*L1*L12; std::complex<double> L3L1L12t=L3L1L12.trace();
Eigen::MatrixXcd L3L1H123=L3*L1*H123; std::complex<double> L3L1H123t=L3L1H123.trace();
Eigen::MatrixXcd L3L1L13=L3*L1*L13; std::complex<double> L3L1L13t=L3L1L13.trace();
Eigen::MatrixXcd L3L1L2=L3*L1*L2; std::complex<double> L3L1L2t=L3L1L2.trace();
Eigen::MatrixXcd L3L1L23=L3*L1*L23; std::complex<double> L3L1L23t=L3L1L23.trace();
Eigen::MatrixXcd L3L12H0=L3*L12*H0; std::complex<double> L3L12H0t=L3L12H0.trace();
Eigen::MatrixXcd L3L12L1=L3*L12*L1; std::complex<double> L3L12L1t=L3L12L1.trace();
Eigen::MatrixXcd L3L12H123=L3*L12*H123; std::complex<double> L3L12H123t=L3L12H123.trace();
Eigen::MatrixXcd L3L12L13=L3*L12*L13; std::complex<double> L3L12L13t=L3L12L13.trace();
Eigen::MatrixXcd L3L12L2=L3*L12*L2; std::complex<double> L3L12L2t=L3L12L2.trace();
Eigen::MatrixXcd L3L12L23=L3*L12*L23; std::complex<double> L3L12L23t=L3L12L23.trace();
Eigen::MatrixXcd L3H123L1=L3*H123*L1; std::complex<double> L3H123L1t=L3H123L1.trace();
Eigen::MatrixXcd L3H123L13=L3*H123*L13; std::complex<double> L3H123L13t=L3H123L13.trace();
Eigen::MatrixXcd L3H123L2=L3*H123*L2; std::complex<double> L3H123L2t=L3H123L2.trace();
Eigen::MatrixXcd L3H123L23=L3*H123*L23; std::complex<double> L3H123L23t=L3H123L23.trace();
Eigen::MatrixXcd L3L13H0=L3*L13*H0; std::complex<double> L3L13H0t=L3L13H0.trace();
Eigen::MatrixXcd L3L13L1=L3*L13*L1; std::complex<double> L3L13L1t=L3L13L1.trace();
Eigen::MatrixXcd L3L13H123=L3*L13*H123; std::complex<double> L3L13H123t=L3L13H123.trace();
Eigen::MatrixXcd L3L13L2=L3*L13*L2; std::complex<double> L3L13L2t=L3L13L2.trace();
Eigen::MatrixXcd L3L13L23=L3*L13*L23; std::complex<double> L3L13L23t=L3L13L23.trace();
Eigen::MatrixXcd L3L2H0=L3*L2*H0; std::complex<double> L3L2H0t=L3L2H0.trace();
Eigen::MatrixXcd L3L2L1=L3*L2*L1; std::complex<double> L3L2L1t=L3L2L1.trace();
Eigen::MatrixXcd L3L2L12=L3*L2*L12; std::complex<double> L3L2L12t=L3L2L12.trace();
Eigen::MatrixXcd L3L2L13=L3*L2*L13; std::complex<double> L3L2L13t=L3L2L13.trace();
Eigen::MatrixXcd L3L2L23=L3*L2*L23; std::complex<double> L3L2L23t=L3L2L23.trace();
Eigen::MatrixXcd L3L23H0=L3*L23*H0; std::complex<double> L3L23H0t=L3L23H0.trace();
Eigen::MatrixXcd L3L23H123=L3*L23*H123; std::complex<double> L3L23H123t=L3L23H123.trace();
Eigen::MatrixXcd L3L23L13=L3*L23*L13; std::complex<double> L3L23L13t=L3L23L13.trace();
Eigen::MatrixXcd L3L23L2=L3*L23*L2; std::complex<double> L3L23L2t=L3L23L2.trace();
Eigen::MatrixXcd H0L2L12L1=H0*L2*L12*L1; std::complex<double> H0L2L12L1t=H0L2L12L1.trace();
Eigen::MatrixXcd H0L2L13H123=H0*L2*L13*H123; std::complex<double> H0L2L13H123t=H0L2L13H123.trace();
Eigen::MatrixXcd H0L23H123L1=H0*L23*H123*L1; std::complex<double> H0L23H123L1t=H0L23H123L1.trace();
Eigen::MatrixXcd H0L23L13L12=H0*L23*L13*L12; std::complex<double> H0L23L13L12t=H0L23L13L12.trace();
Eigen::MatrixXcd H0L3H123L12=H0*L3*H123*L12; std::complex<double> H0L3H123L12t=H0L3H123L12.trace();
Eigen::MatrixXcd H0L3L13L1=H0*L3*L13*L1; std::complex<double> H0L3L13L1t=H0L3L13L1.trace();
Eigen::MatrixXcd H0L3L23L2=H0*L3*L23*L2; std::complex<double> H0L3L23L2t=H0L3L23L2.trace();
Eigen::MatrixXcd L1H0L2L12=L1*H0*L2*L12; std::complex<double> L1H0L2L12t=L1H0L2L12.trace();
Eigen::MatrixXcd L1H0L23H123=L1*H0*L23*H123; std::complex<double> L1H0L23H123t=L1H0L23H123.trace();
Eigen::MatrixXcd L1H0L3L13=L1*H0*L3*L13; std::complex<double> L1H0L3L13t=L1H0L3L13.trace();
Eigen::MatrixXcd L1L12L2H0=L1*L12*L2*H0; std::complex<double> L1L12L2H0t=L1L12L2H0.trace();
Eigen::MatrixXcd L1H123L23H0=L1*H123*L23*H0; std::complex<double> L1H123L23H0t=L1H123L23H0.trace();
Eigen::MatrixXcd L1L13H123L12=L1*L13*H123*L12; std::complex<double> L1L13H123L12t=L1L13H123L12.trace();
Eigen::MatrixXcd L1L13L3H0=L1*L13*L3*H0; std::complex<double> L1L13L3H0t=L1L13L3H0.trace();
Eigen::MatrixXcd L1L2L12H0=L1*L2*L12*H0; std::complex<double> L1L2L12H0t=L1L2L12H0.trace();
Eigen::MatrixXcd L1L23H123H0=L1*L23*H123*H0; std::complex<double> L1L23H123H0t=L1L23H123H0.trace();
Eigen::MatrixXcd L1L23L2L13=L1*L23*L2*L13; std::complex<double> L1L23L2L13t=L1L23L2L13.trace();
Eigen::MatrixXcd L1L3L13H0=L1*L3*L13*H0; std::complex<double> L1L3L13H0t=L1L3L13H0.trace();
Eigen::MatrixXcd L1L3L2H123=L1*L3*L2*H123; std::complex<double> L1L3L2H123t=L1L3L2H123.trace();
Eigen::MatrixXcd L1L3L23L12=L1*L3*L23*L12; std::complex<double> L1L3L23L12t=L1L3L23L12.trace();
Eigen::MatrixXcd L12H0L2L1=L12*H0*L2*L1; std::complex<double> L12H0L2L1t=L12H0L2L1.trace();
Eigen::MatrixXcd L12H0L23L13=L12*H0*L23*L13; std::complex<double> L12H0L23L13t=L12H0L23L13.trace();
Eigen::MatrixXcd L12H0L3H123=L12*H0*L3*H123; std::complex<double> L12H0L3H123t=L12H0L3H123.trace();
Eigen::MatrixXcd L12L1H0L2=L12*L1*H0*L2; std::complex<double> L12L1H0L2t=L12L1H0L2.trace();
Eigen::MatrixXcd L12L1L13H123=L12*L1*L13*H123; std::complex<double> L12L1L13H123t=L12L1L13H123.trace();
Eigen::MatrixXcd L12L1L2H0=L12*L1*L2*H0; std::complex<double> L12L1L2H0t=L12L1L2H0.trace();
Eigen::MatrixXcd L12L1L3L23=L12*L1*L3*L23; std::complex<double> L12L1L3L23t=L12L1L3L23.trace();
Eigen::MatrixXcd L12H123L13L1=L12*H123*L13*L1; std::complex<double> L12H123L13L1t=L12H123L13L1.trace();
Eigen::MatrixXcd L12H123L3H0=L12*H123*L3*H0; std::complex<double> L12H123L3H0t=L12H123L3H0.trace();
Eigen::MatrixXcd L12L13H123L1=L12*L13*H123*L1; std::complex<double> L12L13H123L1t=L12L13H123L1.trace();
Eigen::MatrixXcd L12L13L23H0=L12*L13*L23*H0; std::complex<double> L12L13L23H0t=L12L13L23H0.trace();
Eigen::MatrixXcd L12L23L13H0=L12*L23*L13*H0; std::complex<double> L12L23L13H0t=L12L23L13H0.trace();
Eigen::MatrixXcd L12L23L2H123=L12*L23*L2*H123; std::complex<double> L12L23L2H123t=L12L23L2H123.trace();
Eigen::MatrixXcd L12L23L3L1=L12*L23*L3*L1; std::complex<double> L12L23L3L1t=L12L23L3L1.trace();
Eigen::MatrixXcd L12L3H123H0=L12*L3*H123*H0; std::complex<double> L12L3H123H0t=L12L3H123H0.trace();
Eigen::MatrixXcd L12L3L2L13=L12*L3*L2*L13; std::complex<double> L12L3L2L13t=L12L3L2L13.trace();
Eigen::MatrixXcd L12L3L23L1=L12*L3*L23*L1; std::complex<double> L12L3L23L1t=L12L3L23L1.trace();
Eigen::MatrixXcd H123H0L2L13=H123*H0*L2*L13; std::complex<double> H123H0L2L13t=H123H0L2L13.trace();
Eigen::MatrixXcd H123H0L23L1=H123*H0*L23*L1; std::complex<double> H123H0L23L1t=H123H0L23L1.trace();
Eigen::MatrixXcd H123H0L3L12=H123*H0*L3*L12; std::complex<double> H123H0L3L12t=H123H0L3L12.trace();
Eigen::MatrixXcd H123L1H0L23=H123*L1*H0*L23; std::complex<double> H123L1H0L23t=H123L1H0L23.trace();
Eigen::MatrixXcd H123L1L13L12=H123*L1*L13*L12; std::complex<double> H123L1L13L12t=H123L1L13L12.trace();
Eigen::MatrixXcd H123L1L23H0=H123*L1*L23*H0; std::complex<double> H123L1L23H0t=H123L1L23H0.trace();
Eigen::MatrixXcd H123L1L3L2=H123*L1*L3*L2; std::complex<double> H123L1L3L2t=H123L1L3L2.trace();
Eigen::MatrixXcd H123L12H0L3=H123*L12*H0*L3; std::complex<double> H123L12H0L3t=H123L12H0L3.trace();
Eigen::MatrixXcd H123L12L1L13=H123*L12*L1*L13; std::complex<double> H123L12L1L13t=H123L12L1L13.trace();
Eigen::MatrixXcd H123L12L13L1=H123*L12*L13*L1; std::complex<double> H123L12L13L1t=H123L12L13L1.trace();
Eigen::MatrixXcd H123L12L23L2=H123*L12*L23*L2; std::complex<double> H123L12L23L2t=H123L12L23L2.trace();
Eigen::MatrixXcd H123L12L3H0=H123*L12*L3*H0; std::complex<double> H123L12L3H0t=H123L12L3H0.trace();
Eigen::MatrixXcd H123L13L2H0=H123*L13*L2*H0; std::complex<double> H123L13L2H0t=H123L13L2H0.trace();
Eigen::MatrixXcd H123L2L13H0=H123*L2*L13*H0; std::complex<double> H123L2L13H0t=H123L2L13H0.trace();
Eigen::MatrixXcd H123L2L23L12=H123*L2*L23*L12; std::complex<double> H123L2L23L12t=H123L2L23L12.trace();
Eigen::MatrixXcd H123L2L3L1=H123*L2*L3*L1; std::complex<double> H123L2L3L1t=H123L2L3L1.trace();
Eigen::MatrixXcd H123L23L2L12=H123*L23*L2*L12; std::complex<double> H123L23L2L12t=H123L23L2L12.trace();
Eigen::MatrixXcd H123L3L2L1=H123*L3*L2*L1; std::complex<double> H123L3L2L1t=H123L3L2L1.trace();
Eigen::MatrixXcd H123L3L23L13=H123*L3*L23*L13; std::complex<double> H123L3L23L13t=H123L3L23L13.trace();
Eigen::MatrixXcd L13H0L2H123=L13*H0*L2*H123; std::complex<double> L13H0L2H123t=L13H0L2H123.trace();
Eigen::MatrixXcd L13H0L23L12=L13*H0*L23*L12; std::complex<double> L13H0L23L12t=L13H0L23L12.trace();
Eigen::MatrixXcd L13H0L3L1=L13*H0*L3*L1; std::complex<double> L13H0L3L1t=L13H0L3L1.trace();
Eigen::MatrixXcd L13L1H0L3=L13*L1*H0*L3; std::complex<double> L13L1H0L3t=L13L1H0L3.trace();
Eigen::MatrixXcd L13L1L23L2=L13*L1*L23*L2; std::complex<double> L13L1L23L2t=L13L1L23L2.trace();
Eigen::MatrixXcd L13L1L3H0=L13*L1*L3*H0; std::complex<double> L13L1L3H0t=L13L1L3H0.trace();
Eigen::MatrixXcd L13L12H0L23=L13*L12*H0*L23; std::complex<double> L13L12H0L23t=L13L12H0L23.trace();
Eigen::MatrixXcd L13L12L1H123=L13*L12*L1*H123; std::complex<double> L13L12L1H123t=L13L12L1H123.trace();
Eigen::MatrixXcd L13L12L23H0=L13*L12*L23*H0; std::complex<double> L13L12L23H0t=L13L12L23H0.trace();
Eigen::MatrixXcd L13L12L3L2=L13*L12*L3*L2; std::complex<double> L13L12L3L2t=L13L12L3L2.trace();
Eigen::MatrixXcd L13H123H0L2=L13*H123*H0*L2; std::complex<double> L13H123H0L2t=L13H123H0L2.trace();
Eigen::MatrixXcd L13H123L2H0=L13*H123*L2*H0; std::complex<double> L13H123L2H0t=L13H123L2H0.trace();
Eigen::MatrixXcd L13H123L3L23=L13*H123*L3*L23; std::complex<double> L13H123L3L23t=L13H123L3L23.trace();
Eigen::MatrixXcd L13L2L23L1=L13*L2*L23*L1; std::complex<double> L13L2L23L1t=L13L2L23L1.trace();
Eigen::MatrixXcd L13L2L3L12=L13*L2*L3*L12; std::complex<double> L13L2L3L12t=L13L2L3L12.trace();
Eigen::MatrixXcd L13L23L2L1=L13*L23*L2*L1; std::complex<double> L13L23L2L1t=L13L23L2L1.trace();
Eigen::MatrixXcd L13L23L3H123=L13*L23*L3*H123; std::complex<double> L13L23L3H123t=L13L23L3H123.trace();
Eigen::MatrixXcd L13L3L2L12=L13*L3*L2*L12; std::complex<double> L13L3L2L12t=L13L3L2L12.trace();
Eigen::MatrixXcd L13L3L23H123=L13*L3*L23*H123; std::complex<double> L13L3L23H123t=L13L3L23H123.trace();
Eigen::MatrixXcd L2H0L3L23=L2*H0*L3*L23; std::complex<double> L2H0L3L23t=L2H0L3L23.trace();
Eigen::MatrixXcd L2L1H0L12=L2*L1*H0*L12; std::complex<double> L2L1H0L12t=L2L1H0L12.trace();
Eigen::MatrixXcd L2L1L23L13=L2*L1*L23*L13; std::complex<double> L2L1L23L13t=L2L1L23L13.trace();
Eigen::MatrixXcd L2L1L3H123=L2*L1*L3*H123; std::complex<double> L2L1L3H123t=L2L1L3H123.trace();
Eigen::MatrixXcd L2L12L23H123=L2*L12*L23*H123; std::complex<double> L2L12L23H123t=L2L12L23H123.trace();
Eigen::MatrixXcd L2L12L3L13=L2*L12*L3*L13; std::complex<double> L2L12L3L13t=L2L12L3L13.trace();
Eigen::MatrixXcd L2H123H0L13=L2*H123*H0*L13; std::complex<double> L2H123H0L13t=L2H123H0L13.trace();
Eigen::MatrixXcd L2H123L1L3=L2*H123*L1*L3; std::complex<double> L2H123L1L3t=L2H123L1L3.trace();
Eigen::MatrixXcd L2H123L12L23=L2*H123*L12*L23; std::complex<double> L2H123L12L23t=L2H123L12L23.trace();
Eigen::MatrixXcd L2H123L23L12=L2*H123*L23*L12; std::complex<double> L2H123L23L12t=L2H123L23L12.trace();
Eigen::MatrixXcd L2H123L3L1=L2*H123*L3*L1; std::complex<double> L2H123L3L1t=L2H123L3L1.trace();
Eigen::MatrixXcd L2L13L1L23=L2*L13*L1*L23; std::complex<double> L2L13L1L23t=L2L13L1L23.trace();
Eigen::MatrixXcd L2L13L12L3=L2*L13*L12*L3; std::complex<double> L2L13L12L3t=L2L13L12L3.trace();
Eigen::MatrixXcd L2L13L23L1=L2*L13*L23*L1; std::complex<double> L2L13L23L1t=L2L13L23L1.trace();
Eigen::MatrixXcd L2L13L3L12=L2*L13*L3*L12; std::complex<double> L2L13L3L12t=L2L13L3L12.trace();
Eigen::MatrixXcd L2L23L3H0=L2*L23*L3*H0; std::complex<double> L2L23L3H0t=L2L23L3H0.trace();
Eigen::MatrixXcd L2L3L23H0=L2*L3*L23*H0; std::complex<double> L2L3L23H0t=L2L3L23H0.trace();
Eigen::MatrixXcd L23H0L3L2=L23*H0*L3*L2; std::complex<double> L23H0L3L2t=L23H0L3L2.trace();
Eigen::MatrixXcd L23L1H0H123=L23*L1*H0*H123; std::complex<double> L23L1H0H123t=L23L1H0H123.trace();
Eigen::MatrixXcd L23L1L3L12=L23*L1*L3*L12; std::complex<double> L23L1L3L12t=L23L1L3L12.trace();
Eigen::MatrixXcd L23L12H0L13=L23*L12*H0*L13; std::complex<double> L23L12H0L13t=L23L12H0L13.trace();
Eigen::MatrixXcd L23L12L1L3=L23*L12*L1*L3; std::complex<double> L23L12L1L3t=L23L12L1L3.trace();
Eigen::MatrixXcd L23L12L3L1=L23*L12*L3*L1; std::complex<double> L23L12L3L1t=L23L12L3L1.trace();
Eigen::MatrixXcd L23H123L12L2=L23*H123*L12*L2; std::complex<double> L23H123L12L2t=L23H123L12L2.trace();
Eigen::MatrixXcd L23H123L3L13=L23*H123*L3*L13; std::complex<double> L23H123L3L13t=L23H123L3L13.trace();
Eigen::MatrixXcd L23L13L1L2=L23*L13*L1*L2; std::complex<double> L23L13L1L2t=L23L13L1L2.trace();
Eigen::MatrixXcd L23L13H123L3=L23*L13*H123*L3; std::complex<double> L23L13H123L3t=L23L13H123L3.trace();
Eigen::MatrixXcd L23L13L3H123=L23*L13*L3*H123; std::complex<double> L23L13L3H123t=L23L13L3H123.trace();
Eigen::MatrixXcd L23L2H0L3=L23*L2*H0*L3; std::complex<double> L23L2H0L3t=L23L2H0L3.trace();
Eigen::MatrixXcd L23L2L3H0=L23*L2*L3*H0; std::complex<double> L23L2L3H0t=L23L2L3H0.trace();
Eigen::MatrixXcd L3L1H0L13=L3*L1*H0*L13; std::complex<double> L3L1H0L13t=L3L1H0L13.trace();
Eigen::MatrixXcd L3L12H0H123=L3*L12*H0*H123; std::complex<double> L3L12H0H123t=L3L12H0H123.trace();
Eigen::MatrixXcd L3L12L1L23=L3*L12*L1*L23; std::complex<double> L3L12L1L23t=L3L12L1L23.trace();
Eigen::MatrixXcd L3H123L1L2=L3*H123*L1*L2; std::complex<double> L3H123L1L2t=L3H123L1L2.trace();
Eigen::MatrixXcd L3L13L12L2=L3*L13*L12*L2; std::complex<double> L3L13L12L2t=L3L13L12L2.trace();
Eigen::MatrixXcd L3L13H123L23=L3*L13*H123*L23; std::complex<double> L3L13H123L23t=L3L13H123L23.trace();
Eigen::MatrixXcd L3L2H0L23=L3*L2*H0*L23; std::complex<double> L3L2H0L23t=L3L2H0L23.trace();

// after we defined all of these Mathematcia kindly also provides the rest
S4=24.*H0H0t*H0H0t+32.*H0t*H0H0H0t+8.*size*H0H0H0H0t-16.*H0H0t*L1L1t+24.*L1L1t*L1L1t-32.*L1t*L1L1L1t+8.*size*L1L1L1L1t-48.*H0H0t*L12L12t+48.*L1L1t*L12L12t+24.*L12L12t*L12L12t-32.*L12t*L12L12L12t+8.*size*L12L12L12L12t+16.*H0H0t*H123H123t-48.*L1L1t*H123H123t-16.*L12L12t*H123H123t+24.*H123H123t*H123H123t+32.*H123t*H123H123H123t+8.*size*H123H123H123H123t-48.*H0H0t*L13L13t+48.*L1L1t*L13L13t+16.*L12L12t*L13L13t-16.*H123H123t*L13L13t+24.*L13L13t*L13L13t-32.*L13t*L13L13L13t+8.*size*L13L13L13L13t-16.*H0H0t*L2L2t+16.*L1L1t*L2L2t+48.*L12L12t*L2L2t-48.*H123H123t*L2L2t+16.*L13L13t*L2L2t+24.*L2L2t*L2L2t-32.*L2t*L2L2L2t+8.*size*L2L2L2L2t-48.*H0H0t*L23L23t+16.*L1L1t*L23L23t+16.*L12L12t*L23L23t-16.*H123H123t*L23L23t+16.*L13L13t*L23L23t+48.*L2L2t*L23L23t+24.*L23L23t*L23L23t-32.*L23t*L23L23L23t+8.*size*L23L23L23L23t-16.*H0H0t*L3L3t+16.*L1L1t*L3L3t+16.*L12L12t*L3L3t-48.*H123H123t*L3L3t+48.*L13L13t*L3L3t+16.*L2L2t*L3L3t+48.*L23L23t*L3L3t+24.*L3L3t*L3L3t-32.*L3t*L3L3L3t+8.*size*L3L3L3L3t+32.*H0L1t*H0L1t+16.*size*H0L1H0L1t-32.*H0t*H0L1L1t+96.*H0L12t*H0L12t-16.*size*H0L12H0L12t-96.*H0t*H0L12L12t+32.*H0H123t*H0H123t-16.*size*H0H123H0H123t+32.*H0t*H0H123H123t+96.*H0L13t*H0L13t-16.*size*H0L13H0L13t-96.*H0t*H0L13L13t+32.*H0L2t*H0L2t+16.*size*H0L2H0L2t-32.*H0t*H0L2L2t+96.*H0L23t*H0L23t-16.*size*H0L23H0L23t-96.*H0t*H0L23L23t+32.*H0L3t*H0L3t+16.*size*H0L3H0L3t-32.*H0t*H0L3L3t+32.*L1t*H0H0L1t-32.*size*H0H0L1L1t+96.*L12t*H0H0L12t-32.*size*H0H0L12L12t+32.*H123t*H0H0H123t+32.*size*H0H0H123H123t+96.*L13t*H0H0L13t-32.*size*H0H0L13L13t+32.*L2t*H0H0L2t-32.*size*H0H0L2L2t+96.*L23t*H0H0L23t-32.*size*H0H0L23L23t+32.*L3t*H0H0L3t-32.*size*H0H0L3L3t+96.*L1L12t*L1L12t+16.*size*L1L12L1L12t-96.*L1t*L1L12L12t-64.*H0L23t*L1H123t+96.*L1H123t*L1H123t-16.*size*L1H123L1H123t+96.*L1t*L1H123H123t+96.*L1L13t*L1L13t+16.*size*L1L13L1L13t-96.*L1t*L1L13L13t+32.*L1L2t*L1L2t-16.*size*L1L2L1L2t-32.*L1t*L1L2L2t+64.*H0H123t*L1L23t+32.*L1L23t*L1L23t-16.*size*L1L23L1L23t-32.*L1t*L1L23L23t+32.*L1L3t*L1L3t-16.*size*L1L3L1L3t-32.*L1t*L1L3L3t-96.*L12t*L1L1L12t+32.*size*L1L1L12L12t-96.*H123t*L1L1H123t-32.*size*L1L1H123H123t-96.*L13t*L1L1L13t+32.*size*L1L1L13L13t-32.*L2t*L1L1L2t+32.*size*L1L1L2L2t-32.*L23t*L1L1L23t+32.*size*L1L1L23L23t-32.*L3t*L1L1L3t+32.*size*L1L1L3L3t-64.*H0L3t*L12H123t+32.*L12H123t*L12H123t+16.*size*L12H123L12H123t+32.*L12t*L12H123H123t+32.*L12L13t*L12L13t-16.*size*L12L13L12L13t-32.*L12t*L12L13L13t+96.*L12L2t*L12L2t+16.*size*L12L2L12L2t-96.*L12t*L12L2L2t-64.*L1L3t*L12L23t+32.*L12L23t*L12L23t-16.*size*L12L23L12L23t-32.*L12t*L12L23L23t+64.*H0H123t*L12L3t-64.*L1L23t*L12L3t+32.*L12L3t*L12L3t-16.*size*L12L3L12L3t-32.*L12t*L12L3L3t-32.*H123t*L12L12H123t-32.*size*L12L12H123H123t-32.*L13t*L12L12L13t+32.*size*L12L12L13L13t-96.*L2t*L12L12L2t+32.*size*L12L12L2L2t-32.*L23t*L12L12L23t+32.*size*L12L12L23L23t-32.*L3t*L12L12L3t+32.*size*L12L12L3L3t+64.*H0L2t*H123L13t+32.*H123L13t*H123L13t+16.*size*H123L13H123L13t-32.*H123t*H123L13L13t+64.*H0L13t*H123L2t+96.*H123L2t*H123L2t-16.*size*H123L2H123L2t-96.*H123t*H123L2L2t-64.*H0L1t*H123L23t+32.*H123L23t*H123L23t+16.*size*H123L23H123L23t-32.*H123t*H123L23L23t-64.*H0L12t*H123L3t+96.*H123L3t*H123L3t-16.*size*H123L3H123L3t-96.*H123t*H123L3L3t+32.*L13t*H123H123L13t-32.*size*H123H123L13L13t+96.*L2t*H123H123L2t-32.*size*H123H123L2L2t+32.*L23t*H123H123L23t-32.*size*H123H123L23L23t+96.*L3t*H123H123L3t-32.*size*H123H123L3L3t-64.*H0H123t*L13L2t+64.*L1L23t*L13L2t+64.*L12L3t*L13L2t+32.*L13L2t*L13L2t-16.*size*L13L2L13L2t-32.*L13t*L13L2L2t+64.*L1L2t*L13L23t+32.*L13L23t*L13L23t-16.*size*L13L23L13L23t-32.*L13t*L13L23L23t+64.*L12L2t*L13L3t+96.*L13L3t*L13L3t+16.*size*L13L3L13L3t-96.*L13t*L13L3L3t-32.*L2t*L13L13L2t+32.*size*L13L13L2L2t-32.*L23t*L13L13L23t+32.*size*L13L13L23L23t-96.*L3t*L13L13L3t+32.*size*L13L13L3L3t+64.*L1L13t*L2L23t+96.*L2L23t*L2L23t+16.*size*L2L23L2L23t-96.*L2t*L2L23L23t+64.*L12L13t*L2L3t+32.*L2L3t*L2L3t-16.*size*L2L3L2L3t-32.*L2t*L2L3L3t-96.*L23t*L2L2L23t+32.*size*L2L2L23L23t-32.*L3t*L2L2L3t+32.*size*L2L2L3L3t-64.*L1L12t*L23L3t+96.*L23L3t*L23L3t+16.*size*L23L3L23L3t-96.*L23t*L23L3L3t-96.*L3t*L23L23L3t+32.*size*L23L23L3L3t+16.*L2t*H0L1L12t+16.*L3t*H0L1L13t+32.*L12t*H0L1L2t+16.*H123t*H0L1L23t+32.*L13t*H0L1L3t-16.*L3t*H0L12H123t+16.*L1t*H0L12L2t+16.*H123t*H0L12L3t+16.*L2t*H0H123L13t-16.*L1t*H0H123L23t-16.*H123t*H0L13L2t+16.*L1t*H0L13L3t+16.*L3t*H0L2L23t+32.*L23t*H0L2L3t+16.*L2t*H0L23L3t-16.*L2t*L1H0L12t-32.*L23t*L1H0H123t-16.*L3t*L1H0L13t-64.*L12t*L1H0L2t+16.*H123t*L1H0L23t-64.*L13t*L1H0L3t+16.*H0t*L1L12L2t+16.*L3t*L1L12L23t+16.*L23t*L1L12L3t+16.*L3t*L1H123L2t+16.*H0t*L1H123L23t-16.*L2t*L1H123L3t-16.*L23t*L1L13L2t-16.*L2t*L1L13L23t+16.*H0t*L1L13L3t+32.*H123t*L1L2L3t+16.*L2t*L12H0L1t-16.*L3t*L12H0H123t-32.*L23t*L12H0L13t-16.*L1t*L12H0L2t+32.*L13t*L12H0L23t+16.*H123t*L12H0L3t-16.*L2t*L12L1H0t-32.*L13t*L12L1H123t-32.*H123t*L12L1L13t-16.*H0t*L12L1L2t+16.*L3t*L12L1L23t+16.*L23t*L12L1L3t+32.*L1t*L12H123L13t+16.*L23t*L12H123L2t+32.*L2t*L12H123L23t-16.*L3t*L12L13L2t-32.*H0t*L12L13L23t-16.*H123t*L12L2L23t-16.*L13t*L12L2L3t-32.*L23t*H123H0L1t-16.*L3t*H123H0L12t+16.*L2t*H123H0L13t+32.*L13t*H123H0L2t-16.*L1t*H123H0L23t-32.*L12t*H123H0L3t+32.*L13t*H123L1L12t-32.*L12t*H123L1L13t-16.*L3t*H123L1L2t+16.*H0t*H123L1L23t+16.*L2t*H123L1L3t-16.*L3t*H123L12H0t-64.*L1t*H123L12L13t-16.*L23t*H123L12L2t-64.*L2t*H123L12L23t+32.*H0t*H123L12L3t-16.*H0t*H123L13L2t-32.*L3t*H123L13L23t-16.*L23t*H123L13L3t-16.*L12t*H123L2L23t+16.*L13t*H123L23L3t+16.*L3t*L13H0L1t+32.*L23t*L13H0L12t+16.*L2t*L13H0H123t-16.*H123t*L13H0L2t-32.*L12t*L13H0L23t-16.*L1t*L13H0L3t-16.*L3t*L13L1H0t+32.*H123t*L13L1L12t+32.*L12t*L13L1H123t-16.*L23t*L13L1L2t-16.*L2t*L13L1L23t-16.*H0t*L13L1L3t+64.*L1t*L13L12H123t-16.*L3t*L13L12L2t+64.*H0t*L13L12L23t-32.*L2t*L13L12L3t+16.*L2t*L13H123H0t-32.*L1t*L13H123L12t-16.*H0t*L13H123L2t+64.*L3t*L13H123L23t+16.*L23t*L13H123L3t-16.*L12t*L13L2L3t+16.*H123t*L13L23L3t+64.*L12t*L2H0L1t+16.*L1t*L2H0L12t+32.*L13t*L2H0H123t-16.*H123t*L2H0L13t-16.*L3t*L2H0L23t-64.*L23t*L2H0L3t-32.*L12t*L2L1H0t+16.*H0t*L2L1L12t+16.*L3t*L2L1H123t-16.*L23t*L2L1L13t-32.*L13t*L2L1L23t-64.*H123t*L2L1L3t-16.*L1t*L2L12H0t-16.*H0t*L2L12L1t+16.*L23t*L2L12H123t-16.*L3t*L2L12L13t+16.*H123t*L2L12L23t-16.*L13t*L2L12L3t-16.*L3t*L2H123L1t-16.*L23t*L2H123L12t-16.*H0t*L2H123L13t+16.*L12t*L2H123L23t+32.*L1t*L2H123L3t-16.*H123t*L2L13H0t-16.*L23t*L2L13L1t-16.*L3t*L2L13L12t-16.*H0t*L2L13H123t-32.*L1t*L2L13L23t-16.*L12t*L2L13L3t+16.*H0t*L2L23L3t+16.*H123t*L23H0L1t-32.*L13t*L23H0L12t-16.*L1t*L23H0H123t+32.*L12t*L23H0L13t+16.*L3t*L23H0L2t-16.*L2t*L23H0L3t+16.*H123t*L23L1H0t+16.*L3t*L23L1L12t+16.*H0t*L23L1H123t-16.*L2t*L23L1L13t-32.*L13t*L23L1L2t+32.*L12t*L23L1L3t+16.*L3t*L23L12L1t+64.*L2t*L23L12H123t-64.*H0t*L23L12L13t-16.*H123t*L23L12L2t+32.*L1t*L23L12L3t-16.*L1t*L23H123H0t+16.*H0t*L23H123L1t-32.*L2t*L23H123L12t-64.*L3t*L23H123L13t-16.*L12t*L23H123L2t-16.*L13t*L23H123L3t-16.*L2t*L23L13L1t+32.*H0t*L23L13L12t+32.*L3t*L23L13H123t-32.*L1t*L23L13L2t-16.*H123t*L23L13L3t-16.*L3t*L23L2H0t+16.*H123t*L23L2L12t+16.*L12t*L23L2H123t-16.*H0t*L23L2L3t+64.*L13t*L3H0L1t+16.*H123t*L3H0L12t-32.*L12t*L3H0H123t+16.*L1t*L3H0L13t+64.*L23t*L3H0L2t+16.*L2t*L3H0L23t-32.*L13t*L3L1H0t+16.*L23t*L3L1L12t-16.*L2t*L3L1H123t+16.*H0t*L3L1L13t+64.*H123t*L3L1L2t+32.*L12t*L3L1L23t+16.*H123t*L3L12H0t+16.*L23t*L3L12L1t+32.*H0t*L3L12H123t-32.*L2t*L3L12L13t-16.*L13t*L3L12L2t+32.*L1t*L3L12L23t+16.*L2t*L3H123L1t-16.*L23t*L3H123L13t-32.*L1t*L3H123L2t+16.*L13t*L3H123L23t-16.*L1t*L3L13H0t-16.*H0t*L3L13L1t+16.*L23t*L3L13H123t-16.*L12t*L3L13L2t+16.*H123t*L3L13L23t-32.*L23t*L3L2H0t-32.*H123t*L3L2L1t-16.*L13t*L3L2L12t-16.*L12t*L3L2L13t+16.*H0t*L3L2L23t-16.*L2t*L3L23H0t-16.*L13t*L3L23H123t-16.*H123t*L3L23L13t-16.*H0t*L3L23L2t+8.*size*H0L2L12L1t+8.*size*H0L2L13H123t+8.*size*H0L23H123L1t+8.*size*H0L23L13L12t+8.*size*H0L3H123L12t+8.*size*H0L3L13L1t+8.*size*H0L3L23L2t+8.*size*L1H0L2L12t+8.*size*L1H0L23H123t+8.*size*L1H0L3L13t-32.*size*L1L12L2H0t+32.*size*L1H123L23H0t+8.*size*L1L13H123L12t-32.*size*L1L13L3H0t-32.*size*L1L2L12H0t-32.*size*L1L23H123H0t+8.*size*L1L23L2L13t-32.*size*L1L3L13H0t-8.*size*L1L3L2H123t-8.*size*L1L3L23L12t+32.*size*L12H0L2L1t+8.*size*L12H0L23L13t+8.*size*L12H0L3H123t+16.*size*L12L1H0L2t+8.*size*L12L1L13H123t-32.*size*L12L1L2H0t-8.*size*L12L1L3L23t-32.*size*L12H123L13L1t+32.*size*L12H123L3H0t+32.*size*L12L13H123L1t-32.*size*L12L13L23H0t+32.*size*L12L23L13H0t+8.*size*L12L23L2H123t-32.*size*L12L23L3L1t+32.*size*L12L3H123H0t-8.*size*L12L3L2L13t-32.*size*L12L3L23L1t+8.*size*H123H0L2L13t+32.*size*H123H0L23L1t-32.*size*H123H0L3L12t+16.*size*H123L1H0L23t-32.*size*H123L1L13L12t+32.*size*H123L1L23H0t-8.*size*H123L1L3L2t+16.*size*H123L12H0L3t+16.*size*H123L12L1L13t+32.*size*H123L12L13L1t+8.*size*H123L12L23L2t-32.*size*H123L12L3H0t+32.*size*H123L13L2H0t-32.*size*H123L2L13H0t-32.*size*H123L2L23L12t+32.*size*H123L2L3L1t-32.*size*H123L23L2L12t-32.*size*H123L3L2L1t-8.*size*H123L3L23L13t-32.*size*L13H0L2H123t-32.*size*L13H0L23L12t+32.*size*L13H0L3L1t+16.*size*L13L1H0L3t+8.*size*L13L1L23L2t-32.*size*L13L1L3H0t+16.*size*L13L12H0L23t-32.*size*L13L12L1H123t+32.*size*L13L12L23H0t-8.*size*L13L12L3L2t+16.*size*L13H123H0L2t-32.*size*L13H123L2H0t-8.*size*L13H123L3L23t+32.*size*L13L2L23L1t-32.*size*L13L2L3L12t+32.*size*L13L23L2L1t+32.*size*L13L23L3H123t+32.*size*L13L3L2L12t+32.*size*L13L3L23H123t+8.*size*L2H0L3L23t+32.*size*L2L1H0L12t-32.*size*L2L1L23L13t-32.*size*L2L1L3H123t+32.*size*L2L12L23H123t+32.*size*L2L12L3L13t-32.*size*L2H123H0L13t-16.*size*L2H123L1L3t+16.*size*L2H123L12L23t-32.*size*L2H123L23L12t+32.*size*L2H123L3L1t+16.*size*L2L13L1L23t-16.*size*L2L13L12L3t-32.*size*L2L13L23L1t+32.*size*L2L13L3L12t-32.*size*L2L23L3H0t-32.*size*L2L3L23H0t+32.*size*L23H0L3L2t-32.*size*L23L1H0H123t+32.*size*L23L1L3L12t-32.*size*L23L12H0L13t-16.*size*L23L12L1L3t+32.*size*L23L12L3L1t+32.*size*L23H123L12L2t+32.*size*L23H123L3L13t+32.*size*L23L13L1L2t-16.*size*L23L13H123L3t-32.*size*L23L13L3H123t+16.*size*L23L2H0L3t-32.*size*L23L2L3H0t+32.*size*L3L1H0L13t+32.*size*L3L12H0H123t-32.*size*L3L12L1L23t+32.*size*L3H123L1L2t+32.*size*L3L13L12L2t-32.*size*L3L13H123L23t+32.*size*L3L2H0L23t;

	}
	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}

//	std::cout<<"temp1"<< temp1.trace()<<std::endl;
//	std::cout<<"temp2"<< temp2.trace()<<std::endl;

	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();


	return S;
}

double Dirac::action10()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd T2=H0*H0;
	std::complex<double> T2t=T2.trace();
	std::complex<double> H0t=H0.trace();


	S2= 2.*size*T2t+2.*H0t*H0t;
	S4=0;
	if(gD4!=0)
	{
		Eigen::MatrixXcd T3=T2*H0;
		Eigen::MatrixXcd T4=T2*T2;
		S4=2.*size*T4.trace()+8.*H0t*T3.trace()+6.*T2t*T2t;

	}

	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}

//	std::cout<<"temp1"<< temp1.trace()<<std::endl;
//	std::cout<<"temp2"<< temp2.trace()<<std::endl;

	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();


	return S;
}

double Dirac::action01()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd T2=H0*H0;
	std::complex<double> T2t=T2.trace();
	std::complex<double> L1t=H0.trace();


	S2= 2.*size*T2.trace()-2.*L1t*L1t;
	S4=0;
	if(gD4!=0)
	{
		Eigen::MatrixXcd T3=T2*H0;
		Eigen::MatrixXcd T4=T2*T2;
		S4=2.*size*T4.trace()-8.*L1t*T3.trace()+6.*T2t*T2t;

	}

	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}

//	std::cout<<"temp1"<< temp1.trace()<<std::endl;
//	std::cout<<"temp2"<< temp2.trace()<<std::endl;

	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();


	return S;
}


double Dirac::action20()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd H02=H0*H0;
	Eigen::MatrixXcd H1232=H123*H123;
	std::complex<double> H02t=H02.trace();
	std::complex<double> H1232t=H1232.trace();
	std::complex<double> H0t=H0.trace();
	std::complex<double> H123t=H123.trace();

	S2= 4.*size*H02t+4.*size*H1232t+4.*H0t*H0t+4.*H123t*H123t;

	S4=0;
	if(gD4!=0)
	{
		Eigen::MatrixXcd H04=H02*H02;
		Eigen::MatrixXcd H1234=H1232*H1232;
		Eigen::MatrixXcd H0H123=H0*H123;
		Eigen::MatrixXcd H0212=H0H123*H0H123;
		Eigen::MatrixXcd H0122=H02*H1232;
		Eigen::MatrixXcd H03=H02*H0;
		Eigen::MatrixXcd H1232H0=H1232*H0;
		Eigen::MatrixXcd H1233=H1232*H123;
		Eigen::MatrixXcd H02H123=H02*H123;

		S4+=4.*size*(H04.trace()+H1234.trace()
		-2.*H0212.trace()
		+4.*H0122.trace());

		S4+=16.*H0t*(H03.trace()+H1232H0.trace())
		+16.*H123t*(H1233.trace()+H02H123.trace());

		S4+=12.*(H02t*H02t+H1232t*H1232t)
		+8.*H02t*H1232t
		+16.*H0H123.trace()*H0H123.trace();

	}

	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}




	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();


	return S;
}


double Dirac::action11()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd H12=H0*H0;
	Eigen::MatrixXcd L1p2=L1*L1;

	//Eigen::MatrixXcd H13,H14,L13,L14,H12L1p2,H1L1,H1L1p2,L1H12,H1L1H1L1;

	std::complex<double> H1t=H0.trace();
	std::complex<double> H12t=H12.trace();
	std::complex<double> L1t=L1.trace();
	std::complex<double> L1p2t=L1p2.trace();

	//mtemp.setIdentity();


	S2= 4.*size*(H12t-L1p2t)+4.*L1t*L1t+4.*H1t*H1t;
	S4=0;
	if(gD4!=0)
	{
		Eigen::MatrixXcd H14=H12*H12;
		Eigen::MatrixXcd L14=L1p2*L1p2;
		Eigen::MatrixXcd HL=H0*L1;
		Eigen::MatrixXcd HLHL=HL*HL;
		Eigen::MatrixXcd HHLL=H12*L1p2;
		Eigen::MatrixXcd H3=H12*H0;
		Eigen::MatrixXcd H2L=H12*L1;
		Eigen::MatrixXcd L1p3=L1p2*L1;
		Eigen::MatrixXcd L2H=L1p2*H0;

		S4+=4.*size*(H14.trace()+L14.trace()+2.*HLHL.trace()-4.*HHLL.trace());

		S4+=16.*H1t*(H3.trace()-L2H.trace())+16.*L1t*(-L1p3.trace()+H2L.trace());

		S4+=12.*(H12t*H12t+L1p2t*L1p2t)-8.*H12t*L1p2t+16.*HL.trace()*HL.trace();

	}

	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}

//	std::cout<<"temp1"<< temp1.trace()<<std::endl;
//	std::cout<<"temp2"<< temp2.trace()<<std::endl;

	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();


	return S;
}

double Dirac::action02()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd L1p2=L1*L1;
	Eigen::MatrixXcd L2p2=L2*L2;
	std::complex<double> L1p2t=L1p2.trace();
	std::complex<double> L2p2t=L2p2.trace();
	std::complex<double> L1t=L1.trace();
	std::complex<double> L2t=L2.trace();

	S2= -4.*size*L1p2t-4.*size*L2p2t+4.*L1t*L1t+4.*L2t*L2t;

	S4=0;
	if(gD4!=0)
	{
		Eigen::MatrixXcd L14=L1p2*L1p2;
		Eigen::MatrixXcd L24=L2p2*L2p2;
		Eigen::MatrixXcd L1L2=L1*L2;
		Eigen::MatrixXcd L1p212=L1L2*L1L2;
		Eigen::MatrixXcd L1122=L1p2*L2p2;
		Eigen::MatrixXcd L13=L1p2*L1;
		Eigen::MatrixXcd L2p2L1=L2p2*L1;
		Eigen::MatrixXcd L23=L2p2*L2;
		Eigen::MatrixXcd L1p2L2=L1p2*L2;

		S4+=4.*size*(L14.trace()+L24.trace()
		-2.*L1p212.trace()
		+4.*L1122.trace());

		S4+=-16.*L1t*(L13.trace()+L2p2L1.trace())
		-16.*L2t*(L23.trace()+L1p2L2.trace());

		S4+=12.*(L1p2t*L1p2t+L2p2t*L2p2t)
		+8.*L1p2t*L2p2t
		+16.*L1L2.trace()*L1L2.trace();

	}

	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}




	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();



	return S;
}

double Dirac::action03()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd L1p2=L1*L1;
	Eigen::MatrixXcd L2p2=L2*L2;
	Eigen::MatrixXcd L3p2=L3*L3;
	Eigen::MatrixXcd Hp2=H0*H0;

	std::complex<double> L1p2t=L1p2.trace();
	std::complex<double> L2p2t=L2p2.trace();
	std::complex<double> L3p2t=L3p2.trace();
	std::complex<double> Hp2t=Hp2.trace();
	std::complex<double> L1t=L1.trace();
	std::complex<double> L2t=L2.trace();
	std::complex<double> L3t=L3.trace();
	std::complex<double> Ht=H0.trace();

	S2= -4.*size*L1p2t-4.*size*L2p2t-4.*size*L3p2t
		+4.*size*Hp2t
		+4.*L1t*L1t+4.*L2t*L2t+4.*L3t*L3t
		+4.*Ht*Ht;


	S4=0;
	if(gD4!=0)
	{
		Eigen::MatrixXcd L1p4=L1p2*L1p2;
		Eigen::MatrixXcd L2p4=L2p2*L2p2;
		Eigen::MatrixXcd L3p4=L3p2*L3p2;
		Eigen::MatrixXcd Hp4=Hp2*Hp2;
		Eigen::MatrixXcd L1L2=L1*L2;
		Eigen::MatrixXcd L1L3=L1*L3;
		Eigen::MatrixXcd L2L3=L2*L3;
		Eigen::MatrixXcd HL1=H0*L1;
		Eigen::MatrixXcd HL2=H0*L2;
		Eigen::MatrixXcd HL3=H0*L3;

		Eigen::MatrixXcd HL1L2L3=HL1*L2L3;
		Eigen::MatrixXcd HL2L1L3=HL2*L1L3;
		Eigen::MatrixXcd HL3L1L2=HL3*L1L2;
		Eigen::MatrixXcd HL1L3L2=HL1*L3*L2;
		Eigen::MatrixXcd HL2L3L1=HL2*L3*L1;
		Eigen::MatrixXcd HL3L2L1=HL3*L2*L1;

		Eigen::MatrixXcd L1212=L1L2*L1L2;
		Eigen::MatrixXcd L1313=L1L3*L1L3;
		Eigen::MatrixXcd L2323=L2L3*L2L3;
		Eigen::MatrixXcd LH1H1=HL1*HL1;
		Eigen::MatrixXcd LH2H2=HL2*HL2;
		Eigen::MatrixXcd LH3H3=HL3*HL3;
		Eigen::MatrixXcd L1122=L1p2*L2p2;
		Eigen::MatrixXcd L1133=L1p2*L3p2;
		Eigen::MatrixXcd L2233=L3p2*L2p2;
		Eigen::MatrixXcd HH11=Hp2*L1p2;
		Eigen::MatrixXcd HH22=Hp2*L2p2;
		Eigen::MatrixXcd HH33=Hp2*L3p2;
		Eigen::MatrixXcd L1p3=L1p2*L1;
		Eigen::MatrixXcd L2p3=L2p2*L2;
		Eigen::MatrixXcd L3p3=L3p2*L3;
		Eigen::MatrixXcd Hp3=Hp2*H0;
		Eigen::MatrixXcd L2p2L1=L2p2*L1;
		Eigen::MatrixXcd L3p2L1=L3p2*L1;
		Eigen::MatrixXcd Hp2L1=Hp2*L1;
		Eigen::MatrixXcd L1p2L2=L1p2*L2;
		Eigen::MatrixXcd L3p2L2=L3p2*L2;
		Eigen::MatrixXcd Hp2L2=Hp2*L2;
		Eigen::MatrixXcd L1p2L3=L1p2*L3;
		Eigen::MatrixXcd L2p2L3=L2p2*L3;
		Eigen::MatrixXcd Hp2L3=Hp2*L3;
		Eigen::MatrixXcd L1p2H=L1p2*H0;
		Eigen::MatrixXcd L2p2H=L2p2*H0;
		Eigen::MatrixXcd L3p2H=L3p2*H0;
		Eigen::MatrixXcd HcL2L3=H0*(L2*L3-L3*L2);
		Eigen::MatrixXcd HcL1L3=H0*(L1*L3-L3*L1);
		Eigen::MatrixXcd HcL1L2=H0*(L1*L2-L2*L1);
		Eigen::MatrixXcd L1cL2L3=L1*(L2*L3-L3*L2);



		S4+=4.*size*(L1p4.trace()
		+L2p4.trace()
		+L3p4.trace()
		+Hp4.trace()
		-2.*L1212.trace()
		-2.*L1313.trace()
		-2.*L2323.trace()
		-2.*LH1H1.trace()
		-2.*LH2H2.trace()
		-2.*LH3H3.trace()
		+4.*L1122.trace()
		+4.*L1133.trace()
		+4.*L2233.trace()
		-4.*HH11.trace()
		-4.*HH22.trace()
		-4.*HH33.trace()
		-4.*HL1L2L3.trace()
		+4.*HL1L3L2.trace()
		+4.*HL2L1L3.trace()
		-4.*HL2L3L1.trace()
		+4.*HL3L2L1.trace()
		-4.*HL3L1L2.trace()
		);



		S4+=-16.*L1t*(L1p3.trace()+L2p2L1.trace()+L3p2L1.trace()-3.*Hp2L1.trace())
		-16.*L2t*(L2p3.trace()+L1p2L2.trace()+L3p2L2.trace()-3.*Hp2L2.trace())
		-16.*L3t*(L3p3.trace()+L1p2L3.trace()+L2p2L3.trace()-3.*Hp2L3.trace())
		+16.*Ht*(Hp3.trace()-3.*L1p2H.trace()-3.*L2p2H.trace()-3.*L3p2H.trace());


	S4+=16.*(L1t*HcL2L3.trace()-L2t*HcL1L3.trace()+L3t*HcL1L2.trace())-48.*Ht*L1cL2L3.trace();


		S4+=
		 12.*(L1p2t*L1p2t)
		+12.*(L2p2t*L2p2t)
		+12.*(L3p2t*L3p2t)
		+12.*(Hp2t*Hp2t)
		+8.*L1p2t*L2p2t
		+8.*L1p2t*L3p2t
		+8.*L2p2t*L3p2t
		-24.*Hp2t*L1p2t
		-24.*Hp2t*L2p2t
		-24.*Hp2t*L3p2t
		+16.*L1L2.trace()*L1L2.trace()
		+16.*L1L3.trace()*L1L3.trace()
		+16.*L2L3.trace()*L2L3.trace()
		+48.*HL1.trace()*HL1.trace()
		+48.*HL2.trace()*HL2.trace()
		+48.*HL3.trace()*HL3.trace();

	/*	std::cout<< 12.*(L1p2t*L1p2t)<< " " <<
		+12.*(L2p2t*L2p2t)<< " " <<
		+12.*(L3p2t*L3p2t)<< " " <<
		+12.*(Hp2t*Hp2t)<< " " <<
		+8.*L1p2t*L2p2t<< " " <<
		+8.*L1p2t*L3p2t<< " " <<
		+8.*L2p2t*L3p2t<< " " <<
		-24.*Hp2t*L1p2t<< " " <<
		-24.*Hp2t*L2p2t<< " " <<
		-24.*Hp2t*L3p2t<< " " <<
		+16.*L1L2.trace()*L1L2.trace()<< " " <<
		+16.*L1L3.trace()*L1L3.trace()<< " " <<
		+16.*L2L3.trace()*L2L3.trace()<< " " <<
		+48.*HL1.trace()*HL1.trace()<< " " <<
		+48.*HL2.trace()*HL2.trace()<< " " <<
		+48.*HL3.trace()*HL3.trace()<< std::endl; */
	}

	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}




	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();



	return S;
}


double Dirac::action031()
{
	std::complex<double> S2,S4;

	// New simpler action, we just want D^2

	Eigen::MatrixXcd L1p2=L1*L1;
	Eigen::MatrixXcd L2p2=L2*L2;
	Eigen::MatrixXcd L3p2=L3*L3;


	std::complex<double> L1p2t=L1p2.trace();
	std::complex<double> L2p2t=L2p2.trace();
	std::complex<double> L3p2t=L3p2.trace();
	std::complex<double> L1t=L1.trace();
	std::complex<double> L2t=L2.trace();
	std::complex<double> L3t=L3.trace();


	S2= -4.*size*L1p2t-4.*size*L2p2t-4.*size*L3p2t
		+4.*L1t*L1t+4.*L2t*L2t+4.*L3t*L3t;

	S4=0;
	if(gD4!=0)
	{
		Eigen::MatrixXcd L1p4=L1p2*L1p2;
		Eigen::MatrixXcd L2p4=L2p2*L2p2;
		Eigen::MatrixXcd L3p4=L3p2*L3p2;
		Eigen::MatrixXcd L1L2=L1*L2;
		Eigen::MatrixXcd L1L3=L1*L3;
		Eigen::MatrixXcd L2L3=L2*L3;

		Eigen::MatrixXcd L1212=L1L2*L1L2;
		Eigen::MatrixXcd L1313=L1L3*L1L3;
		Eigen::MatrixXcd L2323=L2L3*L2L3;
		Eigen::MatrixXcd L1122=L1p2*L2p2;
		Eigen::MatrixXcd L1133=L1p2*L3p2;
		Eigen::MatrixXcd L2233=L3p2*L2p2;
		Eigen::MatrixXcd L1p3=L1p2*L1;
		Eigen::MatrixXcd L2p3=L2p2*L2;
		Eigen::MatrixXcd L3p3=L3p2*L3;
		Eigen::MatrixXcd L2p2L1=L2p2*L1;
		Eigen::MatrixXcd L3p2L1=L3p2*L1;
		Eigen::MatrixXcd L1p2L2=L1p2*L2;
		Eigen::MatrixXcd L3p2L2=L3p2*L2;
		Eigen::MatrixXcd L1p2L3=L1p2*L3;
		Eigen::MatrixXcd L2p2L3=L2p2*L3;



		S4+=4.*size*(L1p4.trace()
		+L2p4.trace()
		+L3p4.trace()
		-2.*L1212.trace()
		-2.*L1313.trace()
		-2.*L2323.trace()
		+4.*L1122.trace()
		+4.*L1133.trace()
		+4.*L2233.trace()	);



		S4+=-16.*L1t*(L1p3.trace()+L2p2L1.trace()+L3p2L1.trace())
		-16.*L2t*(L2p3.trace()+L1p2L2.trace()+L3p2L2.trace())
		-16.*L3t*(L3p3.trace()+L1p2L3.trace()+L2p2L3.trace());



		S4+=
		 12.*(L1p2t*L1p2t)
		+12.*(L2p2t*L2p2t)
		+12.*(L3p2t*L3p2t)
		+8.*L1p2t*L2p2t
		+8.*L1p2t*L3p2t
		+8.*L2p2t*L3p2t
		+16.*L1L2.trace()*L1L2.trace()
		+16.*L1L3.trace()*L1L3.trace()
		+16.*L2L3.trace()*L2L3.trace();

	/*	std::cout<< 12.*(L1p2t*L1p2t)<< " " <<
		+12.*(L2p2t*L2p2t)<< " " <<
		+12.*(L3p2t*L3p2t)<< " " <<
		+12.*(Hp2t*Hp2t)<< " " <<
		+8.*L1p2t*L2p2t<< " " <<
		+8.*L1p2t*L3p2t<< " " <<
		+8.*L2p2t*L3p2t<< " " <<
		-24.*Hp2t*L1p2t<< " " <<
		-24.*Hp2t*L2p2t<< " " <<
		-24.*Hp2t*L3p2t<< " " <<
		+16.*L1L2.trace()*L1L2.trace()<< " " <<
		+16.*L1L3.trace()*L1L3.trace()<< " " <<
		+16.*L2L3.trace()*L2L3.trace()<< " " <<
		+48.*HL1.trace()*HL1.trace()<< " " <<
		+48.*HL2.trace()*HL2.trace()<< " " <<
		+48.*HL3.trace()*HL3.trace()<< std::endl; */
	}

	if(abs(S2.imag())>1e-10||abs(S4.imag())>1e-10)
	{
		std::cout<<"This sucks, your action isn't real o.O"<<std::endl;
		std::cout<< S2.imag()<<std::endl;
		std::cout<< S4.imag()<<std::endl;
	}




	S= gD2*S2.real()+gD22*S2.real()*S2.real()+gD4*S4.real();



	return S;
}
