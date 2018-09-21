#include <inttypes.h>
#include <fstream>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <vector>
#include "Riostream.h"

#include "jKiss64.hh"

using namespace std;

#define nm 1.;
#define um 1000.*nm;
#define mm 1000.*um;
#define cm 10.*mm;
#define m 1000.*mm;

#define g 1.;
#define kg 1000.*g;

#define ns 1.;
#define us 1000.*ns;
#define ms 1000.*us;
#define s 1000.*ms;

#define eV 1.;
#define keV 1000.*eV;
#define MeV 1000.*keV;

#define V 1.;

#define Gy 1.;

double gap_size;
double HV;
int N_pas;
double T_pas;
double energie_primaire;
double dose_initiale;
double potentiel_ionisation;

struct bucket 
{
	int nature; 				//nature du seau, 0 électron; 1 ion positif; 2 ion négatif
	int intervalle; 		//intervalle dans lequel il se trouve
	double quantum; 	  //charge totale
	double position;	  //position en x
	double dispersion;	//dispersion en x
};

void EntryParameters(int config_simu)
{
	std::vector<string> Variable_init;
	std::vector<double> Value_init;

	stringstream ss;
	ss<<config_simu;
	string indice=ss.str();

	int iint;
	int ichar;
	int ivar;
	char next;
	string data_file;
	string path_file;
	string filename="./Entry/Entry_param_"+indice+".txt";
	// string pathname="./Entry/Entry_path.txt";
	// string filename="./Entry/Entry_param_1.txt";
	string tmp;

	Variable_init.clear();
	Value_init.clear();
	Variable_init.push_back("Taille de gap (mm)"); // 0
	Value_init.push_back(0.);	//	0
	Variable_init.push_back("Nombre de pas"); // 1
	Value_init.push_back(0.);	//	1
	Variable_init.push_back("Haute tension (V)"); // 2
	Value_init.push_back(0.);	//	2
	Variable_init.push_back("Pas temporel (ns):"); // 3
	Value_init.push_back(0.);	//	3
	Variable_init.push_back("Dose (Gy/pulse):"); // 5
	Value_init.push_back(0.);	//	5

	cout<<endl;
	ifstream datafile_param(filename.c_str());
	// ifstream datafile_param("./Entry/Entry_param_1.txt");
	if(!datafile_param)
	{
		cout<<"No entry parameters file for this configuration: "<<config_simu<<endl;
		cout<<"No file: "<<filename<<endl;
		exit(0);
	}
	else
	{
		while(true)
		{
			string buffer;
			string variable;
			double value;
			getline(datafile_param,tmp);
			for(ichar=0;ichar<tmp.size();ichar++)
			{
				next=tmp[ichar];
				if(next!=':')
					variable+=next;
				else
				{	
					for(iint=ichar+2;iint<tmp.size();iint++)
					{	
						next=tmp[iint];
						buffer+=next;
					}
					ichar=tmp.size();	
				}	
			}
			int ind_value=-1;
			for(ivar=0;ivar<Variable_init.size();ivar++)
			{
				if(!variable.compare(Variable_init[ivar]))
					ind_value=ivar;
			}
			if(ind_value==-1)
				cout<<"Variable "<<variable<<" not defined. Maybe check the entry file."<<endl;
			else
			{
				if(buffer.find_first_not_of(' ')!=std::string::npos)
				{
					if(ind_value==0)
					{
						gap_size=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" mm; nouvelle valeur: "<<gap_size<<" mm"<<endl;
					}
					if(ind_value==1)
					{
						N_pas=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<"; nouvelle valeur: "<<N_pas<<endl;
					}
					if(ind_value==2)
					{
						HV=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" V; nouvelle valeur: "<<HV<<" V"<<endl;
					}
					if(ind_value==3)
					{
						T_pas=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" ns; nouvelle valeur: "<<T_pas<<" ns"<<endl;
					}
					if(ind_value==5)
					{
						dose_initiale=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" Gy/pulse; nouvelle valeur: "<<dose_initiale<<" Gy/pulse"<<endl;
					}
				}
				else
					cout<<"Value not define for "<<variable<<". Keeping the default one: "<<Value_init[ind_value]<<". Maybe check the entry file."<<endl;
			}	
			if(datafile_param.eof()) break;
		}
		datafile_param.close();
	}
	gap_size*=mm;
	HV*=V;
	T_pas*=ns;
	dose_initiale*=Gy;


}

double Uniform()
{
	double	uni;
	u64 x;
	x=randk();
	uni=(double)x;
	uni/=pow(2.,64.);
	return	uni;
}

double Gaussian(double mean,double rms)
{
	double	data;
	double	U1,U2,Disp,Fx;
  do{
  	U1=2.*Uniform()-1.;
  	U2=2.*Uniform()-1.;
  	Fx=U1*U1+U2*U2;
  }while(Fx>=1.);
  Fx=sqrt((-2.*log(Fx))/Fx);
  Disp=U1*Fx;
  data=mean+Disp*rms;
	return	data;
}

int Ionisation()
{
	//valeurs pour électrons de 4.5 MeV
	double stopping_power_air=1.812*MeV*cm*cm/g; //MeV.cm2.g-1
	double stopping_power_water=1.882*MeV*cm*cm/g; //MeV.cm2.g-1
	double densite_air=1.225*kg/(m*m*m); //
}

int main()
{
	double determinant=0.;
	double **Matrice=(double**)malloc(N_pas*sizeof(double*));
	double **MatriceInv=(double**)malloc(N_pas*sizeof(double*));
	for(int i=0;i<N_pas;i++)
	{
		Matrice[i]=(double*)malloc(N_pas*sizeof(double));
		MatriceInv[i]=(double*)malloc(N_pas*sizeof(double));
	}

	std::vector<bucket> buck;

	for(int i=0;i<N_pas;i++)
	{
		free(Matrice[i]);  
		free(MatriceInv[i]);  
	}
	free(Matrice);
	free(MatriceInv);
}