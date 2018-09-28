// Simulation of Simple Ionization Chamber
// et ça se prononce "saucisse"

#include <algorithm> 
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
// #include <Python.h>
// #include <matplotlib.h>

#include "jKiss64.hh"
#include "Vector.hh"

using namespace std;

#define e_c 1.60217646E-19
#define electron 	0
#define cation 		1
#define anion 		2

#define air 	2728
#define azote 27

#define mort		0
#define vivant	1
#define pending	2
#define traite	3
#define dehors	4

double nm=1.E-6;
double um=1000.*nm;
double mm=1000.*um;
double cm=10.*mm;
double m=1000.*mm;

double g=1.;
double kg=1000.*g;

double ns=1.;
double us=1000.*ns;
double ms=1000.*us;
double s=1000.*ms;

double eV=1.E-6;
double keV=1000.*eV;
double MeV=1000.*keV;

double V=1.;

double Gy=1./e_c*eV/kg;

int nature_gaz;
double gap_size;
double HV;
int N_pas;
int B_pas;
double T_pas;
double dose_initiale;
double duree_pulse;
double duree_acqui;

bool faisceau;
double mobilite_electron;
double mobilite_cation;
double mobilite_anion;
double temps_attachement;
double recombinaison;

struct bucket 
{
	int condition;
	int nature; 				//nature du seau, 0 électron; 1 ion positif; 2 ion négatif
	int intervalle; 		//intervalle dans lequel il se trouve
	double charge;			//-1 ou 1;
	double quanta; 	  	//charge totale
	double position;	  //position en x
};

void PrintUnits()
{
	if(ns==1.)	cout<<"Temps : ns"<<endl;
	if(us==1.)	cout<<"Temps : us"<<endl;
	if(ms==1.)	cout<<"Temps : ms"<<endl;
	if(s==1.)		cout<<"Temps : s"<<endl;

	if(nm==1.)	cout<<"Dist. : nm"<<endl;
	if(um==1.)	cout<<"Dist. : um"<<endl;
	if(mm==1.)	cout<<"Dist. : mm"<<endl;
	if(cm==1.)	cout<<"Dist. : cm"<<endl;
	if(m==1.)		cout<<"Dist. : m"<<endl;

	if(V==1.)		cout<<"Tens. : V"<<endl;

	if(g==1.)		cout<<"Masse : g"<<endl;
	if(kg==1.)	cout<<"Masse : kg"<<endl;

	if(eV==1.)	cout<<"Ener. : eV"<<endl;
	if(keV==1.)	cout<<"Ener. : keV"<<endl;
	if(MeV==1.)	cout<<"Ener. : MeV"<<endl;
}

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
	Variable_init.push_back("Pas temporel (ns)"); // 3
	Value_init.push_back(0.);	//	3
	Variable_init.push_back("Dose (Gy/pulse)"); // 4
	Value_init.push_back(0.);	//	4
	Variable_init.push_back("Durée du pulse (us)"); // 5
	Value_init.push_back(0.);	//	5
	Variable_init.push_back("Nombre de seaux par intervalle"); // 6
	Value_init.push_back(0.);	//	6
	Variable_init.push_back("Gaz"); // 7
	Value_init.push_back(0.);	//	7
	Variable_init.push_back("Durée acquisition (us)"); // 8
	Value_init.push_back(0.);	//	8

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
					if(ind_value==7)
					{
						if(!buffer.compare("N2")||!buffer.compare("Azote")||!buffer.compare("azote"))
						{
							nature_gaz=azote;
							cout<<"Azote"<<endl;
						}
						if(!buffer.compare("Air")||!buffer.compare("air"))
						{
							nature_gaz=air;
							cout<<"Air"<<endl;
						}
					}
					if(ind_value==0)
					{
						gap_size=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" mm; nouvelle valeur: "<<gap_size<<" mm"<<endl;
						gap_size*=mm;
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
						HV*=V;
					}
					if(ind_value==3)
					{
						T_pas=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" ns; nouvelle valeur: "<<T_pas<<" ns"<<endl;
						T_pas*=ns;
					}
					if(ind_value==4)
					{
						dose_initiale=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" Gy/pulse; nouvelle valeur: "<<dose_initiale<<" Gy/pulse"<<endl;
						dose_initiale*=Gy;
					}
					if(ind_value==5)
					{
						duree_pulse=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" us; nouvelle valeur: "<<duree_pulse<<" us"<<endl;
						duree_pulse*=us;
					}
					if(ind_value==6)
					{
						B_pas=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<"; nouvelle valeur: "<<B_pas<<endl;
					}
					if(ind_value==8)
					{
						duree_acqui=(double)atof(buffer.c_str());
						cout<<Variable_init[ind_value]<<" valeur par défaut: "<<Value_init[ind_value]<<" us; nouvelle valeur: "<<duree_acqui<<" us"<<endl;
						duree_acqui*=us;
					}
				}
				else
					cout<<"Value not define for "<<variable<<". Keeping the default one: "<<Value_init[ind_value]<<". Maybe check the entry file."<<endl;
			}	
			if(datafile_param.eof()) break;
		}
		datafile_param.close();
	}
}

double Uniform()
{
	double uni;
	u64 x;
	x=randk();
	uni=(double)x;
	uni/=pow(2.,64.);
	return	uni;
}

double Gaussian(double mean,double rms)
{
	double data;
	double U1,U2,Disp,Fx;
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

double Ionisation()
{
	double potentiel_ionisation=33.97*eV;
	//valeurs pour électrons de 4.5 MeV
	double stopping_power_air=1.812*MeV*cm*cm/g; //MeV.cm2.g-1
	double stopping_power_water=1.882*MeV*cm*cm/g; //MeV.cm2.g-1
	double densite_eau=1000.*kg/(m*m*m); 
	double densite_air=1.225*kg/(m*m*m); 
	double surface_electrode=1.*cm*cm;
	double volume_electrode=surface_electrode*gap_size;
	double dose_par_volume=dose_initiale*densite_eau;
	double intensite_par_unite_surface=dose_par_volume/(densite_eau*stopping_power_water);
	double nbr_ionisation=intensite_par_unite_surface*stopping_power_air*densite_air*volume_electrode/potentiel_ionisation;

	nbr_ionisation*=T_pas/duree_pulse;	//nombre d'ionisations par pas de temps
	nbr_ionisation/=N_pas;	//nombre d'ionisations par intervalle d'espace
	cout<<"Ionisations par itération de temps : "<<nbr_ionisation<<endl;
	return nbr_ionisation;
}

double AttachementTime(double E)
{
	double Ta0=95.24*ns;
	double E0=258.5*V/mm;
	switch(nature_gaz)
	{
		case(azote):
			return 0.;
		break;	
		case(air):
			return Ta0*(1.-exp(-E/E0));
		break;	
	}
	return 0.;
}

double ElectronSpeed(double E)
{
	double v;
	double v1,v2;
	double E1,E2;
	switch(nature_gaz)
	{
		case azote:
			v1=7.741E7*mm/s;
			v2=2.162E6*mm/s;
			E1=1.076E3*V/mm;
			E2=6.412*V/mm;
		break;
		case air:
			v1=1.717E8*mm/s;
			v2=6.499E6*mm/s;
			E1=2.819E3*V/mm;
			E2=12.72*V/mm;
		break;	
	}
	return v1*(1-exp(-E/E1))+v2*(1-exp(-E/E2));
}

double IonSpeed(int n_part)
{
	switch(nature_gaz)
	{
		case azote:
			return 188.*mm*mm/(V*s);
		break;
		case air:
			switch(n_part)
			{
				case cation:
					return 186*mm*mm/(V*s);
				break;
				case anion:
					return 209*mm*mm/(V*s);
				break;
			}
		break;	
	}
	return 0.;
}

double Dispersion(double mob)
{
	return sqrt(0.052*V*mob*T_pas);
}

int IntervalLocation(double position)
{
	int interv;
	double espace=gap_size/N_pas;
	interv=ceil(position/espace);
	interv=std::max(interv,0);
	interv=std::min(interv,N_pas);
	return interv;
}

vvector ElectricField(vvector rho)
{
	vvector ChElec=vector_create(N_pas);
	for(int ind=0;ind<N_pas;ind++)
		ChElec->coef[ind]=HV/gap_size;
	return ChElec;
}

void Python(std::vector<bucket>& seaux)
{
	int sizeofbucket=seaux.size();
	int tot_elec=0;
	int tot_cati=0;
	int tot_anio=0;
	ofstream output_file_e("./Output_e.txt");
	ofstream output_file_c("./Output_c.txt");
	for(int i=0;i<sizeofbucket;i++)
	{
		// output_file<<seaux[i].quanta<<endl;
		if(seaux[i].nature==electron)
		{
			output_file_e<<seaux[i].position<<endl;
			tot_elec++;
		}
		if(seaux[i].nature==cation)
		{
			output_file_c<<seaux[i].position<<endl;
			tot_cati++;
		}
		if(seaux[i].nature==anion)
		{
			// output_file_c<<seaux[i].position<<endl;
			tot_anio++;
		}
	}
	cout<<"Total seaux : "<<sizeofbucket<<"; electrons : "<<tot_elec<<"; cations : "<<tot_cati<<"; anions : "<<tot_anio<<endl;
	output_file_e.close();
	output_file_c.close();
}

int main()
{
	cout<<endl;
	PrintUnits();
	std::vector<bucket> buck;
	std::vector<bucket> buck_tmp;
	std::vector<bucket> buck_tot;
	EntryParameters(0);

	int i_temps;
	double nb_ionisation=0.;
	double champ_electrique;
	double rms;
	double mean;
	double distance;

	i_temps=ceil(duree_acqui/T_pas);
	mobilite_electron=1E5*mm*mm/(V*s);
	mobilite_cation=200*mm*mm/(V*s);
	mobilite_anion=200*mm*mm/(V*s);
	recombinaison=1.98E-3*mm*mm*mm/s;
	champ_electrique=HV/gap_size;

	nb_ionisation=Ionisation();

	vvector Nbr_elec_avant=vector_null(N_pas);
	vvector Nbr_elec_apres=vector_null(N_pas);
	vvector Nbr_cati_avant=vector_null(N_pas);
	vvector Nbr_cati_apres=vector_null(N_pas);
	vvector Nbr_anio_avant=vector_null(N_pas);
	vvector Nbr_anio_apres=vector_null(N_pas);
	vvector Rho_charge=vector_null(N_pas);
	
	// double determinant=0.;
	// double **Matrice=(double**)malloc(N_pas*sizeof(double*));
	// double **MatriceInv=(double**)malloc(N_pas*sizeof(double*));
	// for(int i=0;i<N_pas;i++)
	// {
	// 	Matrice[i]=(double*)malloc(N_pas*sizeof(double));
	// 	MatriceInv[i]=(double*)malloc(N_pas*sizeof(double));
	// }

	// for(int i=0;i<N_pas;i++)
	// {
	// 	free(Matrice[i]);  
	// 	free(MatriceInv[i]);  
	// }
	// free(Matrice);
	// free(MatriceInv);
	// cout<<nb_ionisation<<endl;
	// cout<<i_temps<<endl;

	faisceau=true;

	for(int rho_ind=0;rho_ind<N_pas;rho_ind++)	// Pour le premier pas, les nombres de charges - et + s'équilibrent
		Rho_charge->coef[rho_ind]=0.;

	// for(int i=0;i<i_temps;i++)
	for(int i=0;i<11;i++)
	{
		if(i*T_pas>duree_pulse&&faisceau==true)
		{
			cout<<"Pulse fini"<<endl;
			faisceau=false;
		}

		buck_tmp.clear();

		Rho_charge=ElectricField(Rho_charge);	// Calcul du champ

		for(int inter_indice=0;inter_indice<N_pas;inter_indice++)
		{
			if(faisceau==true)	// initialisation en début de pas de temps
			{
				for(int k=0;k<B_pas;k++)
				{
					int indice=inter_indice*B_pas+k;
					bucket seau_ini;
					seau_ini.condition=vivant;
					seau_ini.intervalle=inter_indice+1;
					seau_ini.quanta=nb_ionisation/B_pas;
					seau_ini.position=gap_size/(N_pas*B_pas)*(indice+.5);

					seau_ini.nature=electron;
					seau_ini.charge=-1.;
					buck.push_back(seau_ini);	//electron

					seau_ini.nature=cation;
					seau_ini.charge=1.;
					buck.push_back(seau_ini);	//ion+
				}
			}

			champ_electrique=Rho_charge->coef[inter_indice];
			
			mobilite_electron=ElectronSpeed(champ_electrique);
			mobilite_cation=IonSpeed(cation);
			mobilite_anion=IonSpeed(anion);
			temps_attachement=AttachementTime(champ_electrique);



			// if(buck[buck_indice].quanta<=0.)
			// 	buck[buck_indice].condition=mort;

			

			for(int buck_indice=0;buck_indice<buck.size();buck_indice++)
			{
				if(buck[buck_indice].intervalle==inter_indice+1&&buck[buck_indice].condition==vivant)
				{
					buck[buck_indice].condition=pending;

					switch(buck[buck_indice].nature)
					{
						case electron:
							Nbr_elec_avant->coef[inter_indice]+=buck[buck_indice].quanta;
							mean=T_pas*mobilite_electron*champ_electrique;
							rms=Dispersion(mobilite_electron);
						break;
						case cation:
							Nbr_cati_avant->coef[inter_indice]+=buck[buck_indice].quanta;
							mean=T_pas*mobilite_cation*champ_electrique;
							rms=Dispersion(mobilite_cation);
						break;
						case anion:
							Nbr_anio_avant->coef[inter_indice]+=buck[buck_indice].quanta;
							mean=T_pas*mobilite_anion*champ_electrique;
							rms=Dispersion(mobilite_anion);
						break;
					}
					distance=Gaussian(mean,rms)*buck[buck_indice].charge;
					buck[buck_indice].position+=distance;
					if(buck[buck_indice].position>=gap_size)
					{
						buck[buck_indice].position=gap_size;
						buck[buck_indice].condition=dehors;
					}
					if(buck[buck_indice].position<=0.)
					{
						buck[buck_indice].position=0.;
						buck[buck_indice].condition=dehors;
					}

					int inter_tmp=IntervalLocation(buck[buck_indice].position);
					buck[buck_indice].intervalle=IntervalLocation(buck[buck_indice].position);
					Rho_charge->coef[inter_tmp]+=buck[buck_indice].charge*buck[buck_indice].quanta;

					if(buck[buck_indice].condition==pending)				//Purge des seaux obsoletes
					{
						buck_tmp.push_back(buck[buck_indice]);
						buck[buck_indice].condition=traite;
					}
					else
						buck_tot.push_back(buck[buck_indice]);
				}
			}
		}
		// cout<<buck.size()<<" "<<buck_tmp.size()<<" "<<buck_tot.size()<<endl;
		buck.clear();
		for(int buck_tmp_indice=0;buck_tmp_indice<buck_tmp.size();buck_tmp_indice++)
		{
			buck_tmp[buck_tmp_indice].condition=vivant;
			buck.push_back(buck_tmp[buck_tmp_indice]);
		}
		
		if(i%10==1)
		{
			cout<<endl;
			cout<<"Buck size "<<buck.size()<<endl;
			cout<<"Temps écoulé "<<i*T_pas/ns<<" ns"<<endl;
		}

		if(buck.size()==0)
		{
			cout<<"Tous les seaux ont été vidés"<<endl;
			break;
		}
	}

	for(int buck_tmp_indice=0;buck_tmp_indice<buck_tmp.size();buck_tmp_indice++)
		buck_tot.push_back(buck_tmp[buck_tmp_indice]);
	// for(int buck_indice=0;buck_indice<buck.size();buck_indice++)
	// 	buck_tot.push_back(buck[buck_indice]);
	cout<<"Total seaux "<<buck_tot.size()<<endl;
	// double mean=10.;
	// double rms=1.;
	// for(int i=0;i<1E6;i++)
	// {
	// 	buck.push_back(bucket());
	// 	buck[i].quanta=Gaussian(mean,rms);
	// }
	Python(buck_tot);

	vector_free(Nbr_elec_avant);
	vector_free(Nbr_elec_apres);
	vector_free(Nbr_cati_avant);
	vector_free(Nbr_cati_apres);
	vector_free(Nbr_anio_avant);
	vector_free(Nbr_anio_apres);
	vector_free(Rho_charge);
}