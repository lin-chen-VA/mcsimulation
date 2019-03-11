#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include "Molecule.h"
using namespace std;

class System
{
	private: double system_x, system_y, system_z;
		 double wallmin;
		 double hplanck;
		 double kboltzmann;
		 double temperature;
		 double twoPi;
		 double piSqr;
		 double rSqrtPi;
		 int mvec;
		 double moveStep, rotateStep;
		 int lessParticle;
		 bool haveWall, notOnlyMove;
		 int modGz, modGr;
		 double thickOfWall, factorWall;
		 bool boundaryX, boundaryY, boundaryZ;
		 double thickX, thickY, thickZ;
		 double sigWallX, sigWallY, sigWallZ;
		 double epsWallX, epsWallY, epsWallZ;
		 double rhoX, rhoY, rhoZ;
		 double ewaldTempX, ewaldTempY, ewaldTempZ;
		 int kMax, ksqMax, slide;
		 double kappa;
		 int numpts;
		 int idratio;
		 double chemPotential;
		 vector<Molecule> modelContainer, moleculeContainer;
		 double KbT;
		 vector<int> ihx, ihy, ihz;
		 vector<double> fhold;
		 vector<double> ssum, csum, ssum2, csum2;
		 int checkTime;
		 time_t insertReal, insertReciprocal, insertTotal;
		 time_t moveReal, moveReciprocal, moveTotal;
		 time_t deleteReal, deleteReciprocal, deleteTotal;
		 bool grPlot;
		 int grColumn, grRow, grResolution;
		 double **grArray;
		 int grCount, *grCountNum;
	public:
		 System()
		 {
			 hplanck = 6.62608e-34;
			 kboltzmann = 1.38066e-23;
			 twoPi = 6.2831853;
			 piSqr = 9.869604401;
			 rSqrtPi = 0.564189583;
			 mvec = 1000;
			 string str;
			 ifstream inputMolecule;
			 inputMolecule.open("input.in");
			 ifstream inputParameter;
			 inputParameter.open("parameter.in");
			 ifstream inputModel;
			 inputModel.open("model.in");
			 unsigned seed;
			 seed = time(0);
			 srand(seed);
			 read(inputMolecule, inputModel, inputParameter);
			 if(haveWall)
			 {
				 if(!boundaryX)
				 {
					 wallmin = system_x + thickX;
				 }
				 else
				 {
					 wallmin = system_x;
				 }
				 if(!boundaryY)
				 {
					 if(system_y + thickY < wallmin)
						 wallmin = system_y + thickY;
				 }
				 else
				 {
					 if(system_y < wallmin)
						 wallmin = system_y;
				 }
				 if(!boundaryZ)
				 {
					 if(system_z + thickZ < wallmin)
						 wallmin = system_z + thickZ;
				 }
				 else
				 {
					 if(system_z < wallmin)
						 wallmin = system_z;
				 }
			 }
			 else
			 {
				 wallmin = system_x;
				 if(system_y < wallmin)
					 wallmin = system_y;
				 if(system_z < wallmin)
					 wallmin = system_z;
			 }
			 KbT = 0.00831451/4.184*temperature;
			 ksetup();
			 ewald();
			 inputMolecule.close();
			 inputModel.close();
			 inputParameter.close();
			 if(checkTime == 1)
			 {
				 insertReal = 0;
				 insertReciprocal = 0;
				 insertTotal = 0;
				 moveReal = 0;
				 moveReciprocal = 0;
				 moveTotal = 0;
				 deleteReal = 0;
				 deleteReciprocal = 0;
				 deleteTotal = 0;
			 }
			 if(grPlot)
			 {
				 int tempNum;
				 double systemRange;
				 tempNum = modelContainer.size();
				 grColumn = (1+tempNum)*tempNum/2;
				 systemRange = sqrt(system_x*system_x+system_y*system_y);
				 grRow = int(systemRange*grResolution);
				 grArray = new double *[grRow];
				 for(int i = 0; i < grRow; i++)
					 grArray[i] = new double [grColumn];
				 for(int i = 0; i < grRow; i++)
					 for(int j = 0; j < grColumn; j++)
						 grArray[i][j] = 0.0;
				 grCount = 0;
				 grCountNum = new int [modelContainer.size()];
				 for(int i = 0; i < modelContainer.size(); i++)
					 grCountNum[i] = 0;
			 } 

		 }

		 ~System()
		 {
			 ofstream outputMolecule;
			 outputMolecule.open("output.out");
			 save(outputMolecule);
			 outputMolecule.close();
			 if(grPlot)
			 {
				 for(int i = 0; i < grRow; i++)
					 delete [] grArray[i];
				 delete [] grArray;
			 }
		 }

		 //Mutator
		 void read(ifstream &inputMolecule, ifstream &inputModel, ifstream &inputParameter)
		 {
			 string str;
			 inputParameter>>str>>haveWall;
			 if(str != "haveWall") cout<<"Parameter haveWall error!"<<endl;
			 inputParameter>>str>>system_x;
			 if(str != "system_x") cout<<"Parameter system_x error!"<<endl;
			 inputParameter>>str>>system_y;
			 if(str != "system_y") cout<<"Parameter system_y error!"<<endl;
			 inputParameter>>str>>system_z;
			 if(str != "system_z") cout<<"Parameter system_z error!"<<endl;
			 inputParameter>>str>>temperature;
			 if(str != "temperature") cout<<"Parameter temperature error!"<<endl;
			 inputParameter>>str>>moveStep;
			 if(str != "moveStep") cout<<"Parameter moveStep error!"<<endl;
			 inputParameter>>str>>rotateStep;
			 if(str != "rotateStep") cout<<"Parameter rotateStep error!"<<endl;
			 inputParameter>>str>>lessParticle;
			 if(str != "lessParticle") cout<<"Parameter lessParticle error!"<<endl;
			 inputParameter>>str>>notOnlyMove;
			 if(str != "notOnlyMove") cout<<"Parameter notOnlyMove error!"<<endl;
			 inputParameter>>str>>modGz;
			 if(str != "modGz") cout<<"Parameter modGz error!"<<endl;
			 inputParameter>>str>>modGr;
			 if(str != "modGr") cout<<"Parameter modGr error!"<<endl;
			 inputParameter>>str>>boundaryX;
			 if(str != "boundaryX") cout<<"Parameter boundaryX error!"<<endl;
			 inputParameter>>str>>thickX;
			 if(str != "thickX") cout<<"Parameter thickX error!"<<endl;
			 inputParameter>>str>>rhoX;
			 if(str != "rhoX") cout<<"Parameter rhoX error!"<<endl;
			 inputParameter>>str>>sigWallX;
			 if(str != "sigWallX") cout<<"Parameter sigWallX error!"<<endl;
			 inputParameter>>str>>epsWallX;
			 if(str != "epsWallX") cout<<"Parameter epsWallX error!"<<endl;
			 inputParameter>>str>>boundaryY;
			 if(str != "boundaryY") cout<<"Parameter boundaryY error!"<<endl;
			 inputParameter>>str>>thickY;
			 if(str != "thickY") cout<<"Parameter thickY error!"<<endl;
			 inputParameter>>str>>rhoY;
			 if(str != "rhoY") cout<<"Parameter rhoY error!"<<endl;
			 inputParameter>>str>>sigWallY;
			 if(str != "sigWallY") cout<<"Parameter sigWallY error!"<<endl;
			 inputParameter>>str>>epsWallY;
			 if(str != "epsWallY") cout<<"Parameter epsWallY error!"<<endl;
			 inputParameter>>str>>boundaryZ;
			 if(str != "boundaryZ") cout<<"Parameter boundaryZ error!"<<endl;
			 inputParameter>>str>>thickZ;
			 if(str != "thickZ") cout<<"Parameter thickZ error!"<<endl;
			 inputParameter>>str>>rhoZ;
			 if(str != "rhoZ") cout<<"Parameter rhoZ error!"<<endl;
			 inputParameter>>str>>sigWallZ;
			 if(str != "sigWallZ") cout<<"Parameter sigWallZ error!"<<endl;
			 inputParameter>>str>>epsWallZ;
			 if(str != "epsWallZ") cout<<"Parameter epsWallZ error!"<<endl;
			 inputParameter>>str>>kMax;
			 if(str != "kMax") cout<<"Parameter kMax error!"<<endl;
			 inputParameter>>str>>ksqMax;
			 if(str != "ksqMax") cout<<"Parameter ksqMax error!"<<endl;
			 inputParameter>>str>>kappa;
			 if(str != "kappa") cout<<"Parameter kappa error!"<<endl;
			 inputParameter>>str>>slide;
			 if(str != "slide") cout<<"Parameter slide error!"<<endl;
			 inputParameter>>str>>numpts;
			 if(str != "numpts") cout<<"Parameter numpts error!"<<endl;
			 inputParameter>>str>>idratio;
			 if(str != "idratio") cout<<"Parameter idratio error!"<<endl;
			 inputParameter>>str>>chemPotential;
			 if(str != "chemPotential") cout<<"Parameter chemPotential error!"<<endl;
			 inputParameter>>str>>checkTime;
			 if(str != "checkTime") cout<<"Parameter checkTime error!"<<endl;
			 inputParameter>>str>>grPlot;
			 if(str != "grPlot") cout<<"Parameter grPlot error!"<<endl;
			 inputParameter>>str>>grResolution;
			 if(str != "grResolution") cout<<"Parameter grResolution error!"<<endl;
			 Molecule temp;
			 int end;
			 inputModel.seekg(0L, ios::end);
			 end = inputModel.tellg();
			 inputModel.seekg(0L, ios::beg);
			 wallThick();
			 if(end > 0)
			 {
				 while(inputModel.tellg() < end - 1)
				 {
					 temp.clear();
					 temp.read(inputModel, twoPi, ewaldTempX, ewaldTempY, ewaldTempZ);
					 modelContainer.push_back(temp);
				 }
			 }
			 inputMolecule.seekg(0L, ios::end);
			 end = inputMolecule.tellg();
			 inputMolecule.seekg(0L, ios::beg);
			 if(end > 0)
			 {
				 while(inputMolecule.tellg() < end - 1)
				 {
					 temp.clear();
					 temp.read(inputMolecule, twoPi, ewaldTempX, ewaldTempY, ewaldTempZ);
					 moleculeContainer.push_back(temp);
				 }
			 }
			 cout<<"Program start..."<<endl;
			 cout<<"-----------------------------------------------"<<endl;
			 cout<<"Read "<<modelContainer.size()<<" molecule models from \"model.in\"."<<endl;
			 for(int i = 0; i < modelContainer.size(); i++)
				 cout<<"    "<<"No."<<setw(3)<<i+1<<" "<<modelContainer[i].getName()<<endl;
			 cout<<"Read "<<moleculeContainer.size()<<" molecules from \"input.in\"."<<endl;
			 cout<<"-----------------------------------------------"<<endl;
			 if(notOnlyMove == 1)
			 {
			 cout<<"Then will randomly choose model molecule to insert into system."<<endl;
			 cout<<"System can do movement and deletion operation by randomly choosing\n molecule from molecules in the system."<<endl;
			 }
			 else
			 {
				 cout<<"Then will randomly choose molecule to implement movement."<<endl;
			 }
		 }

		 void save(ofstream &outputMolecule)
		 {
			 for(int i = 0; i < moleculeContainer.size(); i++)
			 {
				 moleculeContainer[i].save(outputMolecule);
			 }
		 }

		 void movement()
		 {
			 if(moleculeContainer.size() == 0)
			 {
				 cout<<"Molecle number in the system is zero, can not do movement operation."<<endl;
				 return;
			 }
			 //cout<<"Molecule number: "<<moleculeContainer.size()<<endl;
			 //cout<<"Test movement."<<endl;
			 time_t start, end, startTotal, endTotal;
			 if(checkTime == 1)
				 time(&startTotal);
			 int choice, last;
			 choice = rand()%moleculeContainer.size();

			 if(moleculeContainer[choice].getFixing()) return;

			 last = moleculeContainer.size() - 1;

			 exchange(choice);//exhange the choosed particle and last particle in the container

			 Molecule temp;
			 temp = moleculeContainer[last];

			 double energyBeforeMove;
			 if(checkTime == 1)
				 time(&start);
			 energyBeforeMove = singleEnergy();
			 if(checkTime == 1)
			 {
				 time(&end);
				 moveReal += end - start;
			 }
			 //cout<<"Before real space: "<<energyBeforeMove<<endl;

			 double reciprocalEnergy = 0;
			 bool isBefore;
			 isBefore = true;
			 if(checkTime == 1)
				 time(&start);
			 ewaldMove(isBefore, reciprocalEnergy);
			 if(checkTime == 1)
			 {
				 time(&end);
				 moveReciprocal += end - start;
			 }
			 //cout<<"Reciprocal energy: "<<reciprocalEnergy<<endl;

			 moleculeContainer[last].shift(moveStep, twoPi, ewaldTempX, ewaldTempY, ewaldTempZ);
			 moleculeContainer[last].tinyRotate(twoPi, rotateStep, ewaldTempX, ewaldTempY, ewaldTempZ);

			 double energyAfterMove;
			 if(checkTime == 1)
				 time(&start);
			 energyAfterMove = singleEnergy();
			 if(checkTime == 1)
			 {
				 time(&end);
				 moveReal += end - start;
			 }
			 //cout<<"After real space: "<<energyAfterMove<<endl;
			 isBefore = false;
			 if(checkTime == 1)
				 time(&start);
			 ewaldMove(isBefore, reciprocalEnergy);
			 if(checkTime == 1)
			 {
				 time(&end);
				 moveReciprocal += end - start;
			 }
			 //cout<<"Reciprocal energy: "<<reciprocalEnergy<<endl;

			 energyAfterMove += reciprocalEnergy;

			 double differenceEnergy;
			 differenceEnergy = energyAfterMove - energyBeforeMove;

			 //cout<<"Difference: "<<differenceEnergy<<endl;

			 if(differenceEnergy < 0 || ran() < exp((-1/KbT)*differenceEnergy))
			 {
				 //cout<<"Movement Accepted "<<differenceEnergy<<endl;
				 if(moleculeContainer[last].outsideTell(system_x, system_y, system_z, haveWall, boundaryX, boundaryY, boundaryZ))
				 {
					 ewaldMoveUnexe(choice);
				 }
				 moleculeContainer[last].align(system_x, system_y, system_z, twoPi, ewaldTempX, ewaldTempY, ewaldTempZ);
				 ewaldMoveExe(choice);
			 }
			 else
			 {
				 moleculeContainer[last] = temp;
				 ewaldMoveUnexe(choice);
			 }
			 if(checkTime == 1)
			 {
				 time(&endTotal);
				 moveTotal += endTotal - startTotal;
			 }

		 }

		 void deletion()
		 {
			 if(moleculeContainer.size() == 0)
			 {
				 cout<<"Molecule number in the system is zero, can not do deletion operation."<<endl;
				 return;
			 }
			 time_t start, end, startTotal, endTotal;
			 if(checkTime == 1)
				 time(&startTotal);
			 //cout<<"Test deletion."<<endl;
			 //cout<<"Molecule number: "<<moleculeContainer.size()<<endl;
			 int choice, last;
			 choice = rand()%moleculeContainer.size();

			 if(moleculeContainer[choice].getFixing()) return;

			 last = moleculeContainer.size() - 1;

			 //cout<<"Atom number: "<<moleculeContainer[choice].getSiteNumber()<<endl;
			 //moleculeContainer[choice].showMolecule();

			 double mass;
			 mass = moleculeContainer[choice].getMass();

			 double thermalWaveLength;
			 thermalWaveLength = sqrt(hplanck*hplanck/twoPi/kboltzmann/temperature/mass*6.02214e+26)*1.0e+10;

			 double probScale;
			 probScale = thermalWaveLength*thermalWaveLength*thermalWaveLength*moleculeContainer.size()/system_x/system_y/system_z;

			 exchange(choice);

			 double deleteEnergy;
			 if(checkTime == 1)
				 time(&start);
			 deleteEnergy = singleEnergy();
			 if(checkTime == 1)
			 {
				 time(&end);
				 deleteReal += end - start;
			 }

			// cout<<"RealSpace: "<<deleteEnergy<<endl;

			 if(checkTime == 1)
				 time(&start);
			 deleteEnergy = deleteEnergy - self(moleculeContainer[last]) + ewaldDelete() - moleculeContainer[last].intra(kappa, wallmin);
			 if(checkTime == 1)
			 {
				 time(&end);
				 deleteReciprocal += end - start;
			 }
			 double acceptFactor;
			 acceptFactor = idratio*probScale*exp(-(chemPotential - deleteEnergy)/KbT);
			 if(ran() < acceptFactor)
			 {
				 ewaldDeleteExe(choice);
				 if(acceptFactor > 1)
					 cout<<"Delete Accepted(energy): "<<deleteEnergy<<" moleclue number: "<<moleculeContainer.size()<<endl;
				 else
					 cout<<"Delete Accepted(random number): "<<deleteEnergy<<" molecule number: "<<moleculeContainer.size()<<endl;
			 }
			 else
			 {
				 ewaldDeleteUnexe(choice);
			 }
			 if(checkTime == 1)
			 {
				 time(&endTotal);
				 deleteTotal += endTotal - startTotal;
			 }
		 }

		 void insertion()
		 {
			 Molecule insertMolecule;
	   		 int choice;
			 time_t start, end, startTotal, endTotal;
			 if(checkTime == 1)
				 time(&startTotal);

			 choice = rand()%modelContainer.size();
			 insertMolecule = modelContainer[choice];

			 //if(insertMolecule.getFixed()) return;

			 double px, py, pz;
			 px = (2*ran() - 1) * system_x/2.0;
			 py = (2*ran() - 1) * system_y/2.0;
			 pz = (2*ran() - 1) * system_z/2.0;
			 
			 insertMolecule.move(px, py, pz, twoPi, ewaldTempX, ewaldTempY, ewaldTempZ);
			 insertMolecule.rotate(twoPi, ewaldTempX, ewaldTempY, ewaldTempZ);
			 moleculeContainer.push_back(insertMolecule);

			 int last;
			 last = moleculeContainer.size() - 1;

			 double mass;
			 mass = insertMolecule.getMass();
			 double thermalWaveLength;
			 thermalWaveLength = sqrt(hplanck*hplanck/twoPi/kboltzmann/temperature/mass*6.02214e+26)*1.0e+10;
			 double probScale;
			 probScale = system_x*system_y*system_z/thermalWaveLength/thermalWaveLength/thermalWaveLength/moleculeContainer.size();

			 double insertEnergy;
			 if(checkTime == 1)
				 time(&start);
			 insertEnergy = singleEnergy();
			 if(checkTime == 1)
			 {
				 time(&end);
				 insertReal += end - start;
			 }

			 //cout<<"Test insertion:"<<endl;
			 //cout<<"Molecule number: "<<moleculeContainer.size()<<endl;
			 //cout<<"Realspace: "<<insertEnergy<<endl;
			 if(checkTime == 1)
				 time(&start);
			 insertEnergy = insertEnergy - self(moleculeContainer[last]) + ewaldInsert() - moleculeContainer[last].intra(kappa, wallmin);
			 if(checkTime == 1)
			 {
				 time(&end);
				 insertReciprocal += end - start;
			 }
			 //cout<<"Insertion energy: "<<insertEnergy<<endl;
			 double acceptFactor;
			 acceptFactor = (1./idratio)*probScale*exp(-(insertEnergy - chemPotential)/KbT);

			 if(ran() < acceptFactor)
			 {
				 ewaldInsertExe();
				 if(acceptFactor > 1)
					 cout<<"Insert Accepted (energy): "<<insertEnergy<<" molecule number: "<<moleculeContainer.size()<<endl;
				 else
					 cout<<"Insert Accepted (random number): "<<insertEnergy<<" molecule number: "<<moleculeContainer.size()<<endl;
			 }
			 else
				 ewaldInsertUnexe();

			 if(checkTime == 1)
			 {
				 time(&endTotal);
				 insertTotal += endTotal - startTotal;
			 }
		 }

		 double self(Molecule m)
		 {
			 double energySelf;

			 energySelf = m.chargeSqr();

			 energySelf = energySelf * kappa / wallmin * rSqrtPi;

			 //cout<<"Test self: "<<energySelf<<endl;
			 return energySelf;
		 }

		 void ksetup()
		 {
			 int xmax, ymax, zmax;
			 xmax = int(kMax * (ewaldTempX/wallmin));
			 ymax = int(kMax * (ewaldTempY/wallmin));
			 zmax = int(kMax * (ewaldTempZ/wallmin));

			 double rksq, rksqmax;
			 rksqmax = double(ksqMax);
			 for(int nx = 0; nx <= xmax; nx++)
				 for(int ny = -ymax; ny <= ymax; ny++)
					 for(int nz = -zmax; nz <= zmax; nz++)
					 {
						 if(nx == 0)
						 {
							 if(ny < 0) continue;
							 if(ny ==0 && nz <= 0) continue;
						 }
						 rksq = sqr(double(nx)/ewaldTempX)+sqr(double(ny)/ewaldTempY)+sqr(double(nz)/ewaldTempZ);
						 if(rksq*sqr(wallmin) <= rksqmax)
						 {
							 if(ihx.size() > mvec) 
								 cout<<"gmothers"<<endl;
							 ihx.push_back(nx);
							 ihy.push_back(ny);
							 ihz.push_back(nz);
							 double temp;
							 temp = 1.0/rksq*exp(-piSqr*rksq/sqr(kappa/wallmin));
							 fhold.push_back(temp);
						 }
					 }
		 }

		 void wallThick()
		 {
			 if(haveWall)
			 {
				 if(!boundaryX)
					 ewaldTempX = system_x + thickX;
				 else
					 ewaldTempX = system_x;
				 if(!boundaryY)
					 ewaldTempY = system_y + thickY;
				 else
					 ewaldTempY = system_y;
				 if(!boundaryZ)
					 ewaldTempZ = system_z + thickZ;
				 else
					 ewaldTempZ = system_z;
			 }
			 else
			 {
				 ewaldTempX = system_x;
				 ewaldTempY = system_y;
				 ewaldTempZ = system_z;
			 }
		 }

		 void ewald()
		 {
			 double ndotr;

			 for(int i = 0; i < ihx.size(); i++)
			 {
				 ssum.push_back(0);
				 csum.push_back(0);
			 }

			 for(int i = 0; i < ihx.size(); i++)
				 for(int j = 0; j < moleculeContainer.size(); j++)
				 {
					 vector<double> c, rxt, ryt, rzt;
					 c = moleculeContainer[j].getChargeContainer();
					 rxt = moleculeContainer[j].getRX();
					 ryt = moleculeContainer[j].getRY();
					 rzt = moleculeContainer[j].getRZ();
					 for(int k = 0; k < c.size(); k++)
					 {
						 ndotr = ihx[i]*rxt[k] + ihy[i]*ryt[k] + ihz[i]*rzt[k];
						 ssum[i] += c[k]*sin(ndotr);
						 csum[i] += c[k]*cos(ndotr);
					 }
				 }

			 for(int i = 0; i < ihx.size(); i++)
			 {
				 ssum2.push_back(ssum[i]);
				 csum2.push_back(csum[i]);
			 }
		 }

		 double ewaldDelete()
		 {
			 int last;
			 last = moleculeContainer.size() - 1;

			 vector<double> s, c;
			 moleculeContainer[last].singleEwald(s, c, ihx, ihy, ihz);

			 for(int i = 0; i < ihx.size(); i++)
			 {
				 ssum2[i] -= s[i];
				 csum2[i] -= c[i];
			 }

			 double energyEwaldDelete = 0;

			 for(int i = 0; i < ihx.size(); i++)
			 {
				 energyEwaldDelete += fhold[i]*(csum[i]*csum[i] + ssum[i]*ssum[i] - csum2[i]*csum2[i] - ssum2[i]*ssum2[i]);
			 }

			 energyEwaldDelete = 2.0*energyEwaldDelete/(twoPi*ewaldTempX*ewaldTempY*ewaldTempZ);
			 return energyEwaldDelete;
		 }

		 double ewaldInsert()
		 {
			 int last;
			 last = moleculeContainer.size() - 1;

			 vector<double> s, c;
			 moleculeContainer[last].singleEwald(s, c, ihx, ihy, ihz);
			 for(int i = 0; i < ihx.size(); i++)
			 {
				 ssum2[i] += s[i];
				 csum2[i] += c[i];
			 }

			 double energyEwaldInsert = 0;

			 for(int i = 0; i < ihx.size(); i++)
			 {
				 energyEwaldInsert += fhold[i]*(csum2[i]*csum2[i] + ssum2[i]*ssum2[i] - csum[i]*csum[i] - ssum[i]*ssum[i]);
			 }

			// cout<<"ewald thick: "<<ewaldTempX<<" "<<ewaldTempY<<" "<<ewaldTempZ<<" "<<wallmin<<endl;
			 energyEwaldInsert = 2.0*energyEwaldInsert/(twoPi*ewaldTempX*ewaldTempY*ewaldTempZ);
			 //cout<<"Test ewaldInsert: "<<energyEwaldInsert<<endl;
			 return energyEwaldInsert;
		 }

		 void ewaldMove(bool isBefore, double &energyEwaldMove)
		 {
			 int last;
			 last = moleculeContainer.size() - 1;

			 energyEwaldMove = 0;

			 vector<double> s, c;
			 
			 if(isBefore == true)
			 {
				 moleculeContainer[last].singleEwald(s, c, ihx, ihy, ihz);
				 for(int i = 0; i < ihx.size(); i++)
				 {
					 ssum2[i] -= s[i];
					 csum2[i] -= c[i];
				 }
			 }
			 else
			 {
				 moleculeContainer[last].singleEwald(s, c, ihx, ihy, ihz);
				 for(int i = 0; i < ihx.size(); i++)
				 {
					 ssum2[i] += s[i];
					 csum2[i] += c[i];
				 }
				 
				 for(int i = 0; i < ihx.size(); i++)
				 {
					 energyEwaldMove += fhold[i]*(csum2[i]*csum2[i] + ssum2[i]*ssum2[i] - csum[i]*csum[i] - ssum[i]*ssum[i]);
				 }

				 energyEwaldMove = 2.0*energyEwaldMove/(twoPi*ewaldTempX*ewaldTempY*ewaldTempZ);
			 }
		 }

		 void ewaldMoveExe(int choice)
		 {
			 sinCosSumRenew();
			 exchange(choice);
		 }

		 void ewaldMoveUnexe(int choice)
		 {
			 sinCosSumRestore();
			 exchange(choice);
		 }

		 void ewaldDeleteExe(int choice)
		 {
			 sinCosSumRenew();
			 exchange(choice);
			 moleculeContainer.erase(moleculeContainer.begin() + choice);

		 }

		 void ewaldInsertExe()
		 {
			 sinCosSumRenew();
		 }

		 void ewaldInsertUnexe()
		 {
			 sinCosSumRestore();
			 moleculeContainer.pop_back();
		 }


		 void ewaldDeleteUnexe(int choice)
		 {
			 sinCosSumRestore();
			 exchange(choice);
		 }

		 double ewaldFull()
		 {
			 double ewaldFullTemp;
			 ewaldFullTemp = 0.0;

			 for(int i = 0; i < ihx.size(); i++)
				 ewaldFullTemp += fhold[i]*(csum[i]*csum[i]+ssum[i]*ssum[i]);
			 ewaldFullTemp = 2.0*ewaldFullTemp/(twoPi*ewaldTempX*ewaldTempY*ewaldTempZ);

			 return ewaldFullTemp;
		 }

		 double energyFull()
		 {
			 double energyFullTemp;
			 energyFullTemp = 0.0;
			 int last;
			 last = moleculeContainer.size() - 1;

			 for(int i = 0; i < moleculeContainer.size(); i++)
			 {
				 exchange(i);
				 energyFullTemp += singleEnergy() - self(moleculeContainer[last]) - moleculeContainer[last].intra(kappa, wallmin);
				 exchange(i);
			 }

			 energyFullTemp += ewaldFull();
			 
			 return energyFullTemp/2.0;
		 }

		 double averageFullEnergy()
		 {
			 return energyFull()/moleculeContainer.size();
		 }

		 void sinCosSumRestore()
		 {
			 for(int i = 0; i < ssum.size(); i++)
			 {
				 ssum2[i] = ssum[i];
				 csum2[i] = csum[i];
			 }
		 }

		 void sinCosSumRenew()
		 {
			 for(int i = 0; i < ssum.size(); i++)
			 {
				 ssum[i] = ssum2[i];
				 csum[i] = csum2[i];
			 }
		 }

		 double singleEnergy()
		 {
			 int last;
			 last = moleculeContainer.size() - 1;
			 double energy = 0;
			 for(int i = 0; i < moleculeContainer.size() - 1; i++)
			 {
				 energy += moleculeContainer[last].interEnergy(moleculeContainer[i], kappa, wallmin, system_x, system_y, system_z, haveWall, boundaryX, boundaryY, boundaryZ);
			 }
			 //cout<<"single energy: "<<energy<<endl;
			 //cout<<"system "<<system_x<<" "<<system_y<<" "<<system_z<<endl;
			 //cout<<"ewald "<<ewaldTempX<<" "<<ewaldTempY<<" "<<ewaldTempZ<<endl;
			 if(haveWall)
				 energy += moleculeContainer[last].wallInteraction(system_x, system_y, system_z, twoPi, boundaryX, thickX, boundaryY, thickY, boundaryZ, thickZ, sigWallX, epsWallX, sigWallY, epsWallY, sigWallZ, epsWallZ, rhoX, rhoY, rhoZ);
			 //cout<<"Wall energy: "<<temp<<endl;
                         return energy;
		 }			 

		 double realspace(double r)
		 {
			 double a1 = 0.254829592;
			 double a2 = -0.284496736;
			 double a3 = 1.421413741;
			 double a4 = -1.453152027;
			 double a5 = 1.061405429;
			 double pew = 0.3275911;

			 double t, xsq, tp, real, rkap, erfc;

			 rkap = (kappa/wallmin)*r;
			 t = 1.0/(1.0 + pew * rkap);
			 xsq = rkap * rkap;
			 tp = t * ( a1 + t * ( a2 + t * ( a3 + t * ( a4 + t * a5 ) ) ) );
			 erfc = tp*exp(-xsq);
			 real = erfc/r;

			 return real;
		 }


		 void exchange(int i)//exchange choosed particle with last particle
		 {
			 Molecule temp;
			 int last;
			 last = moleculeContainer.size() - 1;
			 temp = moleculeContainer[i];
			 moleculeContainer[i] = moleculeContainer[last];
			 moleculeContainer[last] = temp;
		 }

		 //return a number within [0,1)
		 double ran()
		 {
			 return double(rand())/(double(RAND_MAX) + double(1));
		 }

		 int sqr(int n) {return n*n;}

		 double sqr(double n) {return n*n;}

		 //analysis the time spend of each part
		 void timeAnalysis()
		 {
			 cout<<"Insertion part spend "<<insertTotal<<" seconds:"<<endl;
			 cout<<"    Real part: "<<insertReal<<" seconds, "<<double(insertReal)/insertTotal*100<<"%"<<endl;
			 cout<<"    Reciprocal part: "<<insertReciprocal<<" seconds, "<<double(insertReciprocal)/insertTotal*100<<"%"<<endl;
			 cout<<"Movement part spend "<<moveTotal<<" seconds:"<<endl;
			 cout<<"    Real part: "<<moveReal<<" seconds, "<<double(moveReal)/moveTotal*100<<"%"<<endl;
			 cout<<"    Reciprocal part: "<<moveReciprocal<<" seconds, "<<double(moveReciprocal)/moveTotal*100<<"%"<<endl;
			 cout<<"Deletion part spend "<<deleteTotal<<" seconds:"<<endl;
			 cout<<"    Real part: "<<deleteReal<<" seconds, "<<double(deleteReal)/deleteTotal*100<<"%"<<endl;
			 cout<<"    Reciprocal part: "<<deleteReciprocal<<" seconds, "<<double(deleteReciprocal)/deleteTotal*100<<"%"<<endl;
			 cout<<"Totally, real part spend "<<double(insertReal + moveReal + deleteReal)/(insertTotal + moveTotal + deleteTotal)*100<<"%, reciprocal part spend "<<double(insertReciprocal + moveReciprocal + deleteReciprocal)/(insertTotal + moveTotal + deleteTotal)*100<<"%."<<endl;
		 }

		 //plot the molecule number vs cycle for each molecule
		 void moleculeNumberPlot(int cycle, ofstream &file)
		 {
			 //cout<<"model number: "<<modelContainer.size()<<endl;
			 
			 int *ptr;
			 ptr = new int [modelContainer.size()];

			 for(int i = 0; i < modelContainer.size(); i++)
				 ptr[i] = 0;

			 for(int i = 0; i < modelContainer.size(); i++)
				 for(int j = 0; j < moleculeContainer.size(); j++)
					 if(modelContainer[i].getName() == moleculeContainer[j].getName())
						 ptr[i]++;

			 file<<setw(10)<<cycle;
			 for(int i = 0; i < modelContainer.size(); i++)
				 file<<setw(10)<<ptr[i];
			 file<<endl;

			 delete [] ptr;
		 }

		 void gr(ofstream &file)
		 {
			 int rowNum, columnNum;
			 int frontNum, backNum;
			 static int countGr = 0;
			 //if the system include three molecule, 0, 1, 2, 3, the column is 00, 01, 02, 03, 11, 12, 13, 22, 23, 33
			 for(int i = 0; i < moleculeContainer.size()-1; i++)
				 for(int j = i+1; j < moleculeContainer.size(); j++)
				 {
					 //find column number
					 for(int k = 0; k < modelContainer.size(); k++)
					 {
						 if(moleculeContainer[i].getName() == modelContainer[k].getName())
							 frontNum = k;
						 if(moleculeContainer[j].getName() == modelContainer[k].getName())
							 backNum = k;
					 }
					 columnNum = 0;
					 for(int k = modelContainer.size(); k > modelContainer.size() - frontNum; k--)
					 {
						 columnNum += k;
					 }
					 columnNum = columnNum + backNum - frontNum;
					 //find row number
					 rowNum =int(distance(moleculeContainer[i], moleculeContainer[j], system_x, system_y, system_z)*grResolution);
					 //renew gr array
					 grArray[rowNum][columnNum] += 2.0;
				 }
			 //countGr += moleculeContainer.size();
			 grCount++;
			 for(int i = 0; i < modelContainer.size(); i++)
				 for(int j = 0; j < moleculeContainer.size(); j++)
					 if(modelContainer[i].getName() == moleculeContainer[j].getName())
						 grCountNum[i]++;
			 double deltV, deltD, radii;
			 double tmp; 
			 for(int i = 0; i < grRow; i++)
			 {
				 radii = double(i)/grResolution;
				 file<<setprecision(3)<<fixed<<setw(5)<<radii;
				 for(int j = 0; j < grColumn; j++)
				 {
					 deltD = double(1)/grResolution;
					 deltV = 4*3.1415926*radii*radii*deltD;
					 //tmp = deltV*grCount*moleculeContainer.size()*moleculeContainer.size()/system_x/system_y/system_z;
					 tmp = deltV*grCount*grCountNum[0]/grCount*grCountNum[0]/grCount/system_x/system_y/system_z;
					 file<<setprecision(4)<<fixed<<" "<<setw(8)<<grArray[i][j]/tmp;
				 }
				 file<<endl;
			 }

		 }


		 double distance(Molecule m1, Molecule m2, double rangeX, double rangeY, double rangeZ)
		 {
			 return m1.centralDistance(m2, rangeX, rangeY, rangeZ);
		 }

		 void generatePdb(){
			 ofstream file;
			 file.open("output.pdb");
			 for(int i = 0; i < moleculeContainer.size(); i++)
			 {
				 moleculeContainer[i].moleculePdb(file, i);
			 }
			 file.close();
		 }

		 //Accessor
		 bool getHaveWall() {return haveWall;}
		 double getSystemX() {return system_x;}
		 double getSystemY() {return system_y;}
		 double getSystemZ() {return system_z;}
		 double getWallmin() {return wallmin;}
		 double getTemperature() {return temperature;}
		 double getMoveStep() {return moveStep;}
		 double getRotateStep() {return rotateStep;}
		 int getLessParticle() {return lessParticle;}
		 bool getNotOnlyMove() {return notOnlyMove;}
		 int getModGz() {return modGz;}
		 int getModGr() {return modGr;}
		 double getThickOfWall() {return thickOfWall;}
		 double getFactorWall() {return factorWall;}
		 bool getHaveBoundaryX() {return boundaryX;}
		 bool getHaveBoundaryY() {return boundaryY;}
		 bool getHaveBoundaryZ() {return boundaryZ;}
		 int getkMax() {return kMax;}
		 int getksqMax() {return ksqMax;}
		 double getKappa() {return kappa;}
		 int getSlide() {return slide;}
		 int getMoleculeNumber() {return moleculeContainer.size();}
		 int getModelNumber() {return modelContainer.size();}
		 int getIdratio() {return idratio;}
};
