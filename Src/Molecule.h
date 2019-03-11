#ifndef MOLECULE_H
#define MOLECULE_H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include "Atom.h"
using namespace std;

class Molecule
{
	private: string name;
		 int siteNumber;
		 double mass;
		 vector<Atom> atomContainer;
		 vector<double> chargeContainer, rx, ry, rz;
		 bool fixing;
	
	public:
		 Molecule()
		 {
		 }

		 Molecule(int n, ifstream &inputFile, double twopi, double system_x, double system_y, double system_z)
		 {
			 siteNumber = n;
			 read(inputFile,twopi, system_x, system_y, system_z);
		 }

		 ~Molecule()
		 {
		 }

		 //Mutator
		 void setSiteNumber(int n)
		 {
			 siteNumber = n;
		 }

		 void setName(string n)
		 {
			 name = n;
		 }

		 void save(ofstream &file)
		 {
			 file<<fixed<<setprecision(4)<<setw(10)<<name<<setw(10)<<siteNumber<<setw(10)<<mass<<setw(10)<<fixing<<endl;
			 for(int i = 0; i < siteNumber; i++)
			 {
				 atomContainer[i].save(file);
			 }
		 }

		 void read(ifstream &file, double twopi, double wallx, double wally, double wallz)
		 {
			 Atom atom;
			 file>>name>>siteNumber>>mass>>fixing;
			 //cout<<name<<" "<<siteNumber<<" "<<mass<<" "<<fixing<<endl;
			 for(int i = 0; i < siteNumber; i++)
			 {
				 atom.read(file);
				 atomContainer.push_back(atom);
			 }

			 if(atomContainer.size() != siteNumber)
			 {
				 cout<<"Error is found in reading molecule!"<<endl;
			 }

			 setR(twopi, wallx, wally, wallz);
		 }

		 void setR(double twopi, double wallx, double wally, double wallz)
		 {
			 chargeContainer.clear();
			 rx.clear();
			 ry.clear();
			 rz.clear();
			 double tempx, tempy, tempz;
			 for(int i = 0; i < siteNumber; i++)
			 {
				 tempx = twopi/wallx*atomContainer[i].getX();
				 tempy = twopi/wallx*atomContainer[i].getY();
				 tempz = twopi/wally*atomContainer[i].getZ();
				 chargeContainer.push_back(atomContainer[i].getCharge());
				 rx.push_back(tempx);
				 ry.push_back(tempy);
				 rz.push_back(tempz);
			 }
		 }


		 void operator= (const Molecule &right)
		 {
			 name = right.name;
			 siteNumber = right.siteNumber;
			 mass = right.mass;
			 fixing = right.fixing;
			 atomContainer = right.atomContainer;
			 chargeContainer = right.chargeContainer;
			 rx = right.rx;
			 ry = right.ry;
			 rz = right.rz;
		 }

		 void shift(double step, double twopi, double wallx, double wally, double wallz)
		 {
			 double dx, dy, dz;
			 
			 dx = (ran() - 0.5)*step;
			 dy = (ran() - 0.5)*step;
			 dz = (ran() - 0.5)*step;

			 for(int i = 0; i < atomContainer.size(); i++)
			 {
				 atomContainer[i].shift(dx, dy, dz);
			 }

			 setR(twopi, wallx, wally, wallz);
		 }

		 void move(double dx, double dy, double dz, double twopi, double wallx, double wally, double wallz)
		 {
			 for(int i = 0; i < atomContainer.size(); i++)
			 {
				 atomContainer[i].shift(dx, dy, dz);
			 }

			 setR(twopi, wallx, wally, wallz);
		 }

		 void tinyRotate(double twoPi, double step, double wallx, double wally, double wallz) //tiny quaternion rotate
		 {
			 double gammar, betar, alphar;
			 Atom *temp, *temp2;
			 temp = new Atom[atomContainer.size()];
			 temp2 = new Atom[atomContainer.size()];

			 double **rotate;
			 rotate = new double *[3];

			 for(int i = 0; i < 3; i++)
				 rotate[i] = new double[3];

			 gammar = twoPi*(ran()-0.5)*step;
			 betar = twoPi*(ran()-0.5)*step;
			 alphar = twoPi*(ran()-0.5)*step;
			 rotate[0][0] = cos(gammar)*cos(betar)*cos(alphar) - sin(gammar)*sin(alphar);
			 rotate[0][1] = cos(gammar)*cos(betar)*sin(alphar) + sin(gammar)*cos(alphar);
			 rotate[0][2] = -cos(gammar)*sin(betar);
			 rotate[1][0] = -sin(gammar)*cos(betar)*cos(alphar) - cos(gammar)*sin(alphar);
			 rotate[1][1] = -sin(gammar)*cos(betar)*sin(alphar) + cos(gammar)*cos(alphar);
			 rotate[1][2] = sin(gammar)*sin(betar);
			 rotate[2][0] = sin(betar)*cos(alphar);
			 rotate[2][1] = sin(betar)*sin(alphar);
			 rotate[2][2] = cos(betar);


			 for(int i = 1; i < atomContainer.size(); i++)
				 temp[i] = atomContainer[i] - atomContainer[0];

			 for(int i = 1; i < atomContainer.size(); i++)
				 temp2[i] = temp[i].rotate(rotate, 3);

			 for(int i = 1; i < atomContainer.size(); i++)
				 atomContainer[i] = temp2[i] + atomContainer[0];

			 for(int i = 0; i < 3; i++)
				 delete [] rotate[i];

			 delete [] rotate;
			 delete [] temp;
			 delete [] temp2;

			 setR(twoPi, wallx, wally, wallz);
		 }

		 void rotate(double twopi, double wallx, double wally, double wallz)//random quaternion rotate
		 {
			 double chi, mu, S1, iota, zeta, S2;
			 Atom *temp, *temp2;
			 temp = new Atom[atomContainer.size()];
			 temp2 = new Atom[atomContainer.size()];
			 
			 double **rotate;	
	       		 rotate = new double *[3];
			 for(int i = 0; i < 3; i++)
				 rotate[i] = new double[3];

			 chi = 2 * ran() - 1;
			 mu = 2 * ran() - 1;
			 S1 = chi * chi + mu * mu;
			 while(S1 >= 1)
			 {
				 chi = 2 * ran() - 1;
				 mu = 2 * ran() - 1;
				 S1 = chi * chi + mu * mu;
			 }

			 iota = 2 * ran() - 1;
			 zeta = 2 * ran() - 1;
			 S2 = iota * iota + zeta * zeta;
			 while(S2 >= 1)
			 {
				 iota = 2 * ran() - 1;
				 zeta = 2 * ran() - 1;
				 S2 = iota * iota + zeta * zeta;
			 }

			 iota = iota * sqrt((1-S1)/S2);
			 zeta = zeta * sqrt((1-S1)/S2);

			 rotate[0][0] = chi*chi - mu*mu - iota*iota + zeta*zeta;
			 rotate[0][1] = 2.0*((chi*mu) - (iota*zeta));
			 rotate[0][2] = 2.0*((iota*chi) + (mu*zeta));
			 rotate[1][0] = 2.0*((chi*mu) + (iota*zeta));
			 rotate[1][1] = mu*mu - iota*iota - chi*chi + zeta*zeta;
			 rotate[1][2] = 2.0*((mu*iota) - (chi*zeta));
			 rotate[2][0] = 2.0*((iota*chi) - (mu*zeta));
			 rotate[2][1] = 2.0*((mu*iota) + (chi*zeta));
			 rotate[2][2] = iota*iota - chi*chi - mu*mu + zeta*zeta;

			 for(int i = 1; i < atomContainer.size(); i++)
				 temp[i] = atomContainer[i] - atomContainer[0];

			 for(int i = 1; i < atomContainer.size(); i++)
				 temp2[i] = temp[i].rotate(rotate, 3);
			
			 for(int i = 1; i < atomContainer.size(); i++)
				 atomContainer[i] = temp2[i] + atomContainer[0];

			 for(int i = 0; i < 3; i++)
				 delete [] rotate[i];

			 delete [] rotate;
			 delete [] temp;
			 delete [] temp2;

			 setR(twopi, wallx, wally, wallz);
		 }

		 void align(double bx, double by, double bz, double twopi, double wallx, double wally, double wallz)
		 {
			 int xIsOver = 0, yIsOver = 0, zIsOver = 0;

			 if(atomContainer[0].getX() > bx/2.0) xIsOver = 1;
			 if(atomContainer[0].getX() < -bx/2.0) xIsOver = -1;

			 if(atomContainer[0].getY() > by/2.0) yIsOver = 1;
			 if(atomContainer[0].getY() < -by/2.0) yIsOver = -1;

			 if(atomContainer[0].getZ() > bz/2.0) zIsOver = 1;
			 if(atomContainer[0].getZ() < -bz/2.0) zIsOver = -1;

			 for(int i = 0; i < atomContainer.size(); i++)
				 atomContainer[i].align(bx,xIsOver,by,yIsOver,bz,zIsOver);
			 
			 setR(twopi, wallx, wally, wallz);
		 }

		 bool outsideTell(double bx, double by, double bz, bool wall, bool tx, bool ty, bool tz)
		 {
			 if(wall)
			 {
			 if(!tx)
			 {
				 if(fabs(atomContainer[0].getX()) > bx/2.0)
					 return true;
			 }
			 if(!ty)
			 {
				 if(fabs(atomContainer[0].getY()) > by/2.0)
					 return true;
			 }
			 if(!tz)
			 {
				 if(fabs(atomContainer[0].getZ()) > bz/2.0)
					 return true;
			 }
			 return false;
			 }
			 else
				 return false;
		 }

		 double chargeSqr()
		 {
			 double charge, chargesqr = 0;

			 for(int i = 0; i < atomContainer.size(); i++)
			 {
				 charge = atomContainer[i].getCharge();
				 chargesqr += charge*charge;
			 }

			 return chargesqr;
		 }

		 void singleEwald(vector<double> &s, vector<double> &c, vector<int> ihx, vector<int> ihy, vector<int> ihz)
		 {
			 double dot;
			 s.clear();
			 c.clear();
			 for(int i = 0; i < ihx.size(); i ++)
			 {
				 s.push_back(0);
				 c.push_back(0);
			 }
			 for(int i = 0; i < siteNumber; i++)
			 {
				 for(int j = 0; j < ihx.size(); j++)
				 {
					 dot =  ihx[j]*rx[i] + ihy[j]*ry[i] + ihz[j]*rz[i];
					 s[j] += chargeContainer[i]*sin(dot);
					 c[j] += chargeContainer[i]*cos(dot);
				 }
			 }
		 }

		 double interEnergy(Molecule &m, double kappa, double wallmin, double wallx, double wally, double wallz, bool haveWall, bool boundaryX, bool boundaryY, bool boundaryZ)
		 {
			 double d;
			 Atom atomA, atomB;
			 double energy = 0;
			 double tempPotential;

			 //cout<<"Molecule a: "<<atomContainer.size()<<endl;
			 //cout<<"Molecule b: "<<m.atomContainer.size()<<endl;
			 for(int i = 0; i < atomContainer.size(); i++)
			 {
				 atomA = atomContainer[i];
				 for(int j = 0; j < m.atomContainer.size(); j++)
				 {
					 atomB = m.atomContainer[j];
					 d = interAtomDistance(atomA, atomB, wallx, wally, wallz, haveWall, boundaryX, boundaryY, boundaryZ);
					 //cout<<"Distance: "<<d<<endl;
					 if((atomA.getEps() == 0) || (atomB.getEps() == 0))
						 tempPotential = 0;
					 else
						 tempPotential = potentialLJ(d, atomA, atomB);
					 energy += tempPotential + atomA.getCharge()*atomB.getCharge()*realspace(d, kappa, wallmin);
				 }
			 }
			return energy;
		 }

		 double potentialLJ(double r, Atom a, Atom b)
		 {
			 double sigma, epsilon, energy;

			 sigma = (a.getSig() + b.getSig())/2.0;
			 epsilon = sqrt(a.getEps()*b.getEps());

			 double ep4, exp1;
			 ep4 = 4.0*epsilon;
			 exp1 = sixTime(sigma/r);
			 energy = ep4*exp1*(exp1 - 1.0);

			 return energy;
		 }

	         double sixTime(double x)
		 {
			 return x*x*x*x*x*x;
		 }

		 double interAtomDistance(Atom a, Atom b, double wallx, double wally, double wallz, bool haveWall, bool boundaryX, bool boundaryY, bool boundaryZ)
		 {
			 double temp, tempx, tempy, tempz;

			 if(haveWall)
			 {
				 if(!boundaryX)
				 {
					 tempx = a.getX() - b.getX();
					 tempx = tempx*tempx;
				 }
				 else
				 {
					 tempx = a.getX() - b.getX();
					 tempx = tempx - wallx*(int(tempx*2/wallx)-int(tempx/wallx));
					 tempx = tempx*tempx;
				 }
				 if(!boundaryY)
				 {
					 tempy = a.getY() - b.getY();
					 tempy = tempy*tempy;
				 }
				 else
				 {
					 tempy = a.getY() - b.getY();
					 tempy = tempy - wally*(int(tempy*2/wally)-int(tempy/wally));
					 tempy = tempy*tempy;
				 }
				 if(!boundaryZ)
				 {
					 tempz = a.getZ() - b.getZ();
					 tempz = tempz*tempz;
				 }
				 else
				 {
					 tempz = a.getZ() - b.getZ();
					 tempz = tempz - wallz*(int(tempz*2/wallz)-int(tempz/wallz));
					 tempz = tempz*tempz;
				 }
			 }
			 else
			 {
				 tempx = a.getX() - b.getX();
				 tempx = tempx - wallx*(int(tempx*2/wallx)-int(tempx/wallx));
				 tempx = tempx*tempx;
				 tempy = a.getY() - b.getY();
				 tempy = tempy - wally*(int(tempy*2/wally)-int(tempy/wally));
				 tempy = tempy*tempy;
				 tempz = a.getZ() - b.getZ();
				 tempz = tempz - wallz*(int(tempz*2/wallz)-int(tempz/wallz));
				 tempz = tempz*tempz;
			 }

			 temp = sqrt(tempx + tempy + tempz);

			 return temp;
		 }

		 double distance(int i, int j)
		 {
			 double dx, dy, dz;

			 dx = atomContainer[i].getX() - atomContainer[j].getX();
			 dy = atomContainer[i].getY() - atomContainer[j].getY();
			 dz = atomContainer[i].getZ() - atomContainer[j].getZ();

			 double d;

			 d = sqrt(dx*dx + dy*dy + dz*dz);

			 return d;
		 }

		 double intra(double kappa, double wallmin)
		 {
			 double energyIntra = 0;

			 for(int i = 0; i < atomContainer.size()-1; i++)
			 {
				 for(int j = i+1; j< atomContainer.size(); j++)
				 {
					 double d;
					 d = distance(i,j);
					 //cout<<"Atom number: "<<atomContainer.size()<<endl;
					 //cout<<"Distance: "<<d<<" "<<i<<" "<<j<<endl;
					 energyIntra += atomContainer[i].getCharge()*atomContainer[j].getCharge()/d - atomContainer[i].getCharge()*atomContainer[j].getCharge()*realspace(d, kappa, wallmin);
				 }
			 }

			 //cout<<"Test intra: "<<energyIntra<<endl;
			 return energyIntra;
		 }

		 double realspace(double r, double kappa, double wallmin)
		 {
			 double real;
			 double a1 = 0.254829592;
			 double a2 = -0.284496736;
			 double a3 = 1.421413741;
			 double a4 = -1.453152027;
			 double a5 = 1.061405429;
			 double pew = 0.3275911;
			 double t, xsq, tp, rkap, erfc;					
			 rkap = (kappa/wallmin)*r;							
			 t = 1.0/(1.0 + pew * rkap);
			 xsq = rkap * rkap;
			 tp = t * ( a1 + t * ( a2 + t * ( a3 + t * ( a4 + t * a5 ) ) ) );
			 erfc = tp*exp(-xsq);
			 real = erfc/r;								
			 return real;
		 }

		 void showAtom(Atom m)
		 {
			 cout<<"Atom: "<<m.getName()<<endl;
			 cout<<m.getX()<<endl<<m.getY()<<endl<<m.getZ()<<endl;
			 cout<<m.getCharge()<<endl<<m.getSig()<<endl<<m.getEps()<<endl;
		 }

		 void showMolecule()
		 {
			 cout<<"Name: "<<name<<endl;
			 cout<<"site: "<<siteNumber<<endl;
			 cout<<"Mass: "<<mass<<endl;
			 cout<<"Atom number: "<<atomContainer.size()<<endl;
			 cout<<"Charge: "<<endl;
			 for(int i = 0; i < chargeContainer.size(); i++)
				 cout<<chargeContainer[i]<<" ";
			 cout<<endl;
			 cout<<"rx       ry       rz"<<endl;
			 for(int i = 0; i < rx.size(); i++)
				 cout<<rx[i]<<" "<<ry[i]<<" "<<rz[i]<<endl;
		 }

		 double ran()
		 {
			 return double(rand())/(double(RAND_MAX) + double(1));
		 }

		 void clear()
		 {
			 atomContainer.clear();
			 rx.clear();
			 ry.clear();
			 rz.clear();
			 chargeContainer.clear();
		 }

		 double wallInteraction(double wallx, double wally, double wallz, double twoPi, bool boundaryX, double thickx, bool boundaryY, double thicky, bool boundaryZ, double thickz, double wallSigX, double wallEpsX, double wallSigY, double wallEpsY, double wallSigZ, double wallEpsZ, double rhox, double rhoy, double rhoz)
		 {
			 double energyWall;
			 energyWall = 0.0;

			 double sigFS, epsFS;
			 double factor;
			 double tempWallEnergy;

			 for(int i = 0; i < atomContainer.size(); i++)
			 {
				 if(atomContainer[i].getEps() == 0)
					 continue;
				 if(!boundaryX)
				 {
					 //cout<<"Boundary x"<<endl;
					 sigFS = (atomContainer[i].getSig() + wallSigX)/2.0;
					 epsFS = sqrt(atomContainer[i].getEps() * wallEpsX)/2.0;
					 factor = twoPi*rhox*epsFS*sigFS*sigFS*thickx;
					 double tempX;
					 tempX = wallx/2.0 - fabs(atomContainer[i].getX());
					 tempWallEnergy = factor*(0.4*pow(sigFS/tempX, 10.0) - pow(sigFS/tempX, 4.0) - pow(sigFS, 4.0)/3/thickx/pow(0.61*thickx+tempX, 3.0));
					 tempX = fabs(wallx - tempX);
					 tempWallEnergy += factor*(0.4*pow(sigFS/tempX, 10.0) - pow(sigFS/tempX, 4.0) - pow(sigFS, 4.0)/3/thickx/pow(0.61*thickx+tempX, 3.0));
					 energyWall += tempWallEnergy;
				 }
				 if(!boundaryY)
				 {
					 //cout<<"Boundary y"<<endl;
					 sigFS = (atomContainer[i].getSig() + wallSigY)/2.0;
					 epsFS = sqrt(atomContainer[i].getEps() * wallEpsY);
					 factor = twoPi*rhoy*epsFS*sigFS*sigFS*thicky;
					 double tempY;
					 tempY = wally/2.0 - fabs(atomContainer[i].getY());
					 tempWallEnergy = factor*(0.4*pow(sigFS/tempY, 10.0) - pow(sigFS/tempY, 4.0) - pow(sigFS, 4.0)/3/thicky/pow(0.61*thicky+tempY, 3.0));
					 tempY = fabs(wally - tempY);
					 tempWallEnergy += factor*(0.4*pow(sigFS/tempY, 10.0) - pow(sigFS/tempY, 4.0) - pow(sigFS, 4.0)/3/thicky/pow(0.61*thicky+tempY, 3.0));
					 energyWall += tempWallEnergy;
				 }
				 if(!boundaryZ)
				 {
					 //cout<<"Boundary z"<<endl;
					 //cout<<atomContainer[i].getName()<<endl;
					 sigFS = (atomContainer[i].getSig() + wallSigZ)/2.0;
					 epsFS = sqrt(atomContainer[i].getEps() * wallEpsZ);
					 factor = twoPi*rhoz*epsFS*sigFS*sigFS*thickz;
					 double tempZ;
					 tempZ = wallz/2.0 - fabs(atomContainer[i].getZ());
					 //cout<<"Z+: "<<tempZ;
					 tempWallEnergy = factor*(0.4*pow(sigFS/tempZ, 10.0) - pow(sigFS/tempZ, 4.0) - pow(sigFS, 4.0)/3/thickz/pow(0.61*thickz+tempZ, 3.0));
					 //cout<<" "<<tempWallEnergy<<endl;
					 tempZ = fabs(wallz - tempZ);
					 //cout<<"Z-: "<<tempZ;
					 tempWallEnergy += factor*(0.4*pow(sigFS/tempZ, 10.0) - pow(sigFS/tempZ, 4.0) - pow(sigFS, 4.0)/3/thickz/pow(0.61*thickz+tempZ, 3.0));
					 //cout<<" "<<temp<<endl;
					 energyWall += tempWallEnergy;
				 }
			 }

			// cout<<"wall energy: "<<energyWall<<endl;
			 return energyWall;
		 }
			 
		 //Accessor
		 int getSiteNumber() {return siteNumber;}

		 string getName() {return name;}

		 double getMass() {return mass;}

		 bool getFixing() {return fixing;}

		 vector<Atom> getContainer() {return atomContainer;}

		 vector<double> getChargeContainer() {return chargeContainer;}

		 vector<double> getRX() {return rx;}

		 vector<double> getRY() {return ry;}

		 vector<double> getRZ() {return rz;}

		 Atom getAtom(int i) {return atomContainer[i];}

		 double centralDistance(const Molecule &m, double rangeX, double rangeY, double rangeZ)
		 {
			 Atom temp1, temp2;
			 temp1 = atomContainer[0];
			 temp2 = m.atomContainer[0];
			 double dist, distx, disty, distz;
			 distx = temp1.getX()-temp2.getX();
			 distx = distx-rangeX*(int(distx*2/rangeX)-int(distx/rangeX));
			 distx = distx*distx;
			 disty = temp1.getY()-temp2.getY();
			 disty = disty-rangeY*(int(disty*2/rangeY)-int(disty/rangeY));
			 disty = disty*disty;
			 distz = temp1.getZ()-temp2.getZ();
			 distz = distz-rangeZ*(int(distz*2/rangeZ)-int(distz/rangeZ));
			 distz = distz*distz;
			 dist = sqrt(distx+disty+distz);
			 return dist;
		 }

		 void moleculePdb(ofstream &file, int num)
		 {
			 for(int i = 0; i < atomContainer.size(); i++)
			 {
				 atomContainer[i].atomPdb(file, num*siteNumber+i);
			 }
		 }
};
#endif
