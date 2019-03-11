#ifndef ATOM_H
#define ATOM_H
#include <iostream>
#include <string>
#include <iomanip>
using namespace std;

class Atom
{
	private: double x, y, z;
		 string name;
		 double charge, sig, eps;
	public:
		 Atom()
		 {
			 x = 0.0; y = 0.0; z = 0.0;
			 name = "none";
			 charge = 0.0;
			 sig = 0;
			 eps = 0;
		 }

		 Atom(string n, double c, double s, double e)
		 {
			 name = n;
			 charge = c;
			 sig = s;
			 eps = e;
		 }

		 Atom(string n, double cx, double cy, double cz, double c, double s, double e)
		 {
			 x = cx; y = cy; z = cz;
			 name = n;
			 charge = c;
			 sig = s;
			 eps = e;
		 }

		 //Mutator
		 void save(ofstream &file)
		 {
			 file<<setw(5)<<name<<" "<<setprecision(10)<<x<<" "<<y<<" "<<z<<" "<<charge<<" "<<sig<<" "<<eps<<endl;
		 }

		 void read(ifstream &file)
		 {
			 file>>name>>x>>y>>z>>charge>>sig>>eps;
		 }

		 void setX(double cx)
		 {
			 x = cx;
		 }

		 void setY(double cy)
		 {
			 y = cy;
		 }

		 void setZ(double cz)
		 {
			 z = cz;
		 }

		 void setCoordination(double cx, double cy, double cz)
		 {
			 x = cx;
			 y = cy;
			 z = cz;
		 }

		 void setName(string n)
		 {
			 name = n;
		 }

		 void setCharge(double c)
		 {
			 charge = c;
		 }

		 void setSig(double s)
		 {
			 sig = s;
		 }

		 void setEps(double e)
		 {
			 eps = e;
		 }

		 void shift(double dx, double dy, double dz)
		 {
			 x += dx;
			 y += dy;
			 z += dz;
		 }

		 Atom operator-(Atom &right)
		 {
			 Atom temp(name,charge,sig,eps);
			 temp.x = x - right.x;
			 temp.y = y - right.y;
			 temp.z = z - right.z;

			 return temp;
		 }

		 Atom operator+(Atom &right)
		 {
			 Atom temp(name,charge,sig,eps);
			 temp.x = x + right.x;
			 temp.y = y + right.y;
			 temp.z = z + right.z;

			 return temp;
		 }

		 const Atom operator=(const Atom &right)
		 {
			 name = right.name;
			 x = right.x;
			 y = right.y;
			 z = right.z;
			 charge = right.charge;
			 sig = right.sig;
			 eps = right.eps;
			 return *this;
		 }


		 Atom rotate(double **rotate, int row)
		 {
			 Atom temp(name,charge,sig,eps);

			 temp.x = x*rotate[0][0] + y*rotate[1][0] + z*rotate[2][0];
			 temp.y = x*rotate[0][1] + y*rotate[1][1] + z*rotate[2][1];
			 temp.z = x*rotate[0][2] + y*rotate[1][2] + z*rotate[2][2];

			 return temp;
		 }

		 void align(double bx, int isX, double by, int isY, double bz, int isZ)
		 {
			 if(isX == 1)
				 x -= bx;
			 else if(isX == -1)
				 x += bx;
			 if(isY == 1)
				 y -= by;
			 else if(isY == -1)
				 y += by;
			 if(isZ == 1)
				 z -= bz;
			 else if(isZ == -1)
				 z += bz;
		 }

		 void showAtom()
		 {
			 cout<<name<<" "<<x<<" "<<y<<" "<<z<<" "<<charge<<" "<<sig<<" "<<eps<<endl;
		 }

		 void atomPdb(ofstream &file, int num)
		 {
			 file<<"HETATH"<<setw(5)<<num+1<<" "<<setw(3)<<left<<name<<"               "<<setw(8)<<setprecision(3)<<fixed<<right<<x<<setw(8)<<y<<setw(8)<<z<<endl;
		 }

		 //Accessor
		 double getX() {return x;}
		 
		 double getY() {return y;}

		 double getZ() {return z;}

		 string getName() {return name;}

		 double getCharge() {return charge;}

		 double getSig() {return sig;}

		 double getEps() {return eps;}

		 ~Atom()
		 {
		 }
};
#endif
