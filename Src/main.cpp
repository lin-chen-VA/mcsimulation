#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include "Molecule.h"
#include "System.h"
#include "functions.h"
using namespace std;

int main()
{
	System s;

	ifstream inputFile("system.in");
	int cycle;
	string str;
	inputFile>>str>>cycle;
	if(str != "cycle")
		cout<<"Parameter idratio error!"<<endl;

	int insertNumber;
	inputFile>>str>>insertNumber;
	if(str != "insertNumber")
		cout<<"Parameter insertNumber error!"<<endl;

	cout<<"Simulation will run "<<cycle<<" cycles, be patient..."<<endl;
	cout<<endl;
	cout<<"--------------------------------------------"<<endl;

	int timeAnalysis;
	inputFile>>str>>timeAnalysis;
	if(str != "timeAnalysis")
		cout<<"Parameter timeAnalysis error!"<<endl;

	int moleculeNumberPlot, moleculeNumberPlotStep;
	inputFile>>str>>moleculeNumberPlot;
	if(str != "moleculeNumberPlot")
		cout<<"Parameter moleculeNumberPlot error!"<<endl;
	inputFile>>str>>moleculeNumberPlotStep;
	if(str != "moleculeNumberPlotStep")
		cout<<"Parameter moleculeNumberPlotStep error!"<<endl;

	int checkFullEnergy, checkFullEnergyStep;
	inputFile>>str>>checkFullEnergy;
	if(str != "checkFullEnergy")
		cout<<"Parameter checkFullEnergy error!"<<endl;
	inputFile>>str>>checkFullEnergyStep;
	if(str != "checkFullEnergyStep")
		cout<<"Parameter checkFullEnergyStep error!"<<endl;

	int mediaOutput, mediaOutputStep;
	inputFile>>str>>mediaOutput;
	if(str != "mediaOutput")
		cout<<"Parameter mediaOutput error!"<<endl;
	inputFile>>str>>mediaOutputStep;
	if(str != "mediaOutputStep")
		cout<<"Parameter mediaOutputStep error!"<<endl;

	int grPlot, grStep;
	inputFile>>str>>grPlot;
	if(str != "grPlot")
		cout<<"Parameter grPlot error!"<<endl;
	inputFile>>str>>grStep;
	if(str != "grStep")
		cout<<"Parameter grStep error!"<<endl;

	double averageEnergy;
	int averageEnergyCount;
	if(checkFullEnergy == 1)
	{
		averageEnergy = 0.0;
		averageEnergyCount = 0;
	}

	ofstream numberPlot;
	if(moleculeNumberPlot == 1)
		numberPlot.open("number.out");

	time_t start;
	if(timeAnalysis == 1)
		time(&start);

	int idratio;
	idratio = s.getIdratio();

	int nfi = 0;
	double randomNumber;

	while(nfi <= cycle*(idratio+2))
	{
		if(s.getNotOnlyMove())
		{
		if(s.getMoleculeNumber() == 0) 
		{
			cout<<"Initial molecule number is zero, insert "<<insertNumber<<" molecules firstly to\n obtain a start configuration..."<<endl;
			while(s.getMoleculeNumber() < insertNumber)
				s.insertion();
		}
		randomNumber = s.ran();
		if(randomNumber < double(idratio)/double(idratio+2))
		{
			//cout<<"Insert trial."<<endl;
			s.insertion();
		}
		else if(randomNumber < double(idratio+1)/double(idratio+2))
		{
			//cout<<"Delete trial."<<endl;
			s.deletion();
		}
		else
		{
			//cout<<"Movement trial."<<endl;
			s.movement();
		}
		nfi++;
		}
		else
		{
			//cout<<"Only movement trial."<<endl;
			s.movement();
			nfi += 3;
		}

		if(moleculeNumberPlot == 1)
		{
			int tempcycle;
			tempcycle = nfi/(idratio+2);
			if(nfi%(moleculeNumberPlotStep*(idratio+2)) == 0)
				s.moleculeNumberPlot(tempcycle, numberPlot);
		}

		if(checkFullEnergy == 1)
		{
			if(nfi%(checkFullEnergyStep*(idratio+2)) == 0)
			{
				averageEnergy += s.averageFullEnergy();
				averageEnergyCount++;
			}
		}

		if(mediaOutput == 1)
		{
			if(nfi%(mediaOutputStep*(idratio+2)) == 0)
			{
				ofstream outputFile;
				outputFile.open("output.out");
				s.save(outputFile);
				outputFile.close();
			}
		}

		if(grPlot == 1)
		{
			if(nfi%(grStep*(idratio+2)) == 0)
			{
				ofstream grOutput;
				grOutput.open("gr.out");
				s.gr(grOutput);
				grOutput.close();
			}
		}
	}

	if(moleculeNumberPlot == 1)
		numberPlot.close();
	//cout<<"Molecule number: "<<s.getMoleculeNumber()<<endl;

	time_t end;
	if(timeAnalysis == 1)
	{
		s.timeAnalysis();
		time(&end);
		cout<<"Program run "<<end - start<<" seconds or "<<double(end - start)/3600<<" hours."<<endl;
	}

	s.generatePdb();

	cout<<"Average energy of system is: "<<averageEnergy/averageEnergyCount<<endl;

	cout<<"Successfully run."<<endl;

	return 0;
}
