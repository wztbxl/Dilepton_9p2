#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "TString.h"

using namespace std;

void sortRunNumber(int nFiles) 
{

	system("mkdir -p RunList");

  char ch[10000];

  ofstream outfile("./minibias_runNumber.dat");

  for (int i = 0; i < nFiles; i++)
  {
    cout << "now at file " << i << endl;
    TString name = Form("/star/u/wangzhen/run20/Dielectron/BadRunList/BadRunList.dat");
    ifstream infile(name.Data());

    while(!infile.eof())
    {
      infile.getline(ch,9999);
      outfile << ch << endl;
    }
    infile.close();
  }
  outfile.close();


    vector<Int_t> runNumber; 
	vector<Int_t> dayNumber;
	runNumber.clear();
	dayNumber.clear();

	ifstream indata("minibias_runNumber.dat");

	Int_t runId;
	while(indata>>runId){ 
    cout << runId << endl;
		vector<Int_t>::iterator iter = find(runNumber.begin(), runNumber.end(), runId);
		if( iter == runNumber.end() ) runNumber.push_back(runId);
	} 
	cout<<"Total size of runNumber list: " << runNumber.size()<<endl;
	sort(runNumber.begin(),runNumber.end());

	Int_t dayId;
	ofstream outdata("RunList/mTotalRunList.dat");
	for(int i=0; i<runNumber.size(); i++){
		dayId = (runNumber[i]/1000)%1000;
		cout<<"RunId:"<<runNumber[i]<<" \t"<<"dayId:"<<dayId<<endl;
		outdata<<runNumber[i]<<endl;

		vector<Int_t>::iterator iter = find(dayNumber.begin(), dayNumber.end(), dayId);
		if( iter == dayNumber.end() ) dayNumber.push_back(dayId);
	}
	outdata.close();
	cout<<"Total size of dayNumber list: " << dayNumber.size()<<endl;
	sort(dayNumber.begin(),dayNumber.end());

	outdata.open("RunList/mTotalDayList.dat");
	for(int i=0; i<dayNumber.size(); i++){
		cout<<"dayId:"<<dayNumber[i]<<endl;
		outdata<<dayNumber[i]<<endl;
	}

	indata.close();
	outdata.close();
	cout<<"end of program"<<endl; 
}
