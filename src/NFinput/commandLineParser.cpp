/*
 * commandLineParser.cpp
 *
 *  Created on: Oct 21, 2008
 *      Author: msneddon
 */

#include "NFinput.hh"





using namespace NFinput;
using namespace std;


bool NFinput::parseArguments(int argc, const char *argv[], map<string,string> &argMap)
{
	for(int a=1; a<argc; a++)
	{
		string s(argv[a]);

		//Find the strings that start with a minus, these are the flags
		if(s.compare(0,1,"-")==0)
		{
			string sFlag = s.substr(1,s.size()-1);
			if(sFlag.empty()) {
				cout<<"   Error in the command line arguments.  You gave a '-' with a space"<<endl;
				cout<<"   directly following. This is not a valid argument."<<endl<<endl;
				return false;
			}

			if(sFlag.compare(0,1,"-")==0){
				sFlag = s.substr(1,sFlag.size()-1);
			}
			if(sFlag.empty()) {
				cout<<"   Error in the command line arguments.  You gave a '--' with a space"<<endl;
				cout<<"   directly following. This is not a valid argument."<<endl<<endl;
				return false;
			}


			//See if the flag has some other input value that follows
			string sVal;
			if((a+1)<argc) {
				sVal=argv[a+1];
			}
			if(sVal.compare(0,1,"-")==0) {
				sVal = "";
			} else {
				a++;
			}

			// cout<<"found:  '"<<sFlag<<"' with arg: '"<<sVal<<"' "<<endl;

			if(argMap.find(sFlag)!=argMap.end()) {
				cout<<"Found two values for the same command line flag: '"<<sFlag<<"' so I am stopping"<<endl;
				return false;
			}

			argMap[sFlag] = sVal;
		}
		else
		{
			cout<<"   Warning when parsing command line arguments.  Valid arguments are preceded by a standard dash, as in '-logo'."<<endl;
			cout<<"\n   This argument: '"<< s <<"'"<<endl;
			cout<<"   did not begin with a proper dash and was ignored."<<endl<<endl;
		}
	}
	return true;
}



int NFinput::parseAsInt(map<string,string> &argMap,string argName,int defaultValue)
{
	if(argMap.find(argName)!=argMap.end()) {
		string strVal = argMap.find(argName)->second;
		try {
			int intVal = NFutil::convertToInt(strVal);
			return intVal;
		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+strVal+"' as an integer."<<endl;
			cout<<"!!  Using the default value of "<<defaultValue<<endl<<endl;;
		}
	}
	return defaultValue;
}


void NFinput::parseAsCommaSeparatedSequence(map<string,string> &argMap,string argName,vector<int> &sequence)
{
	if(argMap.find(argName)!=argMap.end()) {
		string argString = argMap.find(argName)->second;
		try {

			vector <string> numberStrings;
			numberStrings.push_back("");
			for(unsigned int i=0; i<argString.length(); i++)
			{
				if(argString.at(i)=='\"' || argString.at(i)==' ') { continue; }
				if(argString.at(i)==',') { numberStrings.push_back(""); continue; }
				numberStrings.at(numberStrings.size()-1) = numberStrings.at(numberStrings.size()-1) + argString.at(i);
			}

			for(unsigned int i=0; i<numberStrings.size(); i++) {
				sequence.push_back(NFutil::convertToInt(numberStrings.at(i)));
			}

		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+argString+"' as a comma separated "+
					"integer sequence, so I'm quitting."<<endl;
			throw std::runtime_error("Could not parse as a comma separated integer sequence");
		}
	}
}

double NFinput::parseAsDouble(map<string,string> &argMap,string argName,double defaultValue)
{
	if(argMap.find(argName)!=argMap.end()) {
		string strVal = argMap.find(argName)->second;
		try {
			double doubleVal = NFutil::convertToDouble(strVal);
			return doubleVal;
		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+strVal+"' as a double."<<endl;
			cout<<"!!  Using the default value of "<<defaultValue<<endl<<endl;;
		}
	}
	return defaultValue;
}



bool NFinput::parseSequence(const string& numString, vector <double> &outputTimes)
{
	double startVal=0, stepVal=1, endVal=0;
	try {

		string::size_type c1 = numString.find_first_of(':');
		if(c1!=string::npos) {
			string::size_type c2 = numString.find_first_of(':',c1+1);
			if(c2!=string::npos) {
					startVal= NFutil::convertToDouble(numString.substr(0,c1));
					stepVal= NFutil::convertToDouble(numString.substr(c1+1,c2-c1-1));
					endVal= NFutil::convertToDouble(numString.substr(c2+1));

			} else {
					startVal= NFutil::convertToDouble(numString.substr(0,c1));
					endVal= NFutil::convertToDouble(numString.substr(c1+1));
			}
		}

	} catch(std::runtime_error e) {
		return false;
	}

	if(startVal>endVal) {
		cout<<"Error: start value of sequence must be <= end value."<<endl;
		return false;
	} else if(stepVal<=0) {
		cout<<"Error: step value of sequence must be >0."<<endl;
		return false;
	} else {
		if(outputTimes.size()>=1)
			if(startVal<=outputTimes.at(outputTimes.size()-1)) {
				cout<<"\n\nError in NFinput::creatComplexOutputDumper: output times given \n";
				cout<<"must be monotonically increasing without any repeated elements.";
				return false;
			}
		//Only if everything went as planned to we then add the output steps accordingly
		for(double d=startVal; d<=endVal; d+=stepVal) {
			outputTimes.push_back(d);
		}
		return true;
	}
	return true;
}



bool NFinput::createSystemDumper(const string& paramStr, System *s, bool verbose)
