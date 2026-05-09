#include "NFinput.hh"
#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <sstream>





using namespace NFinput;
using namespace std;

namespace {

string tfun_trim_copy(string s) {
	while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
		s.erase(s.begin());
	}
	while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
		s.pop_back();
	}
	return s;
}

string tfun_to_lower_copy(string s) {
	std::transform(s.begin(), s.end(), s.begin(),
		[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return s;
}

bool tfun_is_time_name(const string &name) {
	string normalized = tfun_to_lower_copy(tfun_trim_copy(name));
	return normalized == "time" || normalized == "t" ||
		normalized == "time()" || normalized == "t()";
}

bool tfun_parse_csv_numbers(const string &csv, vector<double> &values, string &error) {
	values.clear();
	std::stringstream ss(csv);
	string token;
	while (std::getline(ss, token, ',')) {
		token = tfun_trim_copy(token);
		if (token.empty()) {
			error = "encountered empty numeric token";
			return false;
		}

		const char *start = token.c_str();
		char *end = NULL;
		errno = 0;
		double v = std::strtod(start, &end);
		if (start == end || errno == ERANGE) {
			error = "invalid numeric token '" + token + "'";
			return false;
		}
		while (*end != '\0' && std::isspace(static_cast<unsigned char>(*end))) {
			++end;
		}
		if (*end != '\0') {
			error = "unexpected trailing characters in token '" + token + "'";
			return false;
		}
		values.push_back(v);
	}

	if (values.empty()) {
		error = "no numeric values were provided";
		return false;
	}
	return true;
}

}  // namespace



bool createLocalFunction(string name,
			string expression,
			vector <string> &argNames,
			vector <string> &refNames,
			vector <string> &refTypes,
			System *s,
			map <string,double> &parameter,
			TiXmlElement * pListOfObservables,
			map<string,int> &allowedStates,
			bool verbose);

bool createCompositeFunction(string name,
			string expression,
			vector <string> &argNames,
			vector <string> &refNames,
			vector <string> &refTypes,
			vector <string> &paramNames,
			System *s,
			bool verbose);




bool createFunction(string name,
		string expression,
		vector <string> &argNames,
		vector <string> &refNames,
		vector <string> &refTypes,
		System *s,
		map <string,double> &parameter,
		TiXmlElement * pListOfObservables,
		map<string,int> &allowedStates,
		bool verbose)
{

	if(expression.size()==0) return true;

	vector <string> varRefNames;
	vector <string> varRefTypes;
	vector <string> paramNames;
	int otherFuncRefCounter=0;


	for(unsigned int rn=0; rn<refNames.size(); rn++) {
		if(refTypes.at(rn)=="Function") {
			otherFuncRefCounter++;
		} else if(refTypes.at(rn)=="Constant") {
			paramNames.push_back(refNames.at(rn));
		} else if(refTypes.at(rn)=="ConstantExpression") {
			paramNames.push_back(refNames.at(rn));
		} else {
			varRefNames.push_back(refNames.at(rn));
			varRefTypes.push_back(refTypes.at(rn));
		}
	}


	if(otherFuncRefCounter!=0) {

		//Must be a composite function, so parse as such
		//cout<<"creating composite function\n";
		createCompositeFunction(name, expression, argNames, refNames, refTypes, paramNames, s, verbose);
		return true;
		//	exit(1);
	}
	else if(argNames.size()==0) //&& otherFuncRefCounter==0)
	{
		//must be a global function, as we have no arguments, so just create it
		//cout<<"creating global function\n";
		GlobalFunction *gf = new GlobalFunction(name, expression,
												varRefNames, varRefTypes, paramNames, s);
		if(!s->addGlobalFunction(gf)) {
			cerr<<"!!!Error:  Function name '"<<name<<"' has already been used.  You can't have two\n";
			cerr<<"functions with the same name, so I'll just stop now."<<endl;
			return false;
		}
		return true;
	}


//	if(expression.size()==0) {
//		cerr<<"!!!Error:  Local function named '"<<name<<"' has no Expression!\n";
//		return false;
//	}
	// else if(argNames.size()>0 && otherFuncRefCounter==0)
	//if we got here, we are creating a local function, so call the create local function function.
	//cout<<"creating local function"<<endl;
	return createLocalFunction(name, expression,
			argNames, refNames, refTypes,
			s, parameter, pListOfObservables, allowedStates, verbose);
}




bool createCompositeFunction(string name,
			string expression,
			vector <string> &argNames,
			vector <string> &refNames,
			vector <string> &refTypes,
			vector <string> &paramNames,
			System *s,
			bool verbose)
{
//	cout<<"must be a composite function..."<<endl;

	for(unsigned int rn=0; rn<refNames.size(); rn++) {
		if(refTypes.at(rn)=="Observable" || refTypes.at(rn)=="MoleculeObservable" || refTypes.at(rn)=="SpeciesObservable") {
			cerr<<"Composite Functions (functions that call other functions) cannot have"<<endl;
			cerr<<"references to observables.  You must put those in base level functions."<<endl;
			exit(1);
		}
	}

	vector <string> functionsCalled;
	for(unsigned int rn=0; rn<refNames.size(); rn++) {
		if(refTypes.at(rn)=="Function") {
			functionsCalled.push_back(refNames.at(rn));
		}
	}


	CompositeFunction *cf = new CompositeFunction(s,name, expression,functionsCalled,argNames,paramNames);
	s->addCompositeFunction(cf);

	return true;
}






bool createLocalFunction(string name,
			string expression,
			vector <string> &argNames,
			vector <string> &refNames,
			vector <string> &refTypes,
			System *s,
			map <string,double> &parameter,
			TiXmlElement * pListOfObservables,
			map<string,int> &allowedStates,
			bool verbose)
{

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

//	cout<<"must be a local function..."<<endl;

	//Remember the original expression (useful for outputting)
	string originalExpression = expression;

	// First, extract out the parameter constants we will use, and make sure
	// that this isn't just a function reference to a local function...
	vector <string> paramNames;
	int otherFuncRefCounter=0;
	for(unsigned int rn=0; rn<refNames.size(); rn++) {
		if(refTypes.at(rn)=="Function") {
			otherFuncRefCounter++;
		} else if(refTypes.at(rn)=="Constant") {
			paramNames.push_back(refNames.at(rn));
		} else if(refTypes.at(rn)=="ConstantExpression") {
			paramNames.push_back(refNames.at(rn));
		}
	}

	if(otherFuncRefCounter==1) {
		//handle the function reference...
		cout<<"local function has an reference to another function!"<<endl;
		cout<<"you should have created this with a composite function constructor!"<<endl;
		return false;
	}
	if(otherFuncRefCounter!=0) {
		cerr<<"!!!Error:  Functions can reference at most one other function!  Quitting."<<endl;
		return false;
	}


	//some vectors to store results
	vector <string> obsUsedName;
	vector <string> obsUsedExpressionRef;
	vector <int> obsUsedScope;
	vector <bool> markForRemoval;

	//First, figure out the scope of the observables that we
	//referenced in the function.  This requires the following annoying
	//set of nested loops:
	for(unsigned int rn=0; rn<refNames.size(); rn++) {
		if(refTypes.at(rn)=="Observable" || refTypes.at(rn)=="MoleculeObservable" || refTypes.at(rn)=="SpeciesObservable") {
			string::size_type sPos=expression.find(refNames.at(rn));
			for( ; sPos!=string::npos; sPos=expression.find(refNames.at(rn),sPos+1)) {

				//first, add this reference to our vector lists
				obsUsedName.push_back(refNames.at(rn));
				obsUsedExpressionRef.push_back(refNames.at(rn));
				obsUsedScope.push_back(-1); //assume global scope (with -1) until we discover otherwise
				markForRemoval.push_back(false); //assume that we aren't removing this guy, yet

				//Now, go about finding the scope by finding the next enclosing parentheses...
				string::size_type openPar = expression.find_first_of('(',sPos);
				string::size_type closePar = expression.find_first_of(')',sPos);
				if(openPar!=string::npos && closePar!=string::npos) {
					if(closePar>openPar) { //if we got here, we found a valid parenthesis to look at

						string possibleArg = expression.substr(openPar+1,closePar-openPar-1);
						NFutil::trim(possibleArg);

						for(unsigned int aIndex=0; aIndex<argNames.size(); aIndex++) {
							if(argNames.at(aIndex)==possibleArg) {
								//hurah!  we found the local scope of this guy

								//Now we have to gently excise the reference, and replace it with
								//a marker so we can uniquely identify it later...
								int lastIndex = obsUsedExpressionRef.size();
								string identifier = "_"+NFutil::toString(lastIndex);
								expression.replace(openPar,closePar-openPar+1,identifier);
								obsUsedExpressionRef.pop_back();
								obsUsedExpressionRef.push_back(obsUsedName.at(lastIndex-1)+identifier);

								obsUsedScope.pop_back();
								obsUsedScope.push_back(aIndex);
								break; //break cause we're done with this scope...
							}
						}
					}
				} //end if statement to find open and closed parentheses
			}//end for loop over sPos
		}
	} //end loop over the possible references


	//as an optimization, certain observables may be used in more than one position
    //with the same scope.  The loops above uniquely identified each instance of each
	//reference - but that means that the function will evaluate separately each
	//instance found in the expression.  If we have the same reference with the same
	//scope, then we can reduce by renaming them as a single reference.

	//so we loop over each used reference...
	for(unsigned int i=0; i<obsUsedExpressionRef.size(); i++) {

		//compare this reference to each of the other used references...
		for(unsigned int j=i+1; j<obsUsedExpressionRef.size(); j++) {

			//If we are already removing element j, then we already found it
			//and replaced it, so we can stop.
			if(markForRemoval.at(j)) continue;

			//cout<<"comparing: "<<obsUsedExpressionRef.at(i)<<" to "<< obsUsedExpressionRef.at(j)<<endl;
			if(obsUsedName.at(i)==obsUsedName.at(j)) {
				if(obsUsedScope.at(i)==obsUsedScope.at(j)) {
					//find it in the expression, replace it with the original, and mark it for removal
					string::size_type refPos =expression.find(obsUsedExpressionRef.at(j));
					expression.replace(refPos, obsUsedExpressionRef.at(j).size(),obsUsedExpressionRef.at(i));
					markForRemoval.at(j)=true;
				}
			}
		}
	}

	////////////////////
	// almost there...

	//Final vectors for storing the list after the removal (it is easier to just copy
	//than risk errors on indexing when we are removing from the original...)
	vector <string> finalObsUsedExpressionRef;
	vector <string> finalObsUsedName;
	vector <int> finalObsUsedScope;
	vector <Observable *> finalLocalObservables;
	//Fill the final vectors
	for(unsigned int i=0; i<obsUsedExpressionRef.size(); i++) {
		if(!markForRemoval.at(i)) {
			finalObsUsedExpressionRef.push_back(obsUsedExpressionRef.at(i));
			finalObsUsedName.push_back(obsUsedName.at(i));
			finalObsUsedScope.push_back(obsUsedScope.at(i));
			finalLocalObservables.push_back(0);
		}
	}

	//now from the list of observable names, we have to go back and
	//actually create new observables for these things, by cloning existing
	//observables.  So we'll do that here.
	for(unsigned int i=0; i<finalObsUsedName.size(); i++) {
		if(finalObsUsedScope.at(i)!=-1) {
			Observable *systemObs=s->getObservableByName(finalObsUsedName.at(i));
			if(systemObs==0) {
				cout.flush();
				cerr<<"Error!! unable to create local function '"<<name<<"' because the reference to\n";
				cerr<<"observable named: '"<<finalObsUsedName.at(i)<<"' was not found."<<endl;
				return false;
			}

			finalLocalObservables.at(i) = systemObs->clone();
		}
	}

//  DEPRECATED - original version for matching up observables, now handled by cloning existing ones
//	//now from the list of observable names, we have to go back and
//	//actually create observables for these things.  So we'll do that here
//	try {
//		TiXmlElement *pObs;
//		for ( pObs = pListOfObservables->FirstChildElement("Observable");
//			pObs != 0; pObs = pObs->NextSiblingElement("Observable")) {
//
//			string observableId="", observableName="", observableType="";
//			if(!pObs->Attribute("id") || !pObs->Attribute("name") || !pObs->Attribute("type")) {
//				cerr<<"Observable tag without a valid 'id' attribute.  Quiting"<<endl;
//				return false;
//			} else {
//				observableId = pObs->Attribute("id");
//				observableName = pObs->Attribute("name");
//				observableType = pObs->Attribute("type");
//			}
//
//			for(unsigned int i=0; i<finalObsUsedExpressionRef.size(); i++) {
//				if(observableName==finalObsUsedName.at(i)) {
//					if(finalObsUsedScope.at(i)!=-1) {
//
//						Observable *systemObs=s->getObservableByName(observableName);
//						if(systemObs==0) {
//							cerr<<"Observable in local function "<<name<<" was not found."<<endl;
//							return false;
//						}
//
//						finalLocalObservables.at(i) = systemObs->clone();
////						if(observableType.compare("Molecules")==0) {
////							cout<<"In local function!!!  Creating an observable!!"<<endl;
////							vector <TemplateMolecule *> tmList;
////							readObservableForTemplateMolecules(pObs,observableName,tmList,
////									s,parameter,allowedStates,false);
////						} else if()
////
//
////						TemplateMolecule *tempmol = 0;
////						if(!readObservable(pObs,observableName,tempmol,
////								s,parameter, allowedStates,false)) {
////							return false;
////						}
////						if(tempmol!=0) {
////							Observable *o  = new Observable(observableName.c_str(),tempmol);
////							finalLocalObservables.at(i)=o;
////							// WE DO NOT ADD THIS OBSERVABLE TO THE MOLECULE HERE!!
////						}
//					}
//					break;
//				}
//			}
//		}
//
//	} catch (...) {
//		//oh no, what happened this time?
//		cout<<"Caught an unexpected error when creating observables for this local function."<<endl;
//		return false;
//	}
//
//
//	//make sure everything went as planned
//	for(unsigned int i=0; i<finalObsUsedExpressionRef.size(); i++) {
//		if(finalObsUsedScope.at(i)!=-1) {
//			if(finalLocalObservables.at(i)==0) {
//				cout<<"Error!! unable to create local function because the reference to\n";
//				cout<<"observable named: '"<<finalObsUsedName.at(i)<<"' was not found."<<endl;
//				return false;
//			}
//		}
//	}
//  END DEPRECATED CODE






	//Ah.  and so we get here.  We now have:
	// 1) original functional expression (originalExpression)
	// 2) new revised functional expression with scope notation removed (expression)
	// 3) a list of each observable referenced and its referenced name (finalObsUsedExpressionRef)
	// 4) a list of what observables those references depend on (finalObsUsedName)
	// 5) a list of the scope of each of those references (finalObsUsedScope)
	// 6) the original list of arguments (argNames)
	// 7) the reduced list of parameter constants (paramNames)

	//so we can finally create our local function...

	LocalFunction * lf = new LocalFunction(s,name,
						originalExpression, expression,
						argNames,
						finalObsUsedExpressionRef,finalObsUsedName,finalLocalObservables,finalObsUsedScope,
						paramNames);
	s->addLocalFunction(lf);
	//cout<<"was added fine."<<endl;
	return true;
}















////  New Function Parser
bool NFinput::initFunctions(
	TiXmlElement * pListOfFunctions,
	System * system,
	map <string,double> &parameter,
	TiXmlElement * pListOfObservables,
	map<string,int> &allowedStates,
	bool verbose)
{
	try {
		vector <string> argNames;
		vector <string> refNames;
		vector <string> refTypes;
		vector <string> refNamesSorted;
		vector <string> refTypesSorted;

		//Loop through the Function tags...
		TiXmlElement *pFunction;
		for ( pFunction = pListOfFunctions->FirstChildElement("Function"); pFunction != 0; pFunction = pFunction->NextSiblingElement("Function"))
		{
			//Check if MoleculeType tag has a name...
			if(!pFunction->Attribute("id")) {
				cerr<<"!!!Error:  Function tag must contain the id attribute.  Quitting."<<endl;
				return false;
			}

			//Read in the Function Name
			string funcName = pFunction->Attribute("id");
			if(verbose) cout<<"\t\tReading Function: "+funcName+"(";


			//Get the list of arguments for this function
			TiXmlElement *pListOfArgs = pFunction->FirstChildElement("ListOfArguments");
			if(pListOfArgs) {
				//Loop through each argument
				TiXmlElement *pArg; bool firstArg = true;
				for ( pArg = pListOfArgs->FirstChildElement("Argument"); pArg != 0; pArg = pArg->NextSiblingElement("Argument"))
				{
					if(verbose && !firstArg) cout<<", ";
					firstArg=false;
					//Check again for errors by making sure the argument has a name
					string argName;
					if(!pArg->Attribute("id")) {
						if(verbose) cout<<" ?? ...\n";
						cerr<<"!!!Error:  Argument tag in Function: '" + funcName + "' must contain the id attribute.  Quitting."<<endl;
						return false;
					}

					argName = pArg->Attribute("id");
					argNames.push_back(argName);
					if(verbose) cout<<argName;
				}
			}
			if(verbose) cout<<")"<<endl;

			//Read in the actual function definition
			string funcExpression = "";
			TiXmlElement *pExpression = pFunction->FirstChildElement("Expression");
			if(pExpression) {
				if(!pExpression->GetText()) {
					if(funcName.substr(0,9).compare("reactant_")!=0) {
						cout<<"\n!!!Warning:  Expression tag for function "+funcName +" does not exist!  Function will not be generated."<<endl;
					}
					continue;
					//return false;
				}
				funcExpression = pExpression->GetText();
				if(verbose) cout<<"\t\t\t = "<<funcExpression<<endl;
			} else {
				cout<<"\n!!!Warning:  Expression tag for a function does not exist!  Function will not be generated."<<endl;
				continue;
				//return false;
			}

			//Get the list of References
			TiXmlElement *pListOfRefs = pFunction->FirstChildElement("ListOfReferences");
			if(pListOfRefs) {

				//Loop through each reference
				TiXmlElement *pRef;
				for ( pRef = pListOfRefs->FirstChildElement("Reference"); pRef != 0; pRef = pRef->NextSiblingElement("Reference"))
				{
					//Check again for errors by making sure the parameter has a name
					if(!pRef->Attribute("name")) {
						cerr<<"!!!Error:  Reference tag in Function: '" + funcName + "' must contain a proper id.  Quitting."<<endl;
						return false;
					} else if(!pRef->Attribute("type")) {
						cerr<<"!!!Error:  Reference tag in Function: '" + funcName + "' must contain a proper type.  Quitting."<<endl;
						return false;
					}

					string refName = pRef->Attribute("name");
					string refType = pRef->Attribute("type");

					if(verbose) cout<<"\t\t\t\tReference: "+refType+" "+refName<<endl;

					refNames.push_back(refName);
					refTypes.push_back(refType);
				}
			}


			// simple sort to order the elements, so that overlapping names do not
			// get parsed wrong when interpreting the functions.
			while(refNames.size()>0)
			{
				unsigned int maxLength = 0;
				int maxIndex = 0;
				string maxName = "";
				string maxType = "";

				for(unsigned int k=0; k<refNames.size(); k++)
				{
					if(refNames.at(k).length()>maxLength) {
						maxName = refNames.at(k);
						maxType = refTypes.at(k);
						maxIndex = k;
						maxLength = refNames.at(k).length();
					}
				}

				// pop off the max value
				refNames.at(maxIndex) = refNames.at(refNames.size()-1);
				refTypes.at(maxIndex) = refTypes.at(refTypes.size()-1);
				refNames.pop_back();
				refTypes.pop_back();

				// add it to the sorted list
				refNamesSorted.push_back(maxName);
				refTypesSorted.push_back(maxType);
			}



			//Here we actually generate the function or the function generator
			if(!createFunction(funcName,
					funcExpression,
					argNames,
					refNamesSorted,
					refTypesSorted,
					system,
					parameter,
					pListOfObservables,
					allowedStates,
					verbose)) {
				return false;
			}

			// Hook up time reference if any
			bool hasTimeRef = false;
			for(unsigned int i=0; i<refTypesSorted.size(); i++) {
				if(refTypesSorted.at(i) == "Time") {
					hasTimeRef = true;
					break;
				}
			}
			if (hasTimeRef) {
				GlobalFunction *gf = system->getGlobalFunctionByName(funcName);
				if (gf) {
					gf->addSystemPointer(system);
				} else {
					CompositeFunction *cf = system->getCompositeFunctionByName(funcName);
					if (cf) {
						cf->addSystemPointer(system);
					}
				}
				system->setHasTimeDependentFunctions(true);
			}

			// AS-2021
			// check to see if it has a type and if yes, if it's of type TFUN
			if(pFunction->Attribute("type")) {
				string funcType = pFunction->Attribute("type");
				if (funcType == "TFUN") {
					const string tfunPlaceholder = "__TFUN_VAL__";
					const string legacyTfunPlaceholder = "__TFUN__VAL__";
					const bool hasNewPlaceholder = (funcExpression.find(tfunPlaceholder) != string::npos);
					const bool hasLegacyPlaceholder = (funcExpression.find(legacyTfunPlaceholder) != string::npos);
					string activePlaceholder;
					string ctrName;
					string ctrType;
					string mode;
					string method;
					string filePath;
					string xDataCsv;
					string yDataCsv;
					vector<double> inlineXs;
					vector<double> inlineYs;
					string csvError;

					if (hasNewPlaceholder && hasLegacyPlaceholder) {
						cerr<<"!!!Error:  TFUN function "<<funcName
						    <<" cannot mix "<<tfunPlaceholder<<" and "<<legacyTfunPlaceholder
						    <<" in the same expression.  Quitting."<<endl;
						return false;
					}
					if (!hasNewPlaceholder && !hasLegacyPlaceholder) {
						cerr<<"!!!Error:  TFUN function "<<funcName
						    <<" expression must contain "<<tfunPlaceholder
						    <<" or "<<legacyTfunPlaceholder<<".  Quitting."<<endl;
						return false;
					}
					activePlaceholder = hasNewPlaceholder ? tfunPlaceholder : legacyTfunPlaceholder;

					if (!pFunction->Attribute("ctrName")) {
						cerr<<"!!!Error:  Can't find counter name for TFUN function "<<funcName<<".  Quitting."<<endl;
						return false;
					}
					ctrName = tfun_trim_copy(pFunction->Attribute("ctrName"));
					if (ctrName.empty()) {
						cerr<<"!!!Error:  TFUN function "<<funcName<<" has empty ctrName.  Quitting."<<endl;
						return false;
					}

					if (pFunction->Attribute("mode")) {
						mode = tfun_to_lower_copy(tfun_trim_copy(pFunction->Attribute("mode")));
					}
					if (pFunction->Attribute("method")) {
						method = tfun_to_lower_copy(tfun_trim_copy(pFunction->Attribute("method")));
					}
					if (method.empty()) {
						method = (activePlaceholder == legacyTfunPlaceholder) ? "step" : "linear";
					}
					if (method != "linear" && method != "step") {
						cerr<<"!!!Error:  TFUN function "<<funcName<<" has unsupported method '"<<method
						    <<"'. Expected linear or step.  Quitting."<<endl;
						return false;
					}

					bool hasFile = pFunction->Attribute("file");
					bool hasXData = pFunction->Attribute("xData");
					bool hasYData = pFunction->Attribute("yData");

					if (hasFile) filePath = tfun_trim_copy(pFunction->Attribute("file"));
					if (hasXData) xDataCsv = pFunction->Attribute("xData");
					if (hasYData) yDataCsv = pFunction->Attribute("yData");

					if (mode.empty() && hasFile && (hasXData || hasYData)) {
						cerr<<"!!!Error:  TFUN function "<<funcName
						    <<" specifies both file and inline data attributes without an explicit mode.  Quitting."<<endl;
						return false;
					}
					if (mode.empty()) {
						if (hasXData || hasYData) mode = "inline";
						else if (hasFile) mode = "file";
					}

					if (mode != "file" && mode != "inline") {
						cerr<<"!!!Error:  TFUN function "<<funcName<<" has invalid mode '"<<mode
						    <<"'. Expected file or inline.  Quitting."<<endl;
						return false;
					}
					if (mode == "file" && (!hasFile || filePath.empty())) {
						cerr<<"!!!Error:  TFUN function "<<funcName
						    <<" in file mode requires non-empty file attribute.  Quitting."<<endl;
						return false;
					}
					if (mode == "file" && (hasXData || hasYData)) {
						cerr<<"!!!Error:  TFUN function "<<funcName
						    <<" in file mode cannot include xData/yData attributes.  Quitting."<<endl;
						return false;
					}
					if (mode == "inline") {
						if (hasFile) {
							cerr<<"!!!Error:  TFUN function "<<funcName
							    <<" in inline mode cannot include file attribute.  Quitting."<<endl;
							return false;
						}
						if (!hasXData || !hasYData) {
							cerr<<"!!!Error:  TFUN function "<<funcName
							    <<" in inline mode requires both xData and yData.  Quitting."<<endl;
							return false;
						}
						if (!tfun_parse_csv_numbers(xDataCsv, inlineXs, csvError)) {
							cerr<<"!!!Error:  TFUN function "<<funcName
							    <<" invalid xData ("<<csvError<<").  Quitting."<<endl;
							return false;
						}
						if (!tfun_parse_csv_numbers(yDataCsv, inlineYs, csvError)) {
							cerr<<"!!!Error:  TFUN function "<<funcName
							    <<" invalid yData ("<<csvError<<").  Quitting."<<endl;
							return false;
						}
						if (inlineXs.size() != inlineYs.size()) {
							cerr<<"!!!Error:  TFUN function "<<funcName<<" has mismatched xData/yData lengths ("
							    <<inlineXs.size()<<" vs "<<inlineYs.size()<<").  Quitting."<<endl;
							return false;
						}
						if (inlineXs.size() < 2) {
							cerr<<"!!!Error:  TFUN function "<<funcName<<" requires at least 2 inline rows.  Quitting."<<endl;
							return false;
						}
						for (size_t i = 1; i < inlineXs.size(); ++i) {
							if (inlineXs[i] <= inlineXs[i - 1]) {
								cerr<<"!!!Error:  TFUN function "<<funcName
								    <<" xData must be strictly increasing.  Quitting."<<endl;
								return false;
							}
						}
					}

					if (tfun_is_time_name(ctrName)) {
						ctrType = "Time";
					} else {
						for (unsigned int r = 0; r < refNamesSorted.size(); ++r) {
							if (refNamesSorted.at(r) != ctrName) continue;
							string refType = refTypesSorted.at(r);
							if (refType == "Observable") ctrType = "Observable";
							else if (refType == "Function") ctrType = "Function";
							else if (refType == "Constant" || refType == "ConstantExpression") ctrType = "Parameter";
							break;
						}
					}

					if (ctrType.empty()) {
						if (parameter.find(ctrName) != parameter.end()) {
							ctrType = "Parameter";
						} else if (system->getGlobalFunctionByName(ctrName) != NULL) {
							ctrType = "Function";
						} else if (system->getObservableByName(ctrName) != NULL) {
							ctrType = "Observable";
						}
					}
					if (ctrType.empty()) {
						cerr<<"!!!Error:  TFUN function "<<funcName<<" has unsupported ctrName '"<<ctrName
						    <<"' (expected time/t, parameter, observable, or function).  Quitting."<<endl;
						return false;
					}

					GlobalFunction *gf = system->getGlobalFunctionByName(funcName);
					CompositeFunction *cf = system->getCompositeFunctionByName(funcName);
					if (!gf && !cf) {
						cerr<<"!!!Error:  Could not find created TFUN function object '"<<funcName<<"'.  Quitting."<<endl;
						return false;
					}

					if (gf) {
						if (mode == "file") gf->enableFileDependency(filePath, method);
						else gf->enableInlineDependency(inlineXs, inlineYs, method);
						gf->setCtrName(activePlaceholder);
					}
					if (cf) {
						if (mode == "file") cf->enableFileDependency(filePath, method);
						else cf->enableInlineDependency(inlineXs, inlineYs, method);
						cf->setCtrName(activePlaceholder);
					}

					if (ctrType == "Observable") {
						Observable *obs = system->getObservableByName(ctrName);
						if (!obs) {
							cerr<<"!!!Error:  TFUN function "<<funcName<<" observable counter '"<<ctrName
							    <<"' was not found.  Quitting."<<endl;
							return false;
						}
						if (gf) obs->addReferenceToGlobalFunction(gf);
						if (cf) obs->addReferenceToCompositeFunction(cf);
					} else if (ctrType == "Time") {
						if (gf) gf->setCounterFromTime(system);
						if (cf) cf->setCounterFromTime(system);
						system->setHasTimeDependentFunctions(true);
					} else if (ctrType == "Parameter") {
						if (gf) gf->setCounterFromParameter(system, ctrName);
						if (cf) cf->setCounterFromParameter(system, ctrName);
					} else if (ctrType == "Function") {
						GlobalFunction *ctrFunc = system->getGlobalFunctionByName(ctrName);
						if (!ctrFunc) {
							cerr<<"!!!Error:  TFUN function "<<funcName<<" counter function '"<<ctrName
							    <<"' was not found.  Quitting."<<endl;
							return false;
						}
						if (!cf) {
							cerr<<"!!!Error:  TFUN function "<<funcName
							    <<" with function counter requires composite-function form.  Quitting."<<endl;
							return false;
						}
						cf->addFunctionPointer(ctrFunc);
					}
				}
			}
			// AS-2021

			//And here we clear our arrays
			argNames.clear();
			refNames.clear();
			refTypes.clear();
			refNamesSorted.clear();
			refTypesSorted.clear();
		}



		//Once we've read in all the functions, we should take care of finalizing
		//the composite functions so they properly reference the other functions
		system->finalizeCompositeFunctions();
	//	cout<<"done reading functions!"<<endl;

		//Getting here means we read everything we could successfully
		return true;
	} catch (...) {
		//Uh oh! we got some unknown exception thrown, so we must abort!
		cerr<<"I caught some unknown error when I was trying to parse out a function.\n";
		cerr<<"I'm at a loss for words right now, so you're on you're own."<<endl;
		return false;
	}
}










//////  Original Function Parser
//bool NFinput::initGlobalFunctions(
//	TiXmlElement * pListOfFunctions,
//	System * system,
//	map <string,double> &parameter,
//	bool verbose)
//{
//	try {
//		vector <string> argNames;
//		vector <string> argTypes;
//		vector <string> paramNames;
//		vector <double> paramValues;
//
//		//Loop through the Function tags...
//		TiXmlElement *pFunction;
//		for ( pFunction = pListOfFunctions->FirstChildElement("Function"); pFunction != 0; pFunction = pFunction->NextSiblingElement("Function"))
//		{
//			//Check if MoleculeType tag has a name...
//			if(!pFunction->Attribute("id")) {
//				cerr<<"!!!Error:  Function tag must contain the id attribute.  Quitting."<<endl;
//				return false;
//			}
//
//			//Read in the Function Name
//			string funcName = pFunction->Attribute("id");
//			if(verbose) cout<<"\t\tReading and Creating Function: "+funcName+"(";
//
//
//			//Get the list of arguments for this function
//			TiXmlElement *pListOfArgs = pFunction->FirstChildElement("ListOfArgs");
//			if(pListOfArgs) {
//				//Loop through each arguement
//				TiXmlElement *pArg; bool firstArg = true;
//				for ( pArg = pListOfArgs->FirstChildElement("Arg"); pArg != 0; pArg = pArg->NextSiblingElement("Arg"))
//				{
//					if(verbose && !firstArg) cout<<", ";
//					firstArg=false;
//					//Check again for errors by making sure the argument has a name
//					string argName, argType;
//					if(!pArg->Attribute("id") || !pArg->Attribute("name") || !pArg->Attribute("type")) {
//						if(verbose) cout<<" ?? ...\n";
//						cerr<<"!!!Error:  Arg tag in Function: '" + funcName + "' must contain the all attributes including id, name, and type.  Quitting."<<endl;
//						return false;
//					}
//
//					argName = pArg->Attribute("name");
//					argType = pArg->Attribute("type");
//					argNames.push_back(argName);
//					argTypes.push_back(argType);
//
//					if(verbose) cout<<argName;
//				}
//			}
//			if(verbose) cout<<")"<<endl;
//
//
//			//Get the list of Parameter Constants for this function
//			TiXmlElement *pListOfParamConst = pFunction->FirstChildElement("ListOfParameterConstants");
//			if(pListOfParamConst) {
//
//				//Loop through each arguement
//				TiXmlElement *pParamConst;
//				for ( pParamConst = pListOfParamConst->FirstChildElement("ParameterConstant"); pParamConst != 0; pParamConst = pParamConst->NextSiblingElement("ParameterConstant"))
//				{
//					//Check again for errors by making sure the parameter has a name
//					if(!pParamConst->Attribute("id")) {
//						cerr<<"!!!Error:  ParameterConstant tag in Function: '" + funcName + "' must contain a proper id.  Quitting."<<endl;
//						return false;
//					}
//
//					string paramName = pParamConst->Attribute("id");
//					if(parameter.find(paramName)==parameter.end()) {
//						cerr<<"Could not find parameter constant: "<<paramName<<" when creating function "<<funcName<<". Quitting"<<endl;
//						return false;
//					}
//					paramNames.push_back(paramName);
//					paramValues.push_back(parameter.find(paramName)->second);
//				}
//			}
//
//
//			//Read in the actual function definition
//			TiXmlElement *pDefinition = pFunction->FirstChildElement("Definition");
//			if(pDefinition) {
//				if(!pDefinition->GetText()) {
//					cerr<<"!!!Error:  Definition tag must actually contain a string for the function.  Quitting."<<endl;
//										return false;
//				}
//				string functionDefintion = pDefinition->GetText();
//				GlobalFunction *gf = new GlobalFunction(funcName, functionDefintion, argNames, argTypes, paramNames, paramValues);
//				if(!system->addGlobalFunction(gf)) {
//					cerr<<"!!!Error:  Function name '"<<funcName<<"' has already been used.  You can't have two\n";
//					cerr<<"functions with the same name, so I'll just stop now."<<endl;
//					return false;
//				}
//				if(verbose) cout<<"\t\t\t= "<<functionDefintion<<endl;
//			} else {
//				cerr<<"!!!Error:  Definition tag for a function must exist!  Quitting."<<endl;
//				return false;
//			}
//
//			argNames.clear();
//			argTypes.clear();
//			paramNames.clear();
//			paramValues.clear();
//		}
//
//
//
//		cout<<"done reading functions!"<<endl;
//		exit(0);
//		//Getting here means we read everything we could successfully
//		return true;
//	} catch (...) {
//		//Uh oh! we got some unknown exception thrown, so we must abort!
//		cerr<<"I caught some unknown error when I was trying to parse out a Global Function.\n";
//		cerr<<"I'm at a loss for words right now, so you're on you're own."<<endl;
//		return false;
//	}
//}
