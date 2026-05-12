#include "NFfunction.hh"
#include <stdexcept>


using namespace std;
using namespace NFcore;
using namespace mu;

double tfun_interpolate_value(
	const vector<double> &xs,
	const vector<double> &ys,
	const string &method,
	double x)
{
	if (xs.empty()) return 0.0;
	if (xs.size() == 1) return ys.front();

	const bool increasing = xs[1] > xs[0];
	if (increasing) {
		if (x <= xs.front()) return ys.front();
		if (x >= xs.back()) return ys.back();
	} else {
		if (x >= xs.front()) return ys.front();
		if (x <= xs.back()) return ys.back();
	}

	size_t i = 0;
	while ((i + 1) < (xs.size() - 1) &&
		(increasing ? (x >= xs[i + 1]) : (x <= xs[i + 1]))) {
		++i;
	}

	if (method == "step") {
		return ys[i];
	}

	double x0 = xs[i];
	double x1 = xs[i + 1];
	double y0 = ys[i];
	double y1 = ys[i + 1];
	return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}


GlobalFunction::GlobalFunction(string name,
		string funcExpression,
		vector <string> &varRefNames,
		vector <string> &varRefTypes,
		vector <string> &paramNames,
		System *s)
{
	if(varRefNames.size()!=varRefTypes.size()) {
		cerr<<"Trying to create a global function, but your variable reference vectors don't match up in size!"<<endl;
		cerr<<"Quitting!"<<endl;
		exit(1);
	}

	this->name = name;
	this->funcExpression = funcExpression;

	this->n_varRefs=varRefNames.size();
	this->varRefNames = new string[n_varRefs];
	this->varRefTypes = new string[n_varRefs];
	for(unsigned int vr=0; vr<n_varRefs; vr++) {
		this->varRefNames[vr]=varRefNames.at(vr);
		this->varRefTypes[vr]=varRefTypes.at(vr);
	}

	this->n_params=paramNames.size();
	this->paramNames = new string[n_params];
	for(unsigned int i=0; i<n_params; i++) {
		this->paramNames[i]=paramNames.at(i);
	}
	p=0;
	this->sysPtr = NULL;

	// AS-2021
	this->fileFunc = false;
	this->interpolationMethod = "linear";
	this->currInd = 0;
	this->dataLen = 0;
	this->counter = NULL;
	this->ctrType = "";
	this->ctrName = "";
	this->counterParamName = "";
	// AS-2021
}



GlobalFunction::~GlobalFunction()
{
	delete [] varRefNames;
	delete [] varRefTypes;
	delete [] paramNames;
	if(p!=NULL) delete p;
}




void GlobalFunction::prepareForSimulation(System *s)
{
	try {
		p=FuncFactory::create();
		for(unsigned int vr=0; vr<n_varRefs; vr++)
		{
			if(varRefTypes[vr]=="Observable") {
				Observable *obs = s->getObservableByName(varRefNames[vr]);
				if(obs==NULL) {
					cout<<"When creating global function: "<<this->name<<endl<<" could not find the observable: ";
					cout<<varRefNames[vr]<<" of type "<<varRefTypes[vr]<<endl;
					cout<<"Quitting."<<endl;
					exit(1);
				}
				obs->addReferenceToMyself(p);
			} else {
				cout<<"here"<<endl;
				cout<<"Uh oh, an unrecognized argType ("<<varRefTypes[vr]<<") for a function! "<<varRefNames[vr]<<endl;
				cout<<"Try using the type: \"MoleculeObservable\""<<endl;
				cout<<"Quitting because this will give unpredicatable results, or just crash."<<endl;
				exit(1);
			}
		}

		for(unsigned int i=0; i<n_params; i++) {
			p->DefineConst(paramNames[i],s->getParameter(paramNames[i]));
		}

		if (this->fileFunc && !this->ctrName.empty()) {
			p->DefineConst(this->ctrName, 0.0);
		}
		p->SetExpr(this->funcExpression);

	}
	catch (mu::Parser::exception_type &e)
	{
		cout<<"Error preparing function "<<name<<" in class GlobalFunction!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		exit(1);
	}
}

void GlobalFunction::updateParameters(System *s) {
	//cout<<"Updating parameters for function: "<<name<<endl;
	for(unsigned int i=0; i<n_params; i++) {
		p->DefineConst(paramNames[i],s->getParameter(paramNames[i]));
	}

}




void GlobalFunction::attatchRxn(ReactionClass *r)
{
	//unsigned int n_rxns;
	//ReactionClass *rxns;

}





void GlobalFunction::printDetails()
{
	cout<<"Global Function: '"<< this->name << "()'"<<endl;
	cout<<" ="<<funcExpression<<endl;
	cout<<"   -Variable References:"<<endl;
	for(unsigned int vr=0; vr<n_varRefs; vr++) {
		cout<<"         "<<varRefTypes[vr]<<":  "<<varRefNames[vr]<<" = " << ""<<endl;
	}
	cout<<"   -Constant Parameters:"<<endl;
	for(unsigned int i=0; i<n_params; i++) {
		cout<<"         "<<paramNames[i]<<endl;
	}





//	// Get the map with the variables
//	mu::Parser::varmap_type variables = p->GetVar();
//	cout << (int)variables.size() << " variables."<<endl;
//	mu::Parser::varmap_type::const_iterator item = variables.begin();
//	// Query the variables
//	for (; item!=variables.end(); ++item)
//	{
//	  cout << "  Name: " << item->first << " Address: [0x" << item->second << "]  Value: "<< *(item->second)<<"\n";
//	}


	if(p!=0) {
		// AS-2021
		if (this->fileFunc==true) {
			this->fileUpdate();
		}
		// AS-2021
		cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
	}
}

// AS-2021
void GlobalFunction::loadParamFile(string filePath) 
{
	string callerName = this->name + " in class GlobalFunction";
	NFutil::TimeSeries ts = NFutil::loadTimeSeries(filePath, callerName);
	this->data.push_back(ts.time);
	this->data.push_back(ts.values);
}

void GlobalFunction::addCounterPointer(double *counter){
	this->ctrType = "Observable";
	this->counter = counter;
}

void GlobalFunction::setCtrName(string name) {
	this->ctrName = name;
}

void GlobalFunction::setInterpolationMethod(string method) {
	string normalized = method;
	std::transform(normalized.begin(), normalized.end(), normalized.begin(),
		[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	if (normalized.empty()) normalized = "linear";
	if (normalized != "linear" && normalized != "step") {
		cerr<<"Error preparing function "<<name<<" in class GlobalFunction!!"<<endl;
		cerr<<"Unsupported TFUN interpolation method '"<<method<<"'."<<endl;
		cerr<<"Quitting."<<endl;
		exit(1);
	}
	this->interpolationMethod = normalized;
}

void GlobalFunction::setCounterFromTime(System *s) {
	this->addSystemPointer(s);
}

void GlobalFunction::setCounterFromParameter(System *s, string paramName) {
	this->ctrType = "Parameter";
	this->sysPtr = s;
	this->counterParamName = paramName;
}

void GlobalFunction::addSystemPointer(System *s) {
	this->ctrType = "System";
	this->sysPtr = s;
}

void GlobalFunction::enableFileDependency(string filePath, string method) {
	try {
		this->loadParamFile(filePath);
	} catch (exception const & e) {
			throw std::runtime_error("Error preparing function " + name + " in class GlobalFunction!!\n" + std::string(e.what()));
	};
	// we just want to keep a record of this
	this->filePath = filePath;
	// this sets it up so that this function knows it's supposed
	// to be pulling values from a file
	this->fileFunc = true;
	// initialize internal index
	this->currInd = 0;
	// pull data lenght so we can reuse it
	this->dataLen = data[0].size();
	// set interpolation method if specified
	if (!method.empty()) {
		this->setInterpolationMethod(method);
	}
}

void GlobalFunction::enableInlineDependency(
	const vector<double> &xs,
	const vector<double> &ys,
	string method)
{
	this->data.clear();
	this->data.push_back(xs);
	this->data.push_back(ys);
	this->filePath = "<inline>";
	this->fileFunc = true;
	this->setInterpolationMethod(method);
	this->currInd = 0;
	this->dataLen = static_cast<int>(xs.size());
}

double GlobalFunction::getCounterValue() {
	double ctrVal = 0.0;
	if (ctrType == "Observable") {
		if (counter == NULL) {
			cerr<<"Error preparing function "<<name<<" in class GlobalFunction!!"<<endl;
			cerr<<"Observable TFUN counter pointer is null."<<endl;
			cerr<<"Quitting."<<endl;
			exit(1);
		}
		ctrVal = (*counter);
	} else if (ctrType == "System") {
		if (this->sysPtr == NULL) {
			cerr<<"Error preparing function "<<name<<" in class GlobalFunction!!"<<endl;
			cerr<<"System TFUN counter pointer is null."<<endl;
			cerr<<"Quitting."<<endl;
			exit(1);
		}
		ctrVal = this->sysPtr->getCurrentTime();
	} else if (ctrType == "Parameter") {
		if (this->sysPtr == NULL || this->counterParamName.empty()) {
			cerr<<"Error preparing function "<<name<<" in class GlobalFunction!!"<<endl;
			cerr<<"Parameter TFUN counter is not configured."<<endl;
			cerr<<"Quitting."<<endl;
			exit(1);
		}
		ctrVal = this->sysPtr->getParameter(counterParamName);
	} else {
		cerr<<"Error preparing function "<<name<<" in class GlobalFunction!!"<<endl;
		cerr<<"TFUN counter type '"<<ctrType<<"' is not supported."<<endl;
		cerr<<"Quitting."<<endl;
		exit(1);
	}
	return ctrVal;
}
void GlobalFunction::fileUpdate() {
	this->fileUpdate(this->getCounterValue());
}

void GlobalFunction::fileUpdate(double ctrVal) {
	if (data.size() < 2 || data[0].size() == 0) {
		cerr << "Error in function " << this->name << " in class GlobalFunction!!" << endl;
		cerr << "Data for file update is empty or malformed." << endl;
		cerr << "Quitting." << endl;
		exit(1);
	}
	double y = tfun_interpolate_value(data[0], data[1], interpolationMethod, ctrVal);
	p->DefineConst(ctrName, y);
	return;
}
// AS-2021

void GlobalFunction::printDetails(System *s)
{
	cout<<"Global Function: '"<< this->name << "()'"<<endl;
	cout<<" ="<<funcExpression<<endl;
	cout<<"   -Variable References:"<<endl;
	for(unsigned int vr=0; vr<n_varRefs; vr++) {
		cout<<"         "<<varRefTypes[vr]<<":  "<<varRefNames[vr]<<" = " << s->getObservableByName(varRefNames[vr])->getCount()<<endl;
	}
	cout<<"   -Constant Parameters:"<<endl;
	for(unsigned int i=0; i<n_params; i++) {
		cout<<"         "<<paramNames[i]<<" = " << s->getParameter(paramNames[i])<<endl;
	}

	
	if(p!=0) {
		// AS-2021
		if (this->fileFunc==true) {
			cout<<"   Function relies on file: "<<this->filePath<<endl;
			this->fileUpdate();
		}
		// AS-2021
		cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
	}
}





StateCounter::StateCounter(string name, MoleculeType *mt, string stateName) {
	this->name=name;
	this->mt = mt;
	this->stateIndex = mt->getCompIndexFromName(stateName);
	this->value=0;

	//Make sure this is a state we can count on!
	if(!mt->isIntegerComponent(stateName)) {
		//if it is not an integer component state, we must abort because
		//we can not evaluate the sum of a non integer component
		cerr<<"Trying to create a stateCounter: '"<<name<<"' on the state: '"<<stateName<<"'\n";
		cerr<<"of MoleculeType: '"<<mt->getName()<<"', but the state you have selected cannot\n";
		cerr<<"have integer values, so summations on this state are undefined.  I am quitting."<<endl;
		exit(1);
	}
}
StateCounter::~StateCounter() {
	mt=0;
}

void StateCounter::add(Molecule *m) {
	if(m->getMoleculeType()==mt) {

	//	cout<<"matched moleculeType"<<endl;
		value+=m->getComponentState(stateIndex);
	//	cout<<"found component state: "<<m->getComponentState(stateIndex)<<endl;
	//	cout<<"updating v`alue to: "<< value<<endl;
	}
}
