#include "NFfunction.hh"

#include <cctype>
#include <math.h>

using namespace NFcore;
#ifndef NFSIM_USE_EXPRTK
using namespace mu;
#endif

namespace {

bool isIdentifierCharacter(char c)
{
	return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}

bool isTimeAliasAt(const std::string &expression, size_t start,
		const std::string &alias, size_t &end)
{
	if (expression.compare(start, alias.size(), alias) != 0) return false;
	if (start > 0 && isIdentifierCharacter(expression[start - 1])) return false;

	size_t cursor = start + alias.size();
	if (cursor < expression.size() && isIdentifierCharacter(expression[cursor])) {
		return false;
	}
	while (cursor < expression.size() &&
			std::isspace(static_cast<unsigned char>(expression[cursor]))) {
		++cursor;
	}
	if (cursor >= expression.size() || expression[cursor] != '(') return false;
	++cursor;
	while (cursor < expression.size() &&
			std::isspace(static_cast<unsigned char>(expression[cursor]))) {
		++cursor;
	}
	if (cursor >= expression.size() || expression[cursor] != ')') return false;

	end = cursor + 1;
	return true;
}

} // namespace

std::string normalizeTimeExpression(const std::string &expression)
{
	std::string normalized;
	normalized.reserve(expression.size());

	for (size_t i = 0; i < expression.size();) {
		size_t end = i;
		if (isTimeAliasAt(expression, i, "time", end) ||
				isTimeAliasAt(expression, i, "t", end)) {
			normalized += "time";
			i = end;
		} else {
			normalized += expression[i++];
		}
	}

	return normalized;
}

bool containsTimeExpression(const std::string &expression)
{
	return normalizeTimeExpression(expression) != expression;
}




mu::Parser * FuncFactory::create(string functionString, vector <string> & variableNames, vector <double *> & variablePtrs)
{
	double PI = 3.14159265358979323846;
	double NA= 6.02214179e23;
	double E = 2.718281828459;

	mu::Parser *p = new mu::Parser();
	try
	{
		p->DefineConst("_PI",PI);
		p->DefineConst("_e",E);
		p->DefineConst("_Na",NA);

		if(variableNames.size()!=variablePtrs.size())  {
			cout<<"Error parsing function in FuncFactory!!  Your variableNames vector and "<<endl;
			cout<<"you variablePtrs vector do not appear to match up!"<<endl;
			cout<<"The function you gave me was: "<<functionString<<endl;
			cout<<"For that, I will quit"<<endl;
			delete p;
			throw std::runtime_error("FuncFactory::create: mismatched variableNames and variablePtrs");
		} else {
			for(unsigned int v=0;v<variableNames.size(); v++) {
				p->DefineVar(variableNames.at(v), variablePtrs.at(v));
			}
		}

		p->SetExpr(functionString);
	}
	catch (mu::Parser::exception_type &e)
	{
		cout<<"Error parsing function in FuncFactory!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		delete p;
		throw std::runtime_error("FuncFactory::create: " + e.GetMsg());
	}

	return p;
}

mu::Parser * FuncFactory::create(bool throw_mock_exception)
{
	double PI = 3.14159265358979323846;
	double NA= 6.02214179e23;
	double E = 2.718281828459;
	mu::Parser *p = new mu::Parser();
	try
	{
		p->DefineConst("_PI",PI);
		p->DefineConst("_e",E);
		p->DefineConst("_Na",NA);
		if (throw_mock_exception) {
			throw mu::Parser::exception_type("Mock exception");
		}
	}
	catch (mu::Parser::exception_type &e)
	{
		cout<<"Error creating function in FuncFactory!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		delete p;
		throw std::runtime_error("FuncFactory::create: " + e.GetMsg());
	}
	return p;
}




double FuncFactory::Eval(mu::Parser *p)
{
	if(p==NULL) {
		cout<<"In FuncFactory: Trying to evaluate a null Parser! You are probably trying\n";
		cout<<"to use a GlobalFunction before it has been prepared! Preparing a GlobalFunction\n";
		cout<<"connects it to Observables, so it must be done before you can use it!\n";
		cout<<"  we've all made this mistake before, but now I'm exiting..."<<endl;
		throw std::runtime_error("FuncFactory::Eval: p is NULL");
	}
	try {
		return p->Eval();
	} catch (mu::Parser::exception_type &e) {
		cout<<"Error evaluating function in FuncFactory!!  "<<endl;
		cout<<"The function was: "<<p->GetExpr()<<endl;
		cout<<"And this is what went wrong:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Terminating your simulation. Better luck next time."<<endl;
		throw std::runtime_error("FuncFactory::Eval: " + e.GetMsg());
	}
	return 0;
}
static double MockThrowException() {
	throw mu::Parser::exception_type("Mock Eval exception");
	return 0;
}


void FuncFactory::test()
{
	double PI = 3.14159265358979323846;
	double NA= 6.02214179e23;
	double E = 2.718281828459;
	cout<<"Beginning diagnostic tests..."<<endl;

	{
	//Test 1: check constants are correct and basic math functions seem right
	cout<<" 1) simple test of constants and predefined functions: ";
	string functionString("sin(_e*cos(3.2/_PI))+ln(_Na*1.14e-11)");//sin(_e*cos(3.2/_PI)+0.11*(1-_Na*1.14)");
	vector <string> variableNames;
	vector <double *> variablePtrs;

	mu::Parser *p = FuncFactory::create(functionString,variableNames,variablePtrs);
	double result = sin(E*cos(3.2/PI))+log(NA*1.14e-11);
	double funcResult = p->Eval();
	if(abs(funcResult - result)<0.0001)
		cout<<"pass."<<endl;
	else
		cout<<"fail! p->Eval() = "<<funcResult<<"  but should be: "<<result<<endl;
	delete p;
	}


	{
	//Test 2a: check that variables are working properly
	cout<<" 2a) test that variable input works properly: ";

	string functionString = "1-(d1/d2)*sin(d1*d2)+1.3*d1^2";

	string d1_name("d1");
	string d2_name("d2");
	vector <string> variableNames;
	variableNames.push_back(d1_name);
	variableNames.push_back(d2_name);

	double d1 = -1.412;
	double d2 = 2.01e1;
	vector <double *> variablePtrs;
	variablePtrs.push_back(&d1);
	variablePtrs.push_back(&d2);

	mu::Parser *p = FuncFactory::create(functionString,variableNames,variablePtrs);
	double result = 1-(d1/d2)*sin(d1*d2)+1.3*pow(d1,2.0);
	double funcResult = FuncFactory::Eval(p);
	if(abs(funcResult - result)<0.0001)
		cout<<"pass."<<endl;
	else
		cout<<"fail! p->Eval() = "<<funcResult<<"  but should be: "<<result<<endl;


	//Test 2b: check that variables are working properly
	cout<<" 2b) another check of the variables: ";
	d1 = 50.1;
	d2 = 3.1;
	result = 1-(d1/d2)*sin(d1*d2)+1.3*pow(d1,2.0);
	funcResult = FuncFactory::Eval(p);
	if(abs(funcResult - result)<0.0001)
		cout<<"pass."<<endl;
	else
		cout<<"fail! p->Eval() = "<<funcResult<<"  but should be: "<<result<<endl;

	delete p;
	}

	{
	//Test 3: Check error path for parameterless create
	cout<<" 3) test parameterless create() error path: ";
	bool threw = false;
	try {
		FuncFactory::create(true);
	} catch(const std::runtime_error& e) {
		threw = true;
	}
	if(threw)
		cout<<"pass."<<endl;
	else {
		cout<<"fail! Exception not thrown."<<endl;
		exit(1);
	}
	}

	{
	//Test 4: Check error path for mismatched vectors in create
	cout<<" 4) test create() error path with mismatched vectors: ";
	string functionString("sin(d1)");
	vector <string> variableNames;
	variableNames.push_back("d1");
	vector <double *> variablePtrs; // Empty, so mismatch
	bool threw = false;
	try {
		FuncFactory::create(functionString, variableNames, variablePtrs);
	} catch(const std::runtime_error& e) {
		threw = true;
	}
	if(threw)
		cout<<"pass."<<endl;
	else {
		cout<<"fail! Exception not thrown."<<endl;
		exit(1);
	}
	}

	{
	//Test 5: Check error path for bad function string evaluation
	cout<<" 5) test Eval() error path with undefined variable string: ";
	string functionString("sin(d1)"); // undefined variable
	vector <string> variableNames;
	vector <double *> variablePtrs;
	bool threw = false;
	try {
		mu::Parser* bad_p = FuncFactory::create(functionString, variableNames, variablePtrs);
		FuncFactory::Eval(bad_p);
	} catch(const std::runtime_error& e) {
		threw = true;
	}
	if(threw)
		cout<<"pass."<<endl;
	else {
		cout<<"fail! Exception not thrown."<<endl;
		exit(1);
	}
	}

	{
	//Test 6: Check error path for Eval with NULL
	cout<<" 6) test Eval() error path with NULL: ";
	bool threw = false;
	try {
		FuncFactory::Eval(NULL);
	} catch(const std::runtime_error& e) {
		threw = true;
	}
	if(threw)
		cout<<"pass."<<endl;
	else {
		cout<<"fail! Exception not thrown."<<endl;
		exit(1);
	}
	}

	{
	//Test 7: Check error path for Eval with mock exception
	cout<<" 7) test Eval() error path with mock exception: ";
	mu::Parser *p = FuncFactory::create();
	p->DefineFun("mock_throw", MockThrowException);
	p->SetExpr("mock_throw()");
	bool threw = false;
	try {
		FuncFactory::Eval(p);
	} catch(const std::runtime_error& e) {
		threw = true;
	}
	if(threw)
		cout<<"pass."<<endl;
	else {
		cout<<"fail! Exception not thrown."<<endl;
		exit(1);
	}
	delete p;
	}

	// Tests 8 to 10 hold for either parser: muParser, the default, and ExprTk
	// (-DNFSIM_USE_EXPRTK=ON) must give the same answers.
	{
	//Test 8: rint rounds a half up, floor(x + 0.5), as muParser and BioNetGen do
	cout<<" 8) test rint() rounds a half up: ";
	const char *expressions[] = {"rint(-2.5)", "rint(-1.5)", "rint(-0.5)", "rint(0.5)", "rint(2.5)"};
	const double expected[] = {-2, -1, 0, 1, 3};
	for(int i=0; i<5; i++) {
		vector <string> variableNames;
		vector <double *> variablePtrs;
		mu::Parser *p = FuncFactory::create(expressions[i],variableNames,variablePtrs);
		double value = FuncFactory::Eval(p);
		delete p;
		if(value != expected[i]) {
			cout<<"fail! "<<expressions[i]<<" = "<<value<<"  but should be: "<<expected[i]<<endl;
			exit(1);
		}
	}
	cout<<"pass."<<endl;
	}

	{
	//Test 9: a literal is read as the double nearest it, as strtod (and the
	//compiler, for the literals below) reads it
	cout<<" 9) test numeric literals are read as strtod reads them: ";
	const char *expressions[] = {"7.77e-05", "0.49999999999999997", "8.2e-02", "2.638e-03",
		"7.6e-07", "5.363e+11", "0.44982514983446460", "0.45723577511449437", "6.23286e-20",
		"2^0.49999999999999997", "-7.77e-05*2"};
	const double expected[] = {7.77e-05, 0.49999999999999997, 8.2e-02, 2.638e-03,
		7.6e-07, 5.363e+11, 0.44982514983446460, 0.45723577511449437, 6.23286e-20,
		pow(2.0, 0.49999999999999997), -7.77e-05*2};
	for(int i=0; i<11; i++) {
		vector <string> variableNames;
		vector <double *> variablePtrs;
		mu::Parser *p = FuncFactory::create(expressions[i],variableNames,variablePtrs);
		double value = FuncFactory::Eval(p);
		delete p;
		if(value != expected[i]) {
			cout.precision(17);
			cout<<"fail! "<<expressions[i]<<" was read as "<<value<<"  but should be: "<<expected[i]<<endl;
			exit(1);
		}
	}
	cout<<"pass."<<endl;
	}

	{
	//Test 10: expressions both parsers refuse. ExprTk read each of the first
	//four as two statements and kept the last one's value: 3*exp(0)(2) gave 2.
	//The fifth puts two numbers side by side, one of them a literal ExprTk
	//misreads, which the shim respells. The last two have a stray ')' and a
	//missing '(', which balance once the expression is put in brackets.
	cout<<" 10) test expressions muParser refuses are refused: ";
	const char *expressions[] = {"3*exp(0)(2)", "3*rint(1.2) (2)", "3*abs(-1) if(1,5,6)",
		"3*exp(0);2", "2 0.45723577511449437", "3*exp(0));(2", "1)(2"};
	for(int i=0; i<7; i++) {
		vector <string> variableNames;
		vector <double *> variablePtrs;
		bool threw = false;
		double value = 0;
		mu::Parser *p = NULL;
		try {
			p = FuncFactory::create(expressions[i],variableNames,variablePtrs);
			value = FuncFactory::Eval(p);
		} catch(const std::runtime_error& e) {
			threw = true;
		}
		delete p;
		if(!threw) {
			cout<<"fail! "<<expressions[i]<<" was not refused; it gave "<<value<<endl;
			exit(1);
		}
	}
	cout<<"pass."<<endl;
	}

	//Thats all the test I can think of!
	cout<<endl<<"Testing complete."<<endl;
}
