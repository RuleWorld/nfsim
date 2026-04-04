#include "NFutil.hh"
#include <fstream>
#include <iostream>


using namespace NFutil;


TimeSeries NFutil::loadTimeSeries(const std::string& filePath, const std::string& callerName)
{
	TimeSeries ts;
	std::ifstream file(filePath.c_str());

	if (!file.good()) {
		std::cerr << "Error preparing function " << callerName << "!!" << std::endl;
		std::cerr << "File doesn't look like it exists" << std::endl;
		std::cerr << "Quitting." << std::endl;
		exit(1);
	}

	try {
		std::string a, b;
		bool hasDirection = false;
		bool isIncreasing = false;
		double prevTime = 0.0;
		bool first = true;

		while (file >> a >> b) {
			double t = NFutil::convertToDouble(a);
			ts.time.push_back(t);

			double v = NFutil::convertToDouble(b);
			ts.values.push_back(v);

			if (first) {
				prevTime = t;
				first = false;
			} else {
				if (t == prevTime) {
					std::cerr << "Error in function " << callerName << "!!" << std::endl;
					std::cerr << "Time values in data file must be strictly monotonic. Found duplicate time: " << t << std::endl;
					std::cerr << "Quitting." << std::endl;
					exit(1);
				}

				if (!hasDirection) {
					isIncreasing = (t > prevTime);
					hasDirection = true;
				} else {
					if ((isIncreasing && t < prevTime) || (!isIncreasing && t > prevTime)) {
						std::cerr << "Error in function " << callerName << "!!" << std::endl;
						std::cerr << "Time values in data file must be strictly monotonic." << std::endl;
						std::cerr << "Quitting." << std::endl;
						exit(1);
					}
				}
				prevTime = t;
			}
		}

		if (ts.time.size() == 0) {
			std::cerr << "Error in function " << callerName << "!!" << std::endl;
			std::cerr << "Data file is empty or invalid format." << std::endl;
			std::cerr << "Quitting." << std::endl;
			exit(1);
		}
	} catch (std::exception const & e) {
		std::cerr << "Error preparing function " << callerName << "!!" << std::endl;
		std::cerr << "Failed to either open or read the file, or invalid number format." << std::endl;
		std::cerr << e.what() << std::endl;
		std::cerr << "Quitting." << std::endl;
		exit(1);
	}

	return ts;
}

double NFutil::convertToDouble(const std::string& s)
{
	bool failIfLeftoverChars = true;
	std::istringstream i(s);
	double x;
	char c;
	if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
		throw std::runtime_error("error in NFutil::convertToDouble(\"" + s + "\")");
	return x;
}
int NFutil::convertToInt(const std::string& s)
{
	bool failIfLeftoverChars = true;
	std::istringstream i(s);
	int x;
	char c;
	if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
		throw std::runtime_error("error in NFutil::convertToInt(\"" + s + "\")");
	return x;
}


string NFutil::toString(double x)
{
	std::ostringstream o;
	if (!(o << x)) {
		cout<<endl; cerr<<"Error converting double to string."<<endl;
		exit(1);
	}
	return o.str();
}
string NFutil::toString(int x)
{
	std::ostringstream o;
	if (!(o << x)) {
		cout<<endl; cerr<<"Error converting double to string."<<endl;
		exit(1);
	}
	return o.str();
}
