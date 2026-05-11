#include "NFutil.hh"
#include <fstream>
#include <iostream>


using namespace NFutil;


TimeSeries NFutil::loadTimeSeries(const std::string& filePath, const std::string& callerName)
{
	TimeSeries ts;
	std::ifstream file(filePath.c_str());

	if (!file.good()) {
		throw std::runtime_error("File doesn't look like it exists: " + filePath);
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
					throw std::runtime_error("Time values in data file must be strictly monotonic. Found duplicate time: " + NFutil::toString(t));
				}

				if (!hasDirection) {
					isIncreasing = (t > prevTime);
					hasDirection = true;
				} else {
					if ((isIncreasing && t < prevTime) || (!isIncreasing && t > prevTime)) {
						throw std::runtime_error("Time values in data file must be strictly monotonic.");
					}
				}
				prevTime = t;
			}
		}

		if (ts.time.size() == 0) {
			throw std::runtime_error("Data file is empty or invalid format.");
		}
	} catch (std::runtime_error const & e) {
		// Re-throw our specifically constructed runtime_errors without wrapping them further
		throw;
	} catch (std::exception const & e) {
		throw std::runtime_error("Failed to either open or read the file, or invalid number format.\n" + std::string(e.what()));
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
