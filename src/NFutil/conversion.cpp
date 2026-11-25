#include "NFutil.hh"
#include <fstream>
#include <iostream>
#include <algorithm>


using namespace NFutil;


namespace {

// One line of a data file, whichever line end it has: "\n", "\r\n" or "\r".
bool getDataLine(std::istream &in, std::string &line)
{
	line.clear();
	bool any = false;
	for (int c = in.get(); c != std::char_traits<char>::eof(); c = in.get()) {
		any = true;
		if (c == '\n') return true;
		if (c == '\r') {
			if (in.peek() == '\n') in.get();
			return true;
		}
		line += static_cast<char>(c);
	}
	return any;
}

} // namespace


TimeSeries NFutil::loadTimeSeries(const std::string& filePath, const std::string& callerName)
{
	TimeSeries ts;
	std::ifstream file(filePath.c_str());

	if (!file.good()) {
		throw std::runtime_error("File doesn't look like it exists: " + filePath);
	}

	try {
		std::string line;
		int lineNumber = 0;
		bool hasDirection = false;
		bool isIncreasing = false;
		double prevTime = 0.0;
		bool first = true;

		// One point per line. A '#' starts a comment that runs to the end of its
		// line: BioNetGen's .tfun files open with a '#' header that names the two
		// columns, and may carry other '#' lines. A line with nothing else on it
		// is skipped.
		while (getDataLine(file, line)) {
			++lineNumber;
			if (lineNumber == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) {
				line.erase(0, 3);  // a UTF-8 byte order mark
			}
			line.erase(std::min(line.find('#'), line.size()));
			std::istringstream row(line);
			std::string a, b, extra;
			if (!(row >> a)) continue;
			double t = 0.0, v = 0.0;
			try {
				if (!(row >> b) || (row >> extra)) throw std::runtime_error("");
				t = NFutil::convertToDouble(a);
				v = NFutil::convertToDouble(b);
			} catch (std::runtime_error const &) {
				NFutil::trim(line);
				throw std::runtime_error("Data file " + filePath + ", line " + NFutil::toString(lineNumber) +
					": expected two numbers, found '" + line + "'");
			}
			ts.time.push_back(t);
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
