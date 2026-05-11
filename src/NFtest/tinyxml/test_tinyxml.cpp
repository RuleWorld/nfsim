#include "test_tinyxml.hh"
#include "../../NFinput/TinyXML/tinyxml.h"
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

void NFtest_tinyxml::run()
{
	cout << "Running TinyXML tests..." << endl;

	// Test case 1: Parse a known good XML string
	cout << "  Testing known good XML string..." << endl;
	const char* good_xml =
		"<?xml version=\"1.0\"?>"
		"<Model>"
		"  <ListOfParameters>"
		"    <Parameter id=\"k1\" value=\"1.5\"/>"
		"  </ListOfParameters>"
		"</Model>";

	TiXmlDocument doc_good;
	doc_good.Parse(good_xml);
	if (doc_good.Error()) {
		throw std::runtime_error("TinyXML failed to parse known good XML string: " + std::string(doc_good.ErrorDesc()));
	}

	TiXmlElement* root = doc_good.RootElement();
	if (!root || std::string(root->Value()) != "Model") {
		throw std::runtime_error("TinyXML parsed good XML but root element is missing or incorrect");
	}

	// Test case 2: Parse a known bad XML string
	cout << "  Testing known bad XML string..." << endl;
	const char* bad_xml =
		"<?xml version=\"1.0\"?>"
		"<Model>"
		"  <ListOfParameters>"
		"    <Parameter id=\"k1\" value=\"1.5\""
		"  </ListOfParameters>"
		"</Model>";

	TiXmlDocument doc_bad;
	doc_bad.Parse(bad_xml);
	if (!doc_bad.Error()) {
		throw std::runtime_error("TinyXML successfully parsed a known bad XML string (it should have failed)");
	}

	cout << "  TinyXML tests passed!" << endl;
}
