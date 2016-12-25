#ifndef TESTING_H
#define TESTING_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <string>
#include <cmath>

class TestSuite
{
	private:
		std::vector< int (*)(TestSuite*) > tests;
		std::fstream testlog;
		std::vector<std::string> names;
	public:
		TestSuite ( const char* logfile ) { testlog.open( logfile, std::ios::out ); }
		~TestSuite () { testlog.close(); }
		void addTest ( int (*test)(TestSuite*), const char* testname ) { tests.push_back ( test ); names.push_back ( std::string(testname) ); }
		int runTests ( void ) {
			int failed(0);
			unsigned int i;
			// std::list< int (*)(TestSuite*) > i;
			// for ( i=tests.begin(); i!=tests.end(); ++i )
			// 	failed += (*i)(this);
			for ( i=0; i<tests.size(); i++ ) {
				std::clog << "=============> " << names[i] << " <============\n";
				failed += tests[i](this);
			}
			std::clog << "\n==> " << failed << " tests failed\n";
			return failed;
		}
		int isequal ( double x, double y, const char* testname, double accuracy=1e-7 ) {
			if ( fabs(x-y)<accuracy ) {
				std::clog << "[ OK ]   " << testname << " value: " << x << "\n";
				return 0;
			} else {
				std::clog << "[FAIL]   " << testname << " value: " << x << " should be: " << y << "\n";
				testlog << testname << " value: " << x << " should be: " << y << "\n";
				return 1;
			}
		}
		int isequal_rel ( double x, double y, const char* testname, double accuracy=1e-7 ) {
			if ( fabs ( log10 ( x/y ) ) < accuracy ) {
				std::clog << "[ OK ]   " << testname << " value: " << x << "\n";
				return 0;
			} else {
				std::clog << "[FAIL]   " << testname << " value: " << x << " should be: " << y << "\n";
				testlog << testname << " value: " << x << " should be: " << y << "\n";
				return 1;
			}
		}
		int isless ( double x, double y, const char* testname ) {
			if ( x<y ) {
				std::clog << "[ OK ]   " << testname << " value: " << x << " < " << y << "\n";
				return 0;
			} else {
				std::clog << "[FAIL]   " << testname << " value: " << x << " should be less than " << y << "\n";
				testlog << testname << " value: " << x << " should be less than " << y << "\n";
				return 1;
			}
		}
		int ismore ( double x, double y, const char* testname ) {
			if ( x>y ) {
				std::clog << "[ OK ]   " << testname << " value: " << x << " > " << y << "\n";
				return 0;
			} else {
				std::clog << "[FAIL]   " << testname << " value: " << x << " should be less than " << y << "\n";
				testlog << testname << " value: " << x << " should be less than " << y << "\n";
				return 1;
			}
		}

		int conditional ( bool condition, const char * testname ) {
			if ( condition ) {
				std::clog << "[ OK ]   " << testname << " good \n";
				return 0;
			} else {
				std::clog << "[FAIL]   " << testname << "\n";
				testlog << testname << "\n";
				return 1;
			}
		}
};

#endif
