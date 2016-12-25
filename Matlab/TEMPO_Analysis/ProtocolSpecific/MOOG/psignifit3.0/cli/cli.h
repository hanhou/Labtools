#ifndef CLI_H
#define CLI_H

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <string.h>
#include <cstdlib>
#include "cli_version.h"

class cli_parser {
	private:
		std::map<std::string,std::string> optargs;
		std::map<std::string,std::string> optargs_help;
		std::map<std::string,bool>        optswitch;
		std::map<std::string,std::string> optswitch_help;
		std::list<std::string>            args;
		std::string                       command;
	public:
		cli_parser ( std::string call ) : command ( call ) {}
		void add_option ( std::string option, std::string helptext, std::string defaultvalue="" ) {
			optargs[option]   = defaultvalue; optargs_help[option] = helptext; }
		void add_switch ( std::string option, std::string helptext, bool defaultvalue=false ) {
			optswitch[option] = defaultvalue; optswitch_help[option] = helptext; }
		void parse_args ( int argc, char ** argv );
		std::string getOptArg ( std::string option ) { return optargs[option]; }
		bool        getOptSet ( std::string option ) { return optswitch[option]; }
		std::string popArg ( void ) { if ( args.empty() ) return ""; std::string out ( args.front() ); args.pop_front(); return out; }
		void print_help ( void );
        void print_version ( void );
};

#endif
