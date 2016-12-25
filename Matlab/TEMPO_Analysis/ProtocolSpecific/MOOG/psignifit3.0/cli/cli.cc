#include "cli.h"

void cli_parser::parse_args ( int argc, char ** argv ) {
	unsigned int arg;

	for (arg=0; arg<argc; arg++) {
		if ( !strcmp ( argv[arg], "--version" ) ) {
			print_version ();
            exit(0);
		}
	}

	for (arg=0; arg<argc; arg++) {
		if ( !strcmp ( argv[arg], "-h" ) ) {
			print_help ();
		}
	}

	for (arg=1; arg<argc; arg++) {
		if ( optargs.count ( argv[arg] ) > 0 && arg<argc-1 ) {
			optargs[argv[arg]] = argv[arg+1]; arg++;
		} else if ( optswitch.count ( argv[arg] ) > 0 ) {
			optswitch[argv[arg]] = !optswitch[argv[arg]];
		} else {
			args.push_front ( argv[arg] );
		}
	}
}

void cli_parser::print_help ( void ) {
    print_version();
	std::map<std::string,std::string>::iterator optargs_i;
	std::map<std::string,bool>::iterator        optswitch_i;

	std::cout << "\nCall: " << command << "\n\n";
	std::cout << "Options:\n";
	std::cout << "--------\n\n";

	for ( optargs_i = optargs.begin(); optargs_i != optargs.end(); optargs_i++ ) {
		std::cout << "  " << (*optargs_i).first << "  default: " << (*optargs_i).second << "\n";
		std::cout << "        " << optargs_help[(*optargs_i).first] << "\n\n";
	}

	std::cout << "Switches:\n";
	std::cout << "---------\n\n";

	for ( optswitch_i = optswitch.begin(); optswitch_i != optswitch.end(); optswitch_i++ ) {
		std::cout << "  " << (*optswitch_i).first << "  default: " << (*optswitch_i).second << "\n";
		std::cout << "        " << optswitch_help[(*optswitch_i).first] << "\n\n";
	}

	exit ( 0 );
}

void cli_parser::print_version ( void ) {
    std::cout << "Psignifit3 cli, version: " << VERSION << "\n";
}
