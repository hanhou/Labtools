#ifndef CLI_UTILITIES_H
#define CLI_UTILITIES_H

#include "../src/psipp.h"
#include <string>
#include <fstream>
#include <cstdio>
#include <list>

PsiData * allocateDataFromFile ( std::string fname, int nafc );

PsiPsychometric * allocatePsychometric ( std::string core, std::string sigmoid, int nafc, PsiData* data, bool verbose );

PsiPrior * allocatePrior ( std::string prior );
void setPriors ( PsiPsychometric * pmf, std::string prior1, std::string prior2, std::string prior3, std::string prior4 );

std::vector<double> getCuts ( std::string cuts );

void print ( std::vector<double> theta, bool matlabformat, std::string varname, FILE *ofile );
void print ( double theta, bool matlabformat, std::string varname, FILE *ofile );
void print ( std::vector< std::vector<double> >& theta, bool matlabformat, std::string varname, FILE *ofile );
void print ( std::vector< std::vector<int> >& theta, bool matlabformat, std::string varname, FILE *ofile );
void print ( std::vector<int> theta, bool matlabformat, std::string varname, FILE *ofile );

void print_fisher ( PsiPsychometric *pmf, std::vector<double> theta, PsiData *data, FILE* ofile, bool matlabformat );

#endif
