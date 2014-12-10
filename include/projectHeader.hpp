#ifndef RASCALS_PROJECT
#define RASCALS_PROJECT

#include <sys/types.h>
#include <sys/stat.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#include <algorithm>
#include <iostream>
#include <ostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <map>
#include <typeinfo>

using namespace std;

// Common matrix types
#define fmatrix matrix<datafloat>
#define imatrix matrix<int>
#define smatrix matrix<string>

#endif
