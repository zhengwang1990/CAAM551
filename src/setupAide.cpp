#include "setupAide.hpp"

/// Default constructor
setupAide::setupAide(){}

/// Constructor
/**
 * @param setupFile file path
 */
setupAide::setupAide(string setupFile){
  read(setupFile);
}

/// Copy constructor
setupAide::setupAide(const setupAide& sa){
  data = sa.data;
}

/// Copies parameters of another setupAide
setupAide& setupAide::operator = (const setupAide& sa){
  data = sa.data;
  return *this;
}

/// Load file into string
/**
 * @param filename file path
 * @return string containing file data
 */
string setupAide::readFile(string filename){
  struct stat statbuf;

  FILE *fh = fopen(filename.c_str(), "r");
  if (!fh){
    printf("Failed to open: %s\n", filename.c_str());
    throw 1;
  }

  stat(filename.c_str(), &statbuf);
  char *source = (char *) malloc(statbuf.st_size + 1);
  size_t countCheck = fread(source, statbuf.st_size, 1, fh);
  if(!countCheck) {
    printf("Failed to read: %s\n", filename.c_str());
    throw 1;
  }
  source[statbuf.st_size] = '\0';

  string ret = source;

  return ret;
}

/// Reads and interprets data from file
/**
 * @param setupFile file path
 */
void setupAide::read(string setupFile){
  string args = readFile(setupFile);

  int size = args.length();
  string current = "", key = "";
  stringstream ss;
  char c;

  string whitespace = " \n\t\v\f\r";

  for(int i=0; i<size; i++){
    c = args[i];

    // Batch strings together
    if(c == '\'' || c == '"'){
      current += c;
      i++;

      while(i < size && args[i] != c)
	current += args[i++];

      if(i >= size)
	break;

      if( i < (size-1) )
	current += args[i];
    }

    // Batch comments
    else if(c == '/' && i < size && args[i+1] == '*'){
      i += 2;

      while( args[i] != '*' || (i < size && args[i+1] != '/') )
	i++;

      if(i >= size)
	break;

      i++;
    }

    // Removing # comments
    else if(c == '#'){
      i++;

      while(i < size && args[i] != '\n')
	i++;
    }

    // Change \[\] to []
    else if(c == '\\' && i < size && (args[i+1] == '[' || args[i+1] == ']')){
      current += args[i+1];
      i += 2;
    }

    // Split keywords []
    else if(c == '['){
      if(!key.empty()){
        int left  = 0;
        int right = current.size();

        while( (left < current.size()) && (whitespace.find( current[left] ) != string::npos) )
          left++;

        while( (0 <= right) && (whitespace.find( current[right - 1] ) != string::npos) )
          right--;

        data[key] = current.substr(left, right - left);
      }

      current = "";
      key     = "";

      i++;

      while(i < size && args[i] != ']')
	key += args[i++];
    }

    // Else add the character
    else
      current += c;

    if(i >= (size-1) && current.length())
      data[key] = current;
  }
}

/// Returns data (string) corresponding to keyword
/**
 * @param key keyword
 * @return data (formatted as string) corresponding to keyword
 */
string setupAide::getArgs(string key){
  return data[key];
}

/// Returns data (string matrix) corresponding to keyword
/**
 * @param m formatted data (string matrix) corresponding to input keyword
 * @param key keyword
 * @param delimeter delimeter to search for and separate entries in data
 */
int setupAide::getArgs(matrix<string>& m, string key, string delimeter){
  string args, current;
  vector<string> argv;
  int argc, size;

  args = data[key];

  size = args.length();

  current = "";

  /// Iterates through all arguments in args, pushing each (separated by a delimeter) onto argv
  for(int i=0; i<size; i++) {

    while( i < size && delimeter.find(args[i]) == string::npos )
      current += args[i++];

    if(current.length())
      argv.push_back(current);

    current = "";
  }

  argc = argv.size();

  /// If there are no arguments, return
  if(!argc)
    return 0;

  m.resize(argc,1);

  for(int i=1; i<=argc; i++)
    m[i] = argv[i-1];

  return 1;
}

void setupAide::append(string setupFile){
  append( setupAide(setupFile) );
}

void setupAide::append(const setupAide &sa){
  for(std::map<string,string>::const_iterator it = sa.data.begin(); it != sa.data.end(); it++)
    data[(string) it->first] = (string) it->second;
}
