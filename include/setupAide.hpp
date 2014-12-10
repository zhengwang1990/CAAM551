#ifndef RASCALS_SETUPAIDE
#define RASCALS_SETUPAIDE

class setupAide;

#include "projectHeader.hpp"

#include "matrix.hpp"

class setupAide {
private:
  map<string,string> data;

public:
  setupAide();
  setupAide(string setupFile);

  setupAide(const setupAide &sa);
  setupAide& operator = (const setupAide &sa);

  string readFile(string setupFile);
  void read(string setupFile);

  inline string operator [] (string key){
    return data[key];
  }

  string getArgs(string key);

  template <class T>
  int getArgs(T &t, string key);

  template <class T>
  int getArgs(matrix<T> &m, string key);

  template <class T>
  int modifyArgs(T &t, string key);

  template <class T>
  int modifyArgs(matrix<T> &m, string key);

  int getArgs(matrix<string> &m, string key, string delimeter);

  void append(string setupFile);
  void append(const setupAide &sa);
};

#include "setupAide.tpp"

#endif
