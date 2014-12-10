/// Returns data (generic type) corresponding to keyword
/**
 * @param t data (generic type) (output)
 * @param key keyword
 */
template <class T>
int setupAide::getArgs(T& t, string key){
  matrix<T> m;

  if( !getArgs(m,key) )
    return 0;

  if(m.size())
    t = m[1];

  return 1;
}

/// Returns data (generic type matrix) corresponding to keyword
/**
 * @param m formatted data (generic type matrix) corresponding to input keyword
 * @param key keyword
 */
template <class T>
int setupAide::getArgs(matrix<T>& m, string key){
  stringstream args;
  vector<T> argv;
  int argc;
  T input;

  args.str( getArgs(key) );

  while(args >> input)
    argv.push_back(input);

  argc = argv.size();

  if(!argc)
    return 0;

  m.resize(argc,1);

  for(int i=1; i<=argc; i++)
    m[i] = argv[i-1];

  return 1;
}

template <class T>
int setupAide::modifyArgs(T &t, string key){
  stringstream ss;
  ss << t;

  string &d         = data[key];
  const int existed = !d.empty();

  d = ss.str();

  return existed;
}

template <class T>
int setupAide::modifyArgs(matrix<T> &m, string key){
  stringstream ss;

  for(int i = 1; i <= m.size(); i++)
    ss << m[i] << " ";

  string &d         = data[key];
  const int existed = !d.empty();

  d = ss.str();

  return existed;
}
