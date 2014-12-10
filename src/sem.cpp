#include "sem.hpp"
#define NThread 10

sem::sem(){ 
}

sem::sem(int _K, int _N){
  K = _K;
  N = _N;
  Ne1 = K;
  Ne2 = Ne1 * Ne1;
  Np1 = N + 1;
  Np2 = Np1 * Np1;
  h = 2.0 / Ne1;
  J = h * h / 4;
}

void sem::init(){
  gNumbering();
  assembleMatrices();  

  cout << "degree of freedom = " << Ni << endl;

  b.resize(Ni);
  uSoln.resize(Ni);
  aAul.resize(Ne2*Np2);

  setupRHS();
  
  iterMax = 10000;
  tol = 1e-3;
}

void sem::gNumbering(){
  nodeId.resize(Ne2, Np2);

  Ni = 0;
  int e = 0, v, enew, vnew;
  for(int ei=1; ei<=Ne1; ei++){
    for(int ej=1; ej<=Ne1; ej++){
      e++;
      v = 0;
      for(int vi=1; vi<=Np1; vi++){
   	for(int vj=1; vj<=Np1; vj++){
   	  v++;
   	  if (((vi==1)&&(ei==1))||((vj==1)&&(ej==1))||((vi==Np1)&&(ei==Ne1))||((vj==Np1)&&(ej==Ne1))){
   	    nodeId(e,v) = 0;
   	    continue;
   	  }
   	  if ((vi==1)&&(ei!=1)){
   	    enew = e - Ne1;
   	    vnew = N*Np1 + vj;
   	    nodeId(e,v) = nodeId(enew,vnew);
   	    continue;
   	  }
   	  if ((vj==1)&&(ej!=1)){
   	    enew = e - 1;
   	    vnew = (vi-1)*Np1 + Np1;
   	    nodeId(e,v) = nodeId(enew,vnew);
   	    continue;
   	  }
   	  Ni++;
   	  nodeId(e,v) = Ni;
   	}
      }
    }
  }

  imatrix listz(Ni);
  int p;
  for(e=1; e<=Ne2; e++){
    for(v=1; v<=Np2; v++){
      p = nodeId(e,v);
      if (p>0)
   	listz(p)++;
    }
  }

  offset.resize(Ni+1);
  offset(1) = 0;
  for(int i=1; i<=Ni; i++)
    offset(i+1) = offset(i) + listz(i);

  list.resize(offset(Ni+1));
  for(int e=1; e<=Ne2; e++){
    for(int v=1; v<=Np2; v++){
      p = nodeId(e,v);
      if (p>0){
   	list(offset(p+1)-listz(p)+1) = (e-1)*Np2+v;
   	listz(p) = listz(p)-1;
      }
    }
  }
}


void sem::assembleMatrices(){

  setupAide gllData("gll/gllData.dat");
  char datanp[BUFSIZ];
  sprintf(datanp, "NODES %2.2d", N);
  gllData.getArgs(xLoc, datanp);
  sprintf(datanp, "WEIGHT %2.2d", N);
  gllData.getArgs(w, datanp);

  // matrix MLoc
  MLoc.resize(Np2);
  int t = 0;
  for(int i=1; i<=Np1; i++){
    for(int j=1; j<=Np1; j++){
      t++;
      MLoc(t) = w(i)*w(j);
    }
  }

  // global matrix M
  M.resize(Ni);
  for(int e=1; e<=Ne2; e++){
    for(int v=1; v<=Np2; v++){
      int p = nodeId(e,v);
      if (p>0)
   	M(p) = M(p) + MLoc(v)*J;
    }
  }
  MJ = J*MLoc;

  // stiffness matrix S
  S.resize(Np2,Np2);
  S = 0;
  // d/dx
  int t1;
  t1 = 0;
  for(int i1=1; i1<=Np1; i1++){
    for (int j = 1; j<=Np1; j++){
      t1++;
      for(int i2=1; i2<=Np1; i2++){
   	int t2 = (i2-1)*Np1 + j;
   	datafloat a = 0;
   	for(int q=1; q<=Np1; q++){
   	  a += partialD(i1,xLoc(q))*partialD(i2,xLoc(q))*w(q);
   	}
   	a *= w(j);
   	S(t1,t2) += a;
      }
    }
  }

  // d/dy
  t1 = 0;
  for(int i=1; i<=Np1; i++){
    for (int j1 = 1; j1<=Np1; j1++){
      t1++;
      for(int j2=1; j2<=Np1; j2++){
   	int t2 = (i-1)*Np1 + j2;
   	datafloat a = 0;
   	for(int q=1; q<=Np1; q++){
   	  a += partialD(j1,xLoc(q))*partialD(j2,xLoc(q))*w(q);
   	}
   	a *= w(i);
   	S(t1,t2) += a;
      }
    }
  }

  xi.resize(Ni);
  yi.resize(Ni);

  for(int i=1; i<=Ni; i++){
    int p = list(offset(i)+1);
    int e = (p+Np2-1)/Np2;
    int ei = (e+K-1)/K;
    int ej = e%K;
    if (ej==0)
      ej = K;
    int v = p-(e-1)*Np2;
    int vi = (v+Np1-1)/Np1;
    int vj = v%Np1;
    if (vj==0)
      vj = Np1;
    xi(i) = (ei-1)*N + vi;
    yi(i) = (ej-1)*N + vj;
  }

  r.resize(K*N);
  for(int e=1; e<=K; e++){
    for(int v=1; v<=N; v++){
      int k = (e-1)*N + v;
      r(k) = h/2*(xLoc(v)+1) + (e-1)*h - 1;
    }
  }

}

datafloat sem::partialD(int n, datafloat x){

  datafloat result = 0;
  datafloat c = 1.0; // coefficient
  for(int i=1; i<=Np1; i++)
    if (i!=n)
      c *= xLoc(n)-xLoc(i);

  datafloat s;
  for(int i=1; i<=Np1; i++){
    if (i!=n){
      s = 1;
      for(int j=1; j<=Np1; j++)
   	if ((j!=n) && (j!=i))
   	  s *= x-xLoc(j);
      result += s;
    }
  }
  result /= c;

  return result;
}


void sem::setupRHS(){

  aAul = 0;
  for(int e=1; e<=Ne2; e++){
    for(int v=1; v<=Np2; v++){
      int id = (e-1)*Np2 + v;
      int p = nodeId(id);
      if (p > 0){
	datafloat x = r(xi(p));
	datafloat y = r(yi(p));
	aAul(id) = source(x,y)*MJ(v);	
      }
    }
  }
  gather(aAul, b);
}

void sem::gather(fmatrix &uAul, fmatrix &u){
#pragma omp parallel for			\
  default(shared) num_threads(NThread)  
  for(int i=1; i<=Ni; i++){
    int start = offset(i) + 1;
    int end = offset(i+1);
    datafloat sum = 0;
    for(int j=start; j<=end; j++){
      sum += uAul(list(j));
    }
    u(i) = sum;
  }
}

datafloat sem::source(datafloat x, datafloat y){
  // return 2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
  return -2*exp((x*x-1)*(y*y-1))*(-2+3*y*y+2*x*x*x*x*y*y+x*x*(3-8*y*y+2*y*y*y*y));
}

void sem::matVec(fmatrix &u, fmatrix &Au){
  aAul = 0;
#pragma omp parallel for			\
  default(shared) num_threads(NThread)
  for(int e=1; e<=Ne2; e++){
    for(int v1=1; v1<=Np2; v1++){
      datafloat sum = 0;
      for(int v2=1; v2<=Np2; v2++){
	int p = nodeId(e,v2);
	if (p > 0){
	  sum += S(v1,v2)*u(p);
	}
      }
      aAul((e-1)*Np2+v1) = sum;
    }
  }
  gather(aAul, Au);
}

void sem::setSchwarzId(int level){

#pragma omp parallel for			\
  default(shared) num_threads(NThread)  
  for(int ei=1; ei<=Ne1; ei++){
    for(int ej=1; ej<=Ne1; ej++){

      int e = (ei-1)*Ne1 + ej;

      // left bottom
      if ((ei>1)&&(ej>1)){
	int en = e - Ne1 - 1;
	int vn = Np2 - Np1 - 1;
	int p = nodeId(en,vn);
	if (p > 0)
	  id[e].push_back(p);
      }

      // left
      if (ei>1){
	int en = e - Ne1;
	for(int vn=Np2-2*Np1+1; vn<=Np2-Np1; vn++){
	  int p = nodeId(en,vn);
	  if (p > 0)
	    id[e].push_back(p);
	}
      }

      // left top
      if ((ei>1)&&(ej<Ne1)){
	int en = e - Ne1 + 1;
	int vn = Np2 - 2*Np1 + 2;
	int p = nodeId(en,vn);
	if (p > 0)
	  id[e].push_back(p);
      }

      // bottom
      if (ej>1){
	int en = e - 1;
	for(int vn=Np1-1; vn<=Np2-1; vn+=Np1){
	  int p = nodeId(en,vn);
	  if (p > 0)
	    id[e].push_back(p);
	}	
      }

      // self
      for(int v=1; v<=Np2; v++){
	int p = nodeId(e,v);
	if (p > 0)
	  id[e].push_back(p);
      }


      // top
      if (ej<Ne1){
	int en = e + 1;
	for(int vn=2; vn<=Np2-Np1+2; vn+=Np1){
	  int p = nodeId(en,vn);
	  if (p > 0)
	    id[e].push_back(p);
	}
      }

      // right bottom
      if ((ei<Ne1)&&(ej>1)){
	int en = e + Ne1 - 1;
	int vn = 2*Np1 - 1;
	int p = nodeId(en,vn);
	if (p > 0)
	  id[e].push_back(p);
      }

      // right
      if (ei<Ne1){
	int en = e + Ne1;
	for(int vn=Np1+1; vn<=2*Np1; vn++){
	  int p = nodeId(en,vn);
	  if (p > 0)
	    id[e].push_back(p);
	}	
      }

      // right top
      if ((ei<Ne1)&&(ej<Ne1)){
	int en = e + Ne1 + 1;
	int vn = Np1 + 2;
	int p = nodeId(en,vn);
	if (p > 0)
	  id[e].push_back(p);
      }      

      //sort(id[e].begin(),id[e].end());
      
    }
  }
  
  if (level==2){
    // Coarse Grid
    int inc = Ne1/2;
    for(int e=1; e<=Ne2; e+=inc){
      int v = (Np2+1)/2;
      int p = nodeId(e,v);
      if (p > 0)
	id[0].push_back(p);       
    }
  }
 
  /*
  for(int i=0; i<=0; i++){
    for(int j=0; j<id[i].size(); j++){
      cout << id[i][j] << " ";
    }
    cout << endl;
  }
  */

}

void sem::setSchwarzMatrix(int level){

  int start;
  if (level==1)
    start = 1;
  else if (level==2)
    start = 0;

#pragma omp parallel for			\
  default(shared) num_threads(NThread)  
  for(int e=start; e<=Ne2; e++){
    int sz = id[e].size();
    A[e].resize(sz,sz);
    uLoc[e].resize(sz);
    A[e] = 0;
    for(int a1=1; a1<=sz; a1++){
      int p1 = id[e][a1-1];
      int start1 = offset(p1) + 1;
      int end1 = offset(p1 + 1);
      for(int i1=start1; i1<=end1; i1++){
	int l1 = list(i1) - 1; // 0-indexed
	int e1 = l1/Np2 + 1;
	int v1 = l1%Np2 + 1;

	for(int a2=1; a2<=sz; a2++){	
	  int p2 = id[e][a2-1];
	  int start2 = offset(p2) + 1;
	  int end2 = offset(p2 + 1);
	  for(int i2=start2; i2<=end2; i2++){
	    int l2 = list(i2) - 1; // 0-indexed
	    int e2 = l2/Np2 + 1;
	    int v2 = l2%Np2 + 1;
	  
	    if (e1==e2)
	      A[e](a1,a2) += S(v1,v2);

	  }
	}
      }
    }
    
    A[e].lu(IPIV[e]);
  }  
}

void sem::additiveSchwarz(fmatrix &u, fmatrix &Mu, int &start, int &end){


  Mu = 0;  
  
#pragma omp parallel for			\
  default(shared) num_threads(NThread)    
  for(int e=start; e<=end; e++){
    int sz = id[e].size();
    
    // gather uLoc from u
    for(int v=1; v<=sz; v++){
      int p = id[e][v-1];
      uLoc[e](v) = u(p);
    }

    //uLoc[e] = A[e] | uLoc[e];
    A[e].solve(IPIV[e], uLoc[e]);

    // scatter uLoc back to u
    for(int v=1; v<=sz; v++){
      int p = id[e][v-1];
#pragma omp atomic
      Mu(p) += uLoc[e](v);
    }      
  }
  
}

void sem::initSchwarz(int level){
  id = (vector<int> *)calloc(Ne2+1,sizeof(vector<int>));
  setSchwarzId(level);
  A = (fmatrix *)calloc(Ne2+1,sizeof(fmatrix));
  IPIV = (imatrix *)calloc(Ne2+1,sizeof(imatrix));
  uLoc = (fmatrix *)calloc(Ne2+1,sizeof(fmatrix));
  setSchwarzMatrix(level);
}

void sem::pcg(int level){

  int start, end;
  if(level == 1){
    start = 1;
    end = Ne2;
  }
  else if(level == 2){
    start = 0;
    end = Ne2;
  }

  cout << "Preconditioned Conjugate Gradient (" << level << "):" << endl;

  tStart = omp_get_wtime();
  
  initSchwarz(level);

  uSoln = 0;

  fmatrix Ap(Ni), p(Ni), r(Ni), z(Ni);
  datafloat resOld, resNew, alpha, pAp;
  datafloat res;
  r = b;

  //z = globalS | r;
  additiveSchwarz(r,z,start,end);

  p = z;

  resOld = 0;
  for(int i=1; i<=Ni; i++)
    resOld += r(i)*z(i);

  for(iter =1; iter<=iterMax; iter++){
    matVec(p,Ap);
    
    pAp = 0;
    for(int i=1; i<=Ni; i++)
      pAp += p(i)*Ap(i);
    
    alpha = resOld/pAp;
    
    uSoln = uSoln + alpha*p;
    
    r = r - alpha*Ap;

    res = 0;
    for(int i=1; i<=Ni; i++)
      res += r(i)*r(i);
    if (res < tol*tol){
      cout << "Norm of Residue = " << sqrt(res) << endl;
      break;
    }
    //z = globalS | r;
    additiveSchwarz(r,z,start,end);
    
    resNew = 0;
    for(int i=1; i<=Ni; i++)
      resNew += r(i)*z(i);
    
    p = z + (resNew/resOld)*p;
      
    resOld = resNew;
  }
  
  tEnd = omp_get_wtime();
  elapsed = ((double)(tEnd-tStart));
  
}

void sem::cg(){

  cout << "Conjugate Gradient: " << endl;

  tStart = omp_get_wtime();
  
  uSoln = 0;

  fmatrix Ap(Ni), p(Ni), r(Ni);
  datafloat resOld, resNew, alpha, pAp;
  r = b;
  p = r;

  resOld = 0;
  for(int i=1; i<=Ni; i++)
    resOld += r(i)*r(i);

  for(iter=1; iter<=iterMax; iter++){
    matVec(p,Ap);
    
    pAp = 0;
    for(int i=1; i<=Ni; i++)
      pAp += p(i)*Ap(i);

    alpha = resOld/pAp;

    uSoln = uSoln + alpha*p;
    
    r = r - alpha*Ap;
    
    resNew = 0;
    for(int i=1; i<=Ni; i++)
      resNew += r(i)*r(i);

    if (resNew < tol*tol){
      cout << "Norm of Residue = " << sqrt(resNew) << endl;
      break;
    }
    
    p = r + (resNew/resOld)*p;

    resOld = resNew;
    
  }

  tEnd = omp_get_wtime();
  elapsed = ((double)(tEnd-tStart));
  
}


void sem::postProcessing(){  
  uExact.resize(Ni);
  for(int i=1; i<=Ni; i++){
    datafloat x = r(xi(i));
    datafloat y = r(yi(i));
    //uExact(i) = sin(M_PI*x)*sin(M_PI*y);
    uExact(i) = exp((x*x-1)*(y*y-1))-1;
  }

  datafloat l2err = 0;
  for(int i=1; i<=Ni; i++){
    l2err = (uExact(i)-uSoln(i))*(uExact(i)-uSoln(i))*M(i);
  }

  cout << "Iteration #     = " << iter << endl;
  cout << "Elapsed Time    = " << elapsed << endl;
  cout << "L2 error        = " << l2err << endl;
}

void sem::matrixGen(){

  globalS.resize(Ni,Ni);
  fmatrix e(Ni), Ae(Ni);
  e = 0;
  e(1) = 1;
  matVec(e,Ae);
  for(int i=1; i<=Ni; i++)
    globalS(i,1) = Ae(i);

  for(int r=2; r<=Ni; r++){
    e(r-1) = 0;
    e(r) = 1;
    matVec(e,Ae);  
    for(int i=1; i<=Ni; i++)
      globalS(i,r) = Ae(i);
  }

}
