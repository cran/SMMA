//// [[Rcpp::depends(RcppArmadillo)]]
//#define TIMING
#include <RcppArmadillo.h>
#include "auxfunc.h"
//#include "/Users/adamlund/Documents/KU/Phd/Project/Computer/Vincent/timer/simple_timer.h"

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////algorithm ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
Rcpp::List pga(arma::mat Phi1, arma::mat Phi2, arma::mat Phi3,
               Rcpp::NumericVector resp,
               std::string penalty,
               double zeta,
               double c,
               arma::vec lambda,
               int nlambda,
               int makelamb,
               double lambdaminratio,
               arma::mat penaltyfactor,
               double reltol,
               int maxiter,
               int steps,
               int btmax,
               double Delta0,
               double nu,
               int alg,
               int ll){

Rcpp::List output;
Rcpp::NumericVector vecY(resp);//??why thevecY??
Rcpp::IntegerVector YDim = vecY.attr("dim");
const arma::cube Z(vecY.begin(), YDim[0], YDim[1], YDim[2], false);//whybrackets?????????

//declare some global variables
int ascent, ascentmax,
    bt, btenter = 0, btiter = 0,
    endmodelno = nlambda,
    n1 = Phi1.n_rows, n2 = Phi2.n_rows, n3 = Phi3.n_rows, Nog = Z.n_slices, ng = n1 * n2 * n3,// n = Nog * ng,
    p1 = Phi1.n_cols, p2 = Phi2.n_cols, p3 = Phi3.n_cols, p = p1 * p2 * p3,
    Stopconv = 0, Stopmaxiter = 0, Stopbt = 0;

double alphamax, ascad = 3.7,
       delta,
       eta = 0.9,
       Lmax, lossBeta, lossBetaprev, lossProp, lossX,
       penBeta, penProp, penX,
       Qdelta,
       relobj,
       tk, tkp1,
       val;

arma::vec df(nlambda),
          eig1, eig2, eig3,
          Iter(nlambda), Loss(maxiter + 1), Pen(maxiter + 1),
          obj(maxiter + 1),
          Stops(3),
          eevX;

arma::mat absBeta(p1, p2 * p3),
          Beta(p1, p2 * p3), Betaprev(p1, p2 * p3), Betas(p, nlambda), BT(nlambda, maxiter + 1),
          dpen(p1, p2 * p3),
          Gamma(p1, p2 * p3), GradlossX(p1, p2 * p3), GradlossX2(p1, p2 * p3),
          Obj(maxiter + 1, nlambda),
          Phi1tPhi1, Phi2tPhi2, Phi3tPhi3, PhitPhiBeta, PhitPhiX, pospart(p1, p2 * p3), Prop(p1, p2 * p3), PhiBeta(n1, n2 * n3), PhiProp(n1, n2 * n3), PhiX(n1, n2 * n3),
          wGamma(p1, p2 * p3),
          X(p1, p2 * p3);

arma::cube PhitZ(p1, p2 * p3, Nog), W(Nog, maxiter + 1, nlambda);

////fill variables
ascentmax = 4;

obj.fill(NA_REAL);
Obj.fill(0);
Betas.fill(42);
Iter.fill(0);

BT.fill(-1);

W.fill(42);
Loss.fill(0);
Pen.fill(0);

////precompute
Phi1tPhi1 = Phi1.t() * Phi1;
Phi2tPhi2 = Phi2.t() * Phi2;
Phi3tPhi3 = Phi3.t() * Phi3;
eig1 = arma::eig_sym(Phi1tPhi1);
eig2 = arma::eig_sym(Phi2tPhi2);
eig3 = arma::eig_sym(Phi3tPhi3);
alphamax = as_scalar(max(kron(eig1, kron(eig2 , eig3))));

////precompute
for(int j = 0; j < Nog; j++){

PhitZ.slice(j) = RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), Z.slice(j), n2, n3), n3, p1), p1, p2);

}

////proximal step size
Lmax = 2 * alphamax * (1 + Delta0 * zeta) / ng; //upper bound on Lipschitz?!?!?!?!?!?!?!?!?!?!?!?!?

//initial step size
if(nu > 0){delta = 1.9 / (nu * Lmax);}else{delta = Delta0;}

////initialize
Betaprev.fill(0);
Beta = Betaprev;
PhiBeta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);
PhitPhiBeta = RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, Beta, p2, p3), p3, p1), p1, p2);

lossBeta = softmaxloss(-eev(PhiBeta, Z, ng), zeta, ll);

////make lambda sequence
if(makelamb == 1){

arma::mat Ze = zeros<mat>(n1, n2 * n3);
arma::mat absgradzeroall(p1, p2 * p2);

absgradzeroall = abs(gradloss(PhitZ, PhitPhiBeta, -eev(PhiBeta, Z, ng), ng, zeta, ll));


arma::mat absgradzeropencoef = absgradzeroall % (penaltyfactor > 0);
arma::mat penaltyfactorpencoef = (penaltyfactor == 0) * 1 + penaltyfactor;
double lambdamax = as_scalar(max(max(absgradzeropencoef / penaltyfactorpencoef)));
double m = log(lambdaminratio);
double M = 0;
double difflamb = abs(M - m) / (nlambda - 1);
double l = 0;

for(int i = 0; i < nlambda ; i++){

lambda(i) = lambdamax * exp(l);
l = l - difflamb;

}

}else{std::sort(lambda.begin(), lambda.end(), std::greater<int>());}

///////////start lambda loop
for (int j = 0; j < nlambda; j++){

Gamma = penaltyfactor * lambda(j);

ascent = 0;

//start MSA loop
for (int s = 0; s < steps; s++){

if(s == 0){

if(penalty != "lasso"){wGamma = Gamma / lambda(j);}else{wGamma = Gamma;}

}else{

if(penalty == "scad"){

absBeta =abs(Beta);
pospart = ((ascad * Gamma - absBeta) + (ascad * Gamma - absBeta)) / 2;
dpen = sign(Beta) % Gamma % ((absBeta <= Gamma) + pospart / (ascad - 1) % (absBeta > Gamma));
wGamma = abs(dpen) % Gamma / lambda(j) % (Beta != 0) + lambda(j) * (Beta == 0);

}

}

/////start proximal loop
if(alg == 1){/////////////////NPG algorithm from chen2016

for (int k = 0; k < maxiter; k++){

if(k == 0){

Betaprev = Beta;
obj(k) = lossBeta + l1penalty(wGamma, Beta);
Obj(k, j) = obj(k);

// W.slice(j).col(k) = wBeta;

BT(j, k) = 1; //force initial backtracking (if deltamin < delta)

}else{//if not the first iteration

X = Beta + (k - 2) / (k + 1) * (Beta - Betaprev);

PhiX = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, X, p2, p3), p3, n1), n1, n2);
eevX = -eev(PhiX, Z, ng);
//wX = exp(-zeta * -eev(PhiX, Z, ng));

PhitPhiX = RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, X, p2, p3), p3, p1), p1, p2);
GradlossX = gradloss(PhitZ, PhitPhiX, eevX, ng, zeta, ll);

//if(ll == 1){GradlossX = GradlossX / accu(wX);}

////check if proximal backtracking occurred last iteration
if(BT(j, k - 1) > 0){bt = 1;}else{bt = 0;}

lossX = softmaxloss(eevX, zeta, ll);
//lossX = accu(wX);

////proximal backtracking from chen2016
BT(j, k) = 0;

while (BT(j, k) < btmax){

Prop = prox_l1(X - delta * GradlossX, delta * wGamma);
PhiProp = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Prop, p2, p3), p3, n1), n1, n2);
lossProp = softmaxloss(-eev(PhiProp, Z, ng), zeta, ll);

val = as_scalar(max(obj(span(0, k - 1))) - c / (2 * delta) * sum_square(Prop - Betaprev));
penProp = l1penalty(wGamma, Prop);

if (lossProp + penProp <= val + 0.0000001){ //need to add a little due to numerical issues??

break;

}else{

delta = eta * delta;
BT(j, k) = BT(j, k) + 1;

//if(delta < deltamin){delta = deltamin;}

}

}//end line search

////check if maximum number of proximal backtraking step is reached
if(BT(j, k) == btmax){Stopbt = 1;}

Betaprev = Beta;
Beta = Prop;
// W.slice(j).col(k) = wBeta;
lossBeta = lossProp;
penBeta = penProp;
obj(k) = lossBeta + penBeta;
Loss(k) = lossBeta;
Pen(k) = penBeta;
Obj(k, j) = obj(k);
Iter(j) = k;

////proximal convergence check
relobj = abs(obj(k) - obj(k - 1)) / (reltol + abs(obj(k - 1)));

if(k < maxiter && relobj < reltol){

df(j) = p - accu((Beta == 0));
Betas.col(j) = vectorise(Beta);
obj.fill(NA_REAL);
Stopconv = 1;
break;

}else if(k == maxiter){

df(j) = p - accu((Beta == 0));
Betas.col(j) = vectorise(Beta);
obj.fill(NA_REAL);
Stopmaxiter = 1;
break;

}

}

////break proximal loop if maximum number of proximal backtraking step is reached
if(Stopbt == 1){

Betas.col(j) = vectorise(Beta);
break;

}

}//end proximal loop

}else{ //MFISTA algortihm plus a backtracking  from parihkandboyd///////////////////////////////////////////////////

for (int k = 0; k < maxiter; k++){

if(k == 0){

Betaprev = Beta;
X = Beta;
obj(k) = lossBeta + l1penalty(wGamma, Beta);
Obj(k, j) = obj(k);
// W.slice(j).col(k) = wBeta;
BT(j, k) = 1; //force initial backtracking
tk = 1;

}else{//if not the first iteration

PhiX = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, X, p2, p3), p3, n1), n1, n2);
//wX = exp(-zeta * -eev(PhiX, Z, ng));
eevX = -eev(PhiX, Z, ng);
PhitPhiX = RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, X, p2, p3), p3, p1), p1, p2);
GradlossX = gradloss(PhitZ, PhitPhiX, eevX, ng, zeta, ll);

////check if backtracking occurred last iteration
if(BT(j, k - 1) > 0){bt = 1;}else{bt = 0;}

if(bt == 1 || nu == 0){// if backtrack

  lossX = softmaxloss(eevX, zeta, ll);
//lossX = accu(wX);
penX = l1penalty(wGamma, X);
BT(j, k) = 0;

while (BT(j, k) < btmax){//line search

Prop = prox_l1(X - delta * GradlossX, delta * wGamma);
PhiProp = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Prop, p2, p3), p3, n1), n1, n2);
lossProp = softmaxloss(-eev(PhiProp, Z, ng), zeta, ll);

Qdelta = as_scalar(lossX + accu(GradlossX % (Prop - X)) + 1 / (2 * delta) * sum_square(Prop - X) + 0 * penX);
penProp = l1penalty(wGamma, Prop);

if(lossProp + 0 * penProp <= Qdelta + 0.0000001){ //need to add a little ??

break;

}else{

delta = eta * delta;
BT(j, k) = BT(j, k) + 1;

}

}//end line search

////check if maximum number of proximal backtraking step is reached
if(BT(j, k) == btmax){Stopbt = 1;}

}else{ //if no backtracking

Prop = prox_l1(X - delta * GradlossX, delta * wGamma);
PhiProp = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Prop, p2, p3), p3, n1), n1, n2);
lossProp = softmaxloss(-eev(PhiProp, Z, ng), zeta, ll);

}

tkp1 = (1 + sqrt(1 + 4 * pow(tk, 2))) / 2;

if(lossProp + penProp <= obj(k - 1)){

Betaprev = Beta;
lossBetaprev = lossBeta;
Beta = Prop;
lossBeta = lossProp;
penBeta = penProp;
obj(k) = lossBeta + penBeta;
// W.slice(j).col(k) = wBeta;
Loss(k) = lossBeta;
Pen(k) = penBeta;

}else{

Betaprev = Beta;
lossBetaprev = lossBeta;
obj(k) = lossBeta + penBeta;
// W.slice(j).col(k) = wBeta;
Loss(k) = lossBeta;
Pen(k) = penBeta;

}

X = Beta + tk / tkp1 * (Prop - Beta) + (tk - 1) / tkp1 * (Beta - Betaprev); //yk+1

Obj(k, j) = obj(k);
Iter(j) = k;

////proximal convergence check
relobj = abs(obj(k) - obj(k - 1)) / (reltol + abs(obj(k - 1)));
 // diffbet = accu(abs(Beta - Betaprev));

if(k < maxiter && relobj < reltol //&& diffbet != 0 //|| k == 0
){

df(j) = p - accu((Beta == 0));
Betas.col(j) = vectorise(Beta);
obj.fill(NA_REAL);
Stopconv = 1;
break;

}else if(k == maxiter){

df(j) = p - accu((Beta == 0));
Betas.col(j) = vectorise(Beta);
obj.fill(NA_REAL);
Stopmaxiter = 1;
break;

}

}

////break proximal loop if maximum number of proximal backtraking step is reached
if(Stopbt == 1){

Betas.col(j) = vectorise(Beta);
break;

}

}//end proximal loop

df(j) = p - accu((Beta == 0));

}

//Stop program if maximum number of backtracking steps or maxiter is reached
if(Stopbt == 1){

endmodelno = j;
break;

}

}//end MSA loop

}//end lambda loop

Stops(0) = Stopconv;
Stops(1) = Stopmaxiter;
Stops(2) = Stopbt;
btenter = accu((BT > -1));
btiter = accu((BT > 0) % BT);

output = Rcpp::List::create(Rcpp::Named("Beta") = Betas,
                            Rcpp::Named("df") = df,
                            Rcpp::Named("btenter") = btenter,
                            Rcpp::Named("btiter") = btiter,
                            Rcpp::Named("Obj") = Obj,
                            Rcpp::Named("Iter") = Iter,
                            Rcpp::Named("endmodelno") = endmodelno,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("BT") = BT,
                            Rcpp::Named("Stops") = Stops);

return output;

}
