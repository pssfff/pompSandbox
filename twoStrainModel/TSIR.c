/* pomp model file: TSIR */

#include <pomp.h>
#include <R_ext/Rdynload.h>

#define N	(__p[__parindex[0]])
#define mu	(__p[__parindex[1]])
#define rho1	(__p[__parindex[2]])
#define rho2	(__p[__parindex[3]])
#define beta1	(__p[__parindex[4]])
#define alpha1	(__p[__parindex[5]])
#define alpha2	(__p[__parindex[6]])
#define iota	(__p[__parindex[7]])
#define lambda	(__p[__parindex[8]])
#define S1_0	(__p[__parindex[9]])
#define I1_0	(__p[__parindex[10]])
#define C1_0	(__p[__parindex[11]])
#define S2_0	(__p[__parindex[12]])
#define I2_0	(__p[__parindex[13]])
#define C2_0	(__p[__parindex[14]])
#define S1	(__x[__stateindex[0]])
#define I1	(__x[__stateindex[1]])
#define C1	(__x[__stateindex[2]])
#define S2	(__x[__stateindex[3]])
#define I2	(__x[__stateindex[4]])
#define C2	(__x[__stateindex[5]])
#define cases1	(__y[__obsindex[0]])
#define cases2	(__y[__obsindex[1]])
#define DS1	(__f[__stateindex[0]])
#define DI1	(__f[__stateindex[1]])
#define DC1	(__f[__stateindex[2]])
#define DS2	(__f[__stateindex[3]])
#define DI2	(__f[__stateindex[4]])
#define DC2	(__f[__stateindex[5]])
#define TN	(__pt[__parindex[0]])
#define Tmu	(__pt[__parindex[1]])
#define Trho1	(__pt[__parindex[2]])
#define Trho2	(__pt[__parindex[3]])
#define Tbeta1	(__pt[__parindex[4]])
#define Talpha1	(__pt[__parindex[5]])
#define Talpha2	(__pt[__parindex[6]])
#define Tiota	(__pt[__parindex[7]])
#define Tlambda	(__pt[__parindex[8]])
#define TS1_0	(__pt[__parindex[9]])
#define TI1_0	(__pt[__parindex[10]])
#define TC1_0	(__pt[__parindex[11]])
#define TS2_0	(__pt[__parindex[12]])
#define TI2_0	(__pt[__parindex[13]])
#define TC2_0	(__pt[__parindex[14]])
#define lik	(__lik[0])

void TSIR_par_trans (double *__pt, double *__p, int *__parindex)
{

        TN = exp(N);
        Tmu = exp(mu);
        Tbeta1 = exp(beta1);
        Talpha1 = exp(alpha1);
        Talpha2 = exp(alpha2);
        Tiota = exp(iota);
        Tlambda = exp(lambda);
        TS1_0 = exp(S1_0);
        TI1_0 = exp(I1_0);
        TC1_0 = exp(C1_0);
        TS2_0 = exp(S2_0);
        TI2_0 = exp(I2_0);
        TC2_0 = exp(C2_0);
        Trho1 = expit(rho1);
        Trho2 = expit(rho2);
 
}


void TSIR_par_untrans (double *__pt, double *__p, int *__parindex)
{

        TN = log(N);
        Tmu = log(mu);
        Tbeta1 = log(beta1);
        Talpha1 = log(alpha1);
        Talpha2 = log(alpha2);
        Tiota = log(iota);
        Tlambda = log(lambda);
        TS1_0 = log(S1_0);
        TI1_0 = log(I1_0);
        TC1_0 = log(C1_0);
        TS2_0 = log(S2_0);
        TI2_0 = log(I2_0);
        TC2_0 = log(C2_0);
        Trho1 = logit(rho1);
        Trho2 = logit(rho2);
 
}


void TSIR_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 
        cases1 = rbinom(I1, rho1);
        cases2 = rbinom(I2, rho2);
 
}


void TSIR_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 
        double lik1 = dbinom(cases1, I1, rho1, give_log);
        double lik2 = dbinom(cases2, I2, rho2, give_log);
        lik = lik1*lik2;
 
}


void TSIR_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{

        // new infections
        I1 = rpois(beta1*pow(I1/N+iota, alpha1)*pow(S1, alpha2));
        I2 = rpois(beta1*pow(I2/N+iota, alpha1)*pow(S2, alpha2));
        // Poisson approximation to exponential losses
        int C1_loss = rpois(lambda*C1);
        int C2_loss = rpois(lambda*C2);
        // udpate counts
        S1 += N*mu - I1 - I2 + C1_loss;
        S2 += N*mu - I2 - I1 + C2_loss;
        C1 += I2 - C1_loss;
        C2 += I1 - C2_loss;
 
}


void TSIR_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{

return;
 
}

#undef N
#undef mu
#undef rho1
#undef rho2
#undef beta1
#undef alpha1
#undef alpha2
#undef iota
#undef lambda
#undef S1_0
#undef I1_0
#undef C1_0
#undef S2_0
#undef I2_0
#undef C2_0
#undef S1
#undef I1
#undef C1
#undef S2
#undef I2
#undef C2
#undef cases1
#undef cases2
#undef DS1
#undef DI1
#undef DC1
#undef DS2
#undef DI2
#undef DC2
#undef TN
#undef Tmu
#undef Trho1
#undef Trho2
#undef Tbeta1
#undef Talpha1
#undef Talpha2
#undef Tiota
#undef Tlambda
#undef TS1_0
#undef TI1_0
#undef TC1_0
#undef TS2_0
#undef TI2_0
#undef TC2_0
