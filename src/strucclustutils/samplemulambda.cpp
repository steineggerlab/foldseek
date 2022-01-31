#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"
#include "StructureUtil.h"

#include <random>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

/************************************************************
 * This Code is copied from HMMer3 and modified to suit our needs
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *
 ************************************************************/

/* histogram.c
 * SRE, Sat Jan 20 16:16:17 1996
 *
 * Accumulation, printing, and fitting of score histograms
 * from database searches.
 *
 * CVS $Id: histogram.c,v 1.18 2003/04/14 16:00:16 eddy Exp $
 ************************************************************
 * Basic API:
 *
 * struct histogram_s *h;
 *
 * h = AllocHistogram(min_hint, max_hint, lumpsize);
 *
 * while (getting scores x) AddToHistogram(h, x);
 *
 * ExtremeValueFitHistogram(h, high_hint);
 * PrintASCIIHistogram(fp, h);
 * FreeHistogram(h);
 */

#ifndef MIN
#define MIN(a,b)         (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)         (((a)>(b))?(a):(b))
#endif

/* Function: AllocHistogram()
 *
 * Purpose:  Allocate and return a histogram structure.
 *           min and max are your best guess. They need
 *           not be absolutely correct; the histogram
 *           will expand dynamically to accomodate scores
 *           that exceed these suggested bounds. The amount
 *           that the histogram grows by is set by "lumpsize".
 *
 * Args:     min:      minimum score (integer)
 *           max:      maximum score (integer)
 *           lumpsize: when reallocating histogram, pad the reallocation
 *                     by this much (saves excessive reallocation)
 */

struct histogram_s {
    int *histogram;		/* counts of hits                     */
    int  min;			/* elem 0 of histogram == min         */
    int  max;                     /* last elem of histogram == max      */
    int  highscore;		/* highest active elem has this score */
    int  lowscore;		/* lowest active elem has this score  */
    int  lumpsize;		/* when resizing, overalloc by this   */
    int  total;			/* total # of hits counted            */

    float *expect;		/* expected counts of hits            */
    int    fit_type;		/* flag indicating distribution type  */
    float  param[3];		/* parameters used for fits           */
    float  chisq;			/* chi-squared val for goodness of fit*/
    float  chip;			/* P value for chisquared             */
};
#define HISTFIT_NONE     0	/* no fit done yet               */
#define HISTFIT_EVD      1	/* fit type = extreme value dist */
#define HISTFIT_GAUSSIAN 2	/* fit type = Gaussian           */
#define EVD_MU		 0	/* EVD fit parameter mu          */
#define EVD_LAMBDA       1	/* EVD fit parameter lambda      */
#define EVD_WONKA        2      /* EVD fit fudge factor          */
#define GAUSS_MEAN       0	/* Gaussian parameter mean       */
#define GAUSS_SD         1	/* Gaussian parameter std. dev.  */


/* Function: ExtremeValueP()
 *
 * Purpose:  Calculate P(S>x) according to an extreme
 *           value distribution, given x and the parameters
 *           of the distribution (characteristic
 *           value mu, decay constant lambda).
 *
 *           This function is exquisitely prone to
 *           floating point exceptions if it isn't coded
 *           carefully.
 *
 * Args:     x      = score
 *           mu     = characteristic value of extreme value distribution
 *           lambda = decay constant of extreme value distribution
 *
 * Return:   P(S>x)
 */
double
ExtremeValueP(float x, float mu, float lambda)
{
    double y;
    /* avoid exceptions near P=1.0 */
    /* typical 32-bit sys: if () < -3.6, return 1.0 */
    if ((lambda * (x - mu)) <= -1. * log(-1. * log(DBL_EPSILON))) return 1.0;
    /* avoid underflow fp exceptions near P=0.0*/
    if ((lambda * (x - mu)) >= 2.3 * (double) DBL_MAX_10_EXP)     return 0.0;
    /* a roundoff issue arises; use 1 - e^-x --> x for small x */
    y = exp(-1. * lambda * (x - mu));
    if       (y < 1e-7) return y;
    else     return (1.0 - exp(-1. * y));
}




static int sre_randseed = 42;	/* default seed for sre_random()   */

/* Function: sre_random()
 *
 * Purpose:  Return a uniform deviate x, 0.0 <= x < 1.0.
 *
 *           sre_randseed is a static variable, set
 *           by sre_srandom(). When it is non-zero,
 *           we re-seed.
 *
 *           Implements L'Ecuyer's algorithm for combining output
 *           of two linear congruential generators, plus a Bays-Durham
 *           shuffle. This is essentially ran2() from Numerical Recipes,
 *           sans their nonhelpful Rand/McNally-esque code obfuscation.
 *
 *           Overflow errors are avoided by Schrage's algorithm:
 *               az % m = a(z%q) - r(z/q) (+m if <0)
 *           where q=m/a, r=m%a
 *
 *           Requires that long int's have at least 32 bits.
 *           This function uses statics and is NOT THREADSAFE.
 *
 * Reference: Press et al. Numerical Recipes in C, 1992.
 *
 * Reliable and portable, but slow. Benchmarks on wrasse,
 * using Linux gcc and Linux glibc rand() (see randspeed, in Testsuite):
 *     sre_random():    0.5 usec/call
 *     rand():          0.2 usec/call
 */
double
sre_random(void)
{
    static long  rnd1;		/* random number from LCG1 */
    static long  rnd2;            /* random number from LCG2 */
    static long  rnd;             /* random number we return */
    static long  tbl[64];		/* table for Bays/Durham shuffle */
    long x,y;
    int i;

    /* Magic numbers a1,m1, a2,m2 from L'Ecuyer, for 2 LCGs.
     * q,r derive from them (q=m/a, r=m%a) and are needed for Schrage's algorithm.
     */
    long a1 = 40014;
    long m1 = 2147483563;
    long q1 = 53668;
    long r1 = 12211;

    long a2 = 40692;
    long m2 = 2147483399;
    long q2 = 52774;
    long r2 = 3791;

    if (sre_randseed > 0)
    {
        rnd1 = sre_randseed;
        rnd2 = sre_randseed;
        /* Fill the table for Bays/Durham */
        for (i = 0; i < 64; i++) {
            x    = a1*(rnd1%q1);   /* LCG1 in action... */
            y    = r1*(rnd1/q1);
            rnd1 = x-y;
            if (rnd1 < 0) rnd1 += m1;

            x    = a2*(rnd2%q2);   /* LCG2 in action... */
            y    = r2*(rnd2/q2);
            rnd2 = x-y;
            if (rnd2 < 0) rnd2 += m2;

            tbl[i] = rnd1-rnd2;
            if (tbl[i] < 0) tbl[i] += m1;
        }
        sre_randseed = 0;		/* drop the flag. */
    }/* end of initialization*/


    x    = a1*(rnd1%q1);   /* LCG1 in action... */
    y    = r1*(rnd1/q1);
    rnd1 = x-y;
    if (rnd1 < 0) rnd1 += m1;

    x    = a2*(rnd2%q2);   /* LCG2 in action... */
    y    = r2*(rnd2/q2);
    rnd2 = x-y;
    if (rnd2 < 0) rnd2 += m2;

    /* Choose our random number from the table... */
    i   = (int) (((double) rnd / (double) m1) * 64.);
    rnd = tbl[i];
    /* and replace with a new number by L'Ecuyer. */
    tbl[i] = rnd1-rnd2;
    if (tbl[i] < 0) tbl[i] += m1;

    return ((double) rnd / (double) m1);
}



/* Function: FreeHistogram()
 *
 * Purpose:  free a histogram structure.
 */
void
FreeHistogram(struct histogram_s *h)
{
    free(h->histogram);
    if (h->expect != NULL) free(h->expect);
    free(h);
}

/* Function: UnfitHistogram()
 *
 * Purpose:  Free only the theoretical fit part of a histogram.
 */
void
UnfitHistogram(struct histogram_s *h)
{
    if (h->expect != NULL) free(h->expect);
    h->expect   = NULL;
    h->fit_type = HISTFIT_NONE;
}



/* Function: Lawless416()
 * Date:     SRE, Thu Nov 13 11:48:50 1997 [St. Louis]
 *
 * Purpose:  Equation 4.1.6 from [Lawless82], pg. 143, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to EVD lambda parameter.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *
 *           Can either deal with a histogram or an array.
 *
 *           Warning: beware overflow/underflow issues! not bulletproof.
 *
 * Args:     x      - array of sample values (or x-axis of a histogram)
 *           y      - NULL (or y-axis of a histogram)
 *           n      - number of samples (or number of histogram bins)
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.1.6 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.1.6 evaluated at lambda
 *
 * Return:   (void)
 */
void
Lawless416(float *x, int *y, int n, float lambda, float *ret_f, float *ret_df)
{

    double esum;			/* \sum e^(-lambda xi)      */
    double xesum;			/* \sum xi e^(-lambda xi)   */
    double xxesum;		/* \sum xi^2 e^(-lambda xi) */
    double xsum;			/* \sum xi                  */
    double mult;			/* histogram count multiplier */
    double total;			/* total samples            */
    int i;


    esum = xesum = xsum  = xxesum = total = 0.;
    for (i = 0; i < n; i++)
    {
        mult = (y == NULL) ? 1. : (double) y[i];
        xsum   += mult * x[i];
        xesum  += mult * x[i] * exp(-1. * lambda * x[i]);
        xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
        esum   += mult * exp(-1. * lambda * x[i]);
        total  += mult;
    }
    *ret_f  = 1./lambda - xsum / total + xesum / esum;
    *ret_df = ((xesum / esum) * (xesum / esum))
              - (xxesum / esum)
              - (1. / (lambda * lambda));

    return;
}


/* Function: Lawless422()
 * Date:     SRE, Mon Nov 17 09:42:48 1997 [St. Louis]
 *
 * Purpose:  Equation 4.2.2 from [Lawless82], pg. 169, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to EVD lambda parameter
 *           for Type I censored data.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *
 *           Can either deal with a histogram or an array.
 *
 *           Warning: beware overflow/underflow issues! not bulletproof.
 *
 * Args:     x      - array of sample values (or x-axis of a histogram)
 *           y      - NULL (or y-axis of a histogram)
 *           n      - number of observed samples (or number of histogram bins)
 *           z      - number of censored samples
 *           c      - censoring value; all observed x_i >= c
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.2.2 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.2.2 evaluated at lambda
 *
 * Return:   (void)
 */
void
Lawless422(float *x, int *y, int n, int z, float c,
           float lambda, float *ret_f, float *ret_df)
{
    double esum;			/* \sum e^(-lambda xi)      + z term    */
    double xesum;			/* \sum xi e^(-lambda xi)   + z term    */
    double xxesum;		/* \sum xi^2 e^(-lambda xi) + z term    */
    double xsum;			/* \sum xi                  (no z term) */
    double mult;			/* histogram count multiplier */
    double total;			/* total samples            */
    int i;

    esum = xesum = xsum  = xxesum = total = 0.;
    for (i = 0; i < n; i++)
    {
        mult = (y == NULL) ? 1. : (double) y[i];
        xsum   += mult * x[i];
        esum   += mult *               exp(-1. * lambda * x[i]);
        xesum  += mult * x[i] *        exp(-1. * lambda * x[i]);
        xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
        total  += mult;
    }

    /* Add z terms for censored data
     */
    esum   += (double) z *         exp(-1. * lambda * c);
    xesum  += (double) z * c *     exp(-1. * lambda * c);
    xxesum += (double) z * c * c * exp(-1. * lambda * c);

    *ret_f  = 1./lambda - xsum / total + xesum / esum;
    *ret_df = ((xesum / esum) * (xesum / esum))
              - (xxesum / esum)
              - (1. / (lambda * lambda));

    return;
}






/* Function: EVDMaxLikelyFit()
 * Date:     SRE, Fri Nov 14 07:56:29 1997 [St. Louis]
 *
 * Purpose:  Given a list or a histogram of EVD-distributed samples,
 *           find maximum likelihood parameters lambda and
 *           mu.
 *
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *           for lambda using Newton/Raphson iterations;
 *           then substitutes lambda into Lawless' equation 4.1.5
 *           to get mu.
 *
 *           Newton/Raphson algorithm developed from description in
 *           Numerical Recipes in C [Press88].
 *
 * Args:     x          - list of EVD distributed samples or x-axis of histogram
 *           c          - NULL, or y-axis of histogram
 *           n          - number of samples, or number of histogram bins
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *
 * Return:   1 on success; 0 on any failure
 */
int
EVDMaxLikelyFit(float *x, int *c, int n, float *ret_mu, float *ret_lambda)
{
    float  lambda, mu;
    float  fx;			/* f(x)  */
    float  dfx;			/* f'(x) */
    double esum;                  /* \sum e^(-lambda xi) */
    double mult;
    double total;
    float  tol = 1e-5;
    int    i;

    /* 1. Find an initial guess at lambda: linear regression here?
     */
    lambda = 0.2;

    /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
     */
    for (i = 0; i < 100; i++)
    {
        Lawless416(x, c, n, lambda, &fx, &dfx);
        if (fabs(fx) < tol) break;             /* success */
        lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
        if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

    /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
     *      Resort to a bisection search. Worse convergence speed
     *      but guaranteed to converge (unlike Newton/Raphson).
     *      We assume (!?) that fx is a monotonically decreasing function of x;
     *      i.e. fx > 0 if we are left of the root, fx < 0 if we
     *      are right of the root.
     */
    if (i == 100)
    {
        float left, right, mid;
        printf(("EVDMaxLikelyFit(): Newton/Raphson failed; switchover to bisection"));

        /* First we need to bracket the root */
        lambda = right = left = 0.2;
        Lawless416(x, c, n, lambda, &fx, &dfx);
        if (fx < 0.)
        {			/* fix right; search left. */
            do {
                left -= 0.1;
                if (left < 0.) {
                    printf(("EVDMaxLikelyFit(): failed to bracket root"));
                    return 0;
                }
                Lawless416(x, c, n, left, &fx, &dfx);
            } while (fx < 0.);
        }
        else
        {			/* fix left; search right. */
            do {
                right += 0.1;
                Lawless416(x, c, n, right, &fx, &dfx);
                if (right > 100.) {
                    printf(("EVDMaxLikelyFit(): failed to bracket root"));
                    return 0;
                }
            } while (fx > 0.);
        }
        /* now we bisection search in left/right interval */
        for (i = 0; i < 100; i++)
        {
            mid = (left + right) / 2.;
            Lawless416(x, c, n, mid, &fx, &dfx);
            if (fabs(fx) < tol) break;             /* success */
            if (fx > 0.)	left = mid;
            else          right = mid;
        }
        if (i == 100) {
            printf(("EVDMaxLikelyFit(): even the bisection search failed"));
            return 0;
        }
        lambda = mid;
    }

    /* 3. Substitute into Lawless 4.1.5 to find mu
     */
    esum = 0.;
    total = 0.;
    for (i = 0; i < n; i++)
    {
        mult   = (c == NULL) ? 1. : (double) c[i];
        esum  += mult * exp(-1 * lambda * x[i]);
        total += mult;
    }
    mu = -1. * log(esum / total) / lambda;

    *ret_lambda = lambda;
    *ret_mu     = mu;
    return 1;
}

/* Function: Gammln()
 *
 * Returns the natural log of the gamma function of x.
 * x is > 0.0.
 *
 * Adapted from a public domain implementation in the
 * NCBI core math library. Thanks to John Spouge and
 * the NCBI. (According to the NCBI, that's Dr. John
 * "Gammas Galore" Spouge to you, pal.)
 */
double
Gammln(double x)
{
    int i;
    double xx, tx;
    double tmp, value;
    static double cof[11] = {
            4.694580336184385e+04,
            -1.560605207784446e+05,
            2.065049568014106e+05,
            -1.388934775095388e+05,
            5.031796415085709e+04,
            -9.601592329182778e+03,
            8.785855930895250e+02,
            -3.155153906098611e+01,
            2.908143421162229e-01,
            -2.319827630494973e-04,
            1.251639670050933e-10
    };

    /* Protect against x=0. We see this in Dirichlet code,
     * for terms alpha = 0. This is a severe hack but it is effective
     * and (we think?) safe. (due to GJM)
     */
    if (x <= 0.0) return 999999.;

    xx       = x - 1.0;
    tx = tmp = xx + 11.0;
    value    = 1.0;
    for (i = 10; i >= 0; i--)	/* sum least significant terms first */
    {
        value += cof[i] / tmp;
        tmp   -= 1.0;
    }
    value  = log(value);
    tx    += 0.5;
    value += 0.918938533 + (xx+0.5)*log(tx) - tx;
    return value;
}

/* Function: IncompleteGamma()
 *
 * Purpose:  Returns 1 - P(a,x) where:
 *           P(a,x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt
 *                  = \frac{\gamma(a,x)}{\Gamma(a)}
 *                  = 1 - \frac{\Gamma(a,x)}{\Gamma(a)}
 *
 *           Used in a chi-squared test: for a X^2 statistic x
 *           with v degrees of freedom, call:
 *                  p = IncompleteGamma(v/2., x/2.)
 *           to get the probability p that a chi-squared value
 *           greater than x could be obtained by chance even for
 *           a correct model. (i.e. p should be large, say
 *           0.95 or more).
 *
 * Method:   Based on ideas from Numerical Recipes in C, Press et al.,
 *           Cambridge University Press, 1988.
 *
 * Args:     a  - for instance, degrees of freedom / 2     [a > 0]
 *           x  - for instance, chi-squared statistic / 2  [x >= 0]
 *
 * Return:   1 - P(a,x).
 */
double
IncompleteGamma(double a, double x)
{
    int iter;			/* iteration counter */

    if (a <= 0.) exit(-1);
    if (x <  0.) exit(-1);

    /* For x > a + 1 the following gives rapid convergence;
     * calculate 1 - P(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)}:
     *     use a continued fraction development for \Gamma(a,x).
     */
    if (x > a+1)
    {
        double oldp;		/* previous value of p    */
        double nu0, nu1;		/* numerators for continued fraction calc   */
        double de0, de1;		/* denominators for continued fraction calc */

        nu0 = 0.;			/* A_0 = 0       */
        de0 = 1.;			/* B_0 = 1       */
        nu1 = 1.;			/* A_1 = 1       */
        de1 = x;			/* B_1 = x       */

        oldp = nu1;
        for (iter = 1; iter < 100; iter++)
        {
            /* Continued fraction development:
             * set A_j = b_j A_j-1 + a_j A_j-2
             *     B_j = b_j B_j-1 + a_j B_j-2
                 * We start with A_2, B_2.
             */
            /* j = even: a_j = iter-a, b_j = 1 */
            /* A,B_j-2 are in nu0, de0; A,B_j-1 are in nu1,de1 */
            nu0 = nu1 + ((double)iter - a) * nu0;
            de0 = de1 + ((double)iter - a) * de0;

            /* j = odd: a_j = iter, b_j = x */
            /* A,B_j-2 are in nu1, de1; A,B_j-1 in nu0,de0 */
            nu1 = x * nu0 + (double) iter * nu1;
            de1 = x * de0 + (double) iter * de1;

            /* rescale */
            if (de1 != 0.)
            {
                nu0 /= de1;
                de0 /= de1;
                nu1 /= de1;
                de1 =  1.;
            }
            /* check for convergence */
            if (fabs((nu1-oldp)/nu1) < 1.e-7)
                return nu1 * exp(a * log(x) - x - Gammln(a));

            oldp = nu1;
        }
        exit(-1);
    }
    else /* x <= a+1 */
    {
        double p;			/* current sum               */
        double val;		/* current value used in sum */

        /* For x <= a+1 we use a convergent series instead:
         *   P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)},
         * where
         *   \gamma(a,x) = e^{-x}x^a \sum_{n=0}{\infty} \frac{\Gamma{a}}{\Gamma{a+1+n}} x^n
         * which looks appalling but the sum is in fact rearrangeable to
         * a simple series without the \Gamma functions:
         *   = \frac{1}{a} + \frac{x}{a(a+1)} + \frac{x^2}{a(a+1)(a+2)} ...
         * and it's obvious that this should converge nicely for x <= a+1.
         */

        p = val = 1. / a;
        for (iter = 1; iter < 10000; iter++)
        {
            val *= x / (a+(double)iter);
            p   += val;

            if (fabs(val/p) < 1.e-7)
                return 1. - p * exp(a * log(x) - x - Gammln(a));
        }
        exit(-1);
    }
    /*NOTREACHED*/
    return 0.;
}


int samplemulambda(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbrAA(par.db1, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader qdbr(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);

    IndexReader *t3DiDbr = NULL;
    IndexReader *tAADbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &qdbr;
        tAADbr = &qdbrAA;
    } else {
        tAADbr = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    }


    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    SubstitutionMatrix subMatAA(blosum.c_str(), 1.4, par.scoreBias);

    //temporary output file
    Debug::Progress progress(tAADbr->sequenceReader->getSize());

    // sub. mat needed for query profile
    int8_t * tinySubMatAA = (int8_t*) mem_align(ALIGN_INT, subMatAA.alphabetSize * 32);
    int8_t * tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat3Di.alphabetSize * 32);

    for (int i = 0; i < subMat3Di.alphabetSize; i++) {
        for (int j = 0; j < subMat3Di.alphabetSize; j++) {
            tinySubMat3Di[i * subMat3Di.alphabetSize + j] = subMat3Di.subMatrix[i][j]; // for farrar profile
        }
    }
    for (int i = 0; i < subMatAA.alphabetSize; i++) {
        for (int j = 0; j < subMatAA.alphabetSize; j++) {
            tinySubMatAA[i * subMatAA.alphabetSize + j] = subMatAA.subMatrix[i][j];
        }
    }
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        EvalueNeuralNet evaluer(tAADbr->sequenceReader->getAminoAcidDBSize(), &subMat3Di);
        std::vector<Matcher::result_t> alignmentResult;
        StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale);
        StructureSmithWaterman reverseStructureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale);

        Sequence qSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        std::string backtrace;
        std::string resultBuffer;
        // write output file
        std::mt19937 rnd(0);

#pragma omp for schedule(dynamic, 1)
        for(size_t id = 0; id < qdbrAA.sequenceReader->getSize(); id++) {
            progress.updateProgress();
            unsigned int queryKey = qdbr.sequenceReader->getDbKey(id);
            unsigned int querySeqLen = qdbr.sequenceReader->getSeqLen(id);
            char *querySeqAA = qdbrAA.sequenceReader->getData(id, thread_idx);
            char *querySeq3Di = qdbr.sequenceReader->getData(id, thread_idx);
            qSeq3Di.mapSequence(id, queryKey, querySeq3Di, querySeqLen);
            qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);
            structureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
            qSeq3Di.reverse();
            qSeqAA.reverse();
            reverseStructureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
            for (int sample = 0; sample < par.nsample; sample++) {
                // pick random number between 0 and size of database
                size_t sampleIdx = static_cast<size_t>(t3DiDbr->sequenceReader->getSize() * drand48());
                const unsigned int dbKey = t3DiDbr->sequenceReader->getDbKey(sampleIdx);
                unsigned int targetId = t3DiDbr->sequenceReader->getId(dbKey);

                char * targetSeq3Di = t3DiDbr->sequenceReader->getData(targetId, thread_idx);
                char * targetSeqAA = tAADbr->sequenceReader->getData(targetId, thread_idx);
                const int targetLen = static_cast<int>(t3DiDbr->sequenceReader->getSeqLen(targetId));
                tSeq3Di.mapSequence(targetId, dbKey, targetSeq3Di, targetLen);
                tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetLen);
                // shuffle a vector of integers
                std::vector<int> indices(targetLen);
                for (int i = 0; i < targetLen; i++) {
                    indices[i] = i;
                }
                std::shuffle(indices.begin(), indices.end(), rnd);
                for (int i = 0; i < targetLen; i++) {
                    std::swap(tSeq3Di.numSequence[i], tSeq3Di.numSequence[indices[i]]);
                    std::swap(tSeqAA.numSequence[i], tSeqAA.numSequence[indices[i]]);
                }

                StructureSmithWaterman::s_align align = structureSmithWaterman.alignScoreEndPos(tSeqAA.numSequence, tSeq3Di.numSequence, targetLen, par.gapOpen.values.aminoacid(),
                                                                                                par.gapExtend.values.aminoacid(), querySeqLen / 2);
                StructureSmithWaterman::s_align revAlign = reverseStructureSmithWaterman.alignScoreEndPos(tSeqAA.numSequence, tSeq3Di.numSequence, targetLen, par.gapOpen.values.aminoacid(),
                                                                                                          par.gapExtend.values.aminoacid(), querySeqLen / 2);
                int32_t score = static_cast<int32_t>(align.score1) - static_cast<int32_t>(revAlign.score1);
                align.evalue = 0.0;
                //align.evalue = evaluer.computeEvalue(score, muLambda.first, muLambda.second);

                //align = structureSmithWaterman.alignStartPosBacktrace(tSeqAA.numSequence, tSeq3Di.numSequence, targetLen, par.gapOpen.values.aminoacid(),
                //                                                      par.gapExtend.values.aminoacid(), par.alignmentMode, backtrace,  align, par.covMode, par.covThr, querySeqLen / 2);

                unsigned int alnLength = Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1, align.dbStartPos1, align.dbEndPos1);
                float seqId = Util::computeSeqId(par.seqIdMode, align.identicalAACnt, querySeqLen, targetLen, alnLength);
                align.score1 = score;
                Matcher::result_t res(dbKey, align.score1, align.qCov, align.tCov, seqId, align.evalue, alnLength, align.qStartPos1, align.qEndPos1, querySeqLen, align.dbStartPos1, align.dbEndPos1, targetLen, backtrace);
                alignmentResult.emplace_back(res);

            }
            if (alignmentResult.size() > 1) {
                SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), Matcher::compareHits);
            }
            float * scores = new float[alignmentResult.size()];
            for(size_t i = 0; i < alignmentResult.size(); i++) {
                scores[i] = alignmentResult[i].score;
            }
            float mu = 0.0;
	    float lambda = 0.0;
            EVDMaxLikelyFit(scores, NULL, alignmentResult.size(), &mu, &lambda);
            delete [] scores;
            resultBuffer.append(querySeqAA, querySeqLen);
	        resultBuffer.push_back('\t');
            resultBuffer.append(querySeq3Di, querySeqLen);
            resultBuffer.push_back('\t');
            resultBuffer.append(SSTR(mu));
            resultBuffer.push_back('\t');
            resultBuffer.append(SSTR(lambda));
            resultBuffer.push_back('\n');
            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            alignmentResult.clear();
        }
    }

    free(tinySubMatAA);
    free(tinySubMat3Di);
    dbw.close();
    if (sameDB == false) {
        delete t3DiDbr;
        delete tAADbr;
    }
    return EXIT_SUCCESS;
}

