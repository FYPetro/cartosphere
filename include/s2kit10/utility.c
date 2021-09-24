
#define _USE_MATH_DEFINES
#include "utility.h"
#include <math.h>

#define compmult(a,b,c,d,e,f) (e) = ((a)*(c))-((b)*(d)); (f) = ((a)*(d))+((b)*(c))

/*****************************************************************

   Given bandwidth bw, seanindex(m,l,bw) will give the position of the
   coefficient f-hat(m,l) in the one-row array that Sean stores the spherical
   coefficients. This is needed to help preserve the symmetry that the
   coefficients have: (l = degree, m = order, and abs(m) <= l)
   
   f-hat(l,-m) = (-1)^m * conjugate( f-hat(l,m) )
   
   Thanks for your help Mark!
   
   ******************************************************************/


int seanindex(int m,
    int l,
    int bw)
{
    int bigL;

    bigL = bw - 1;

    if (m >= 0)
        return(m * (bigL + 1) - ((m * (m - 1)) / 2) + (l - m));
    else
        return(((bigL * (bigL + 3)) / 2) + 1 +
        ((bigL + m) * (bigL + m + 1) / 2) + (l - abs(m)));
}


/************************************************************************/
/*
   multiplies harmonic coefficients of a function and a filter.
   See convolution theorem of Driscoll and Healy for details.
   
   bw -> bandwidth of problem
   size = 2*bw
   
   datacoeffs should be output of an FST, filtercoeffs the
   output of an FZT.  There should be (bw * bw) datacoeffs,
   and bw filtercoeffs.
   rres and ires should point to arrays of dimension bw * bw.
   
   */

void TransMult(double *rdatacoeffs, double *idatacoeffs,
    double *rfiltercoeffs, double *ifiltercoeffs,
    double *rres, double *ires,
    int bw)
{

    int m, l, size;
    double *rdptr, *idptr, *rrptr, *irptr;

    size = 2 * bw;

    rdptr = rdatacoeffs;
    idptr = idatacoeffs;
    rrptr = rres;
    irptr = ires;

    for (m = 0; m<bw; m++) {
        for (l = m; l<bw; l++) {
            compmult(rfiltercoeffs[l], ifiltercoeffs[l],
                rdptr[l - m], idptr[l - m],
                rrptr[l - m], irptr[l - m]);

            rrptr[l - m] *= sqrt(4 * M_PI / (2 * l + 1));
            irptr[l - m] *= sqrt(4 * M_PI / (2 * l + 1));

        }
        rdptr += bw - m; idptr += bw - m;
        rrptr += bw - m; irptr += bw - m;
    }
    for (m = bw + 1; m<size; m++) {
        for (l = size - m; l<bw; l++) {
            compmult(rfiltercoeffs[l], ifiltercoeffs[l],
                rdptr[l - size + m], idptr[l - size + m],
                rrptr[l - size + m], irptr[l - size + m]);

            rrptr[l - size + m] *= sqrt(4 * M_PI / (2 * l + 1));
            irptr[l - size + m] *= sqrt(4 * M_PI / (2 * l + 1));

        }
        rdptr += m - bw; idptr += m - bw;
        rrptr += m - bw; irptr += m - bw;
    }

}
