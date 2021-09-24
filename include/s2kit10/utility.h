#ifndef _UTILITY_H
#define _UTILITY_H

#if defined __cplusplus
/// Public APIs will have C linkage if __cplusplus is present.
extern "C" {
#endif

int seanindex(int,
        int,
        int);

void TransMult(double *, double *,
        double *, double *,
        double *, double *,
        int);

#if defined __cplusplus
} /// End of C linkage section.
#endif

#endif