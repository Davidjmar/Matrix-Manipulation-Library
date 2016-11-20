#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "matrix.h"

int main(void) {
  matrix * R, * GJ, * Ad, * Ai, * I;
  int err;

  R = newMatrix(6, 12); assert(R);
  err = setElement(R, 2, 1, 2.0); assert(!err);
  err = setElement(R, 1, 2, 1.0); assert(!err);
  err = setElement(R, 3, 2, 3.0); assert(!err);
  err = setElement(R, 2, 3, 2.0); assert(!err);
  err = setElement(R, 4, 3, 1.0); assert(!err);
  err = setElement(R, 3, 4, 1.0); assert(!err);
  err = setElement(R, 5, 4, 3.0); assert(!err);
  err = setElement(R, 4, 5, 2.0); assert(!err);
  err = setElement(R, 6, 5, 2.0); assert(!err);
  err = setElement(R, 5, 6, 1.0); assert(!err);

  err = setElement(R, 1, 7,  1.0); assert(!err);
  err = setElement(R, 2, 8,  1.0); assert(!err);
  err = setElement(R, 3, 9,  1.0); assert(!err);
  err = setElement(R, 4, 10, 1.0); assert(!err);
  err = setElement(R, 5, 11, 1.0); assert(!err);
  err = setElement(R, 6, 12, 1.0); assert(!err);

  GJ = newMatrix(6, 12); assert(GJ);
  err = gaussJordan(R, GJ); assert(!err);
  Ad = newMatrix(6, 6); assert(Ad);
  Ai = newMatrix(6, 6); assert(Ai);
  I = newMatrix(6, 6); assert(I);
  err = submatrix(R, Ad, 1, 1); assert(!err);
  err = submatrix(GJ, Ai, 1, 7); assert(!err);
  printf("\nAd:\n");
  err = printMatrix(Ad); assert(!err);
  printf("\nAi:\n");
  err = printMatrix(Ai); assert(!err);
  err = product(Ad, Ai, I); assert(!err);
  if (isIdentity(I)) {
    printf("\nAi is the inverse of Ad\n");
  } else {
    /* The problem may be rounding. Let's look at the result. */
    printf("\nAi is not the inverse of Ad.  Their product is:\n");
    err = printMatrix(I); assert(!err);
  }
  deleteMatrix(Ad);
  deleteMatrix(Ai);
  deleteMatrix(I);
  deleteMatrix(R);
  deleteMatrix(GJ);

  return 0;
}
