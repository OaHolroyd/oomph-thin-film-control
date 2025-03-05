#ifndef CONTROL_LQR_H
#define CONTROL_LQR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

#include "c-utils.h"
#include "control-core.h"
#include "linalg.h"

static double **LQR_K; /* control operator */

/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* compute the control matrix in the Benney case */
void lqr_benney_compute_K(double **lqr_k) {
  /* Jacobian */
  double **A = malloc_f2d(NX * NY, NX * NY);
  benney_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(NX * NY, CNTL_M);
  benney_actuator(B);

  /* control matrix */
  dlqr(A, B, DX * DY * MU, 1 - MU, NX * NY, CNTL_M, lqr_k);

  free_2d(A);
  free_2d(B);
}

/* compute the control matrix in the weighted-residuals case */
void lqr_wr_compute_K(double **lqr_k) {
  /* Jacobian */
  double **A = malloc_f2d(2 * NX * NY, 2 * NX * NY);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2 * NX * NY, CNTL_M);
  wr_actuator(B);

  /* full control matrix */
  dlqr(A, B, DX * DY * MU, 1 - MU, 2 * NX * NY, CNTL_M, lqr_k);

  free_2d(A);
  free_2d(B);
}

/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void lqr_set(void) {
  /* pick from the available ROMs */
  switch (RT) {
  case BENNEY:
    LQR_K = malloc_f2d(CNTL_M, NX * NY);
    lqr_benney_compute_K(LQR_K);
    break;
  case WR:
    LQR_K = malloc_f2d(CNTL_M, 2 * NX * NY);
    lqr_wr_compute_K(LQR_K);
    break;
  default:
    ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void lqr_free(void) { free(LQR_K); }

/* [REQUIRED] steps the system forward in time given the interfacial height */
void lqr_step(double dt, double *h, double *q) {
  /* f = K * (h-1) */
  for (int k = 0; k < CNTL_M; k++) {
    Amag[k] = 0.0;
    for (int i = 0; i < NY; i++) {
      for (int j = 0; j < NX; j++) {
        Amag[k] += LQR_K[k][IJTOK(i, j)] * (interp(JTOX(j), ITOY(i), h) - 1.0);
      } // j end
    } // i end

    /* only WR uses the flux */
    if (RT == WR) {
      for (int i = 0; i < NY; i++) {
        for (int j = 0; j < NX; j++) {
          Amag[k] += LQR_K[k][NX * NY + IJTOK(i, j)] *
                     (interp(JTOX(j), ITOY(i), q) - 2.0 / 3.0);
        } // j end
      } // i end
    }
  }
}

/* [REQUIRED] returns the estimator as a function of x */
double lqr_estimator(double x, double y) { return 0.0; }

/* [REQUIRED] outputs the internal matrices */
void lqr_output(void) {
  if (RT == BENNEY) {
    output_d2d("out/K.dat", LQR_K, CNTL_M, NX * NY);
  } else {
    output_d2d("out/K.dat", LQR_K, CNTL_M, 2 * NX * NY);
  }
}

/* [REQUIRED] generates the control matrix CM = F*K */
void lqr_matrix(double **CM) {
  /* forcing matrix */
  double **F = malloc_f2d(NX * NY, CNTL_M);
  forcing_matrix(F);

  for (int i = 0; i < NX * NY; i++) {
    for (int j = 0; j < 2 * NX * NY; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < CNTL_M; k++) {
        CM[i][j] += F[i][k] * LQR_K[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#ifdef __cplusplus
}
#endif

#endif
