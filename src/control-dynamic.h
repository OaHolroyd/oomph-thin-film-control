#ifndef CONTROL_DYNAMIC_H
#define CONTROL_DYNAMIC_H

#ifdef __cplusplus
extern "C" {
#endif


#include <math.h>
#include <complex.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"


static double complex **DYNAMIC_A; // dynamics (linearised system)
static double complex **DYNAMIC_B; // control (linearised system)
static double complex **DYNAMIC_K; // control operator
static double complex *DYNAMIC_z; // estimator

static double complex **DYNAMIC_A1; // inversion matrix
static double complex *DYNAMIC_z0; // previous timestep
static int *DYNAMIC_ipiv; // pivot array


static void (*dynamic_update)(double dt, double *h); // update function


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* sets the various matrices required for the Benney system to work */
void dynamic_benney_set(void) {
  int i, j, k;


  /* compute (real) observer matrix (actually the transpose) */
  double **rPhi = malloc_f2d(N, P);
  benney_observer(rPhi);


  /* compute controllable wavenumbers */
  double *K = malloc(M*sizeof(double));
  j = 1;
  for (i = 0; i < M; i++) {
    j *= -1;
    k = j*(i+1)/2; // uses integer rounding
    K[i] = (2*M_PI/LX)*k;
  } // i end


  /* compute spectral Jacobian */
  double complex **J = malloc_z2d(M, M);
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      J[i][j] = 0.0;
    } // j end
    J[i][i] = -2*I*K[i]
              -(2.0/3.0/tan(THETA) - 8.0*RE/15.0)*K[i]*K[i]
              -1.0/3.0/CA*K[i]*K[i]*K[i]*K[i];
  } // i end


  /* compute spectral Psi - actuator */
  double complex **Psi = malloc_z2d(M, M);
  double xk;
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      Psi[i][j] = 0.0;

      /* manually compute the DFT via dot product */
      for (k = 0; k < N; k++) {
        xk = ITOX(k);
        Psi[i][j] += actuator(xk-Aloc[j]) * (cos(K[i]*xk) - I*sin(K[i]*xk))/N;
      } // k end

      /* (1 + 2RE/3 dx) component */
      Psi[i][j] *= 1.0 + 2.0*RE/3.0*I*K[i];
    } // j end
  } // i end


  /* compute spectral K */
  zlqr(J, Psi, MU*LX, (1-MU), M, M, DYNAMIC_K);


  /* compute spectral Phi - observer (this is actually the transpose) */
  double complex **Phi = malloc_z2d(M, P);
  for (i = 0; i < M; i++) {
    for (j = 0; j < P; j++) {
      Phi[i][j] = 0.0;

      /* manually compute the inverse DFT via dot product */
      for (k = 0; k < N; k++) {
        xk = ITOX(k);
        Phi[i][j] += rPhi[k][j] * (cos(K[i]*xk) + I*sin(K[i]*xk));
      } // k end

    } // j end
  } // i end


  /* compute L (this is actually the transpose) */
  const double DYNAMIC_MU_L = 0.5; // TODO: choose this optimally
  double complex **L = malloc_z2d(P, M);
  zlqr(J, Phi, DYNAMIC_MU_L*LX, DYNAMIC_MU_L, M, P, L);


  /* compute A */
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      DYNAMIC_A[i][j] = J[i][j];
      for (k = 0; k < M; k++) {
        DYNAMIC_A[i][j] += Psi[i][k] * DYNAMIC_K[k][j];
      } // k end
      for (k = 0; k < P; k++) {
        DYNAMIC_A[i][j] -= L[k][i] * Phi[j][k];
      } // k end
    } // j end
  } // i end


  /* compute B */
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      DYNAMIC_B[i][j] = 0.0;
      for (k = 0; k < P; k++) {
        DYNAMIC_B[i][j] += L[k][i] * rPhi[j][k];
      } // k end
    } // j end
  } // i end


  free_2d(rPhi);
  free(K);
  free_2d(J);
  free_2d(Psi);
  free_2d(Phi);
  free_2d(L);
}

/* solves the linearised system z_t = Az + B(h-1) forward by one timestep */
// TODO: this currently uses C-N but Euler looks like it would be fine too
void dynamic_benney_update_estimator(double dt, double *h) {
  int i, j;

  /* store current estimator */
  for (i = 0; i < M; i++) {
    DYNAMIC_z0[i] = DYNAMIC_z[i];
  } // i end


  /* compute the explicit half */
  for (i = 0; i < M; i++) {
    /* linearised dynamics */
    for (j = 0; j < M; j++) {
      DYNAMIC_z[i] += 0.5*dt*DYNAMIC_A[i][j]*DYNAMIC_z0[j];
    } // j end

    /* linearised control */
    for (j = 0; j < N; j++) {
      DYNAMIC_z[i] += dt*DYNAMIC_B[i][j]*(interp(ITOX(j), h) - 1.0);
    } // j end
  } // i end


  /* only need to recompute the implicit matrix if dt has changed */
  static double dt_prev = -100000;
  if (fabs(dt-dt_prev) > 1e-15) {
    /* set implicit matrix */
    for (i = 0; i < M; i++) {
      for (j = 0; j < M; j++) {
        DYNAMIC_A1[i][j] = -0.5*dt*DYNAMIC_A[i][j];
      } // j end
      DYNAMIC_A1[i][i] += 1.0;
    } // i end

    /* compute all-in-one LU factorisation */
    zlu(DYNAMIC_A1, M);
  }
  dt_prev = dt;

  /* solve LUx = z */
  zlusv(DYNAMIC_A1, DYNAMIC_z, M);
}

/* sets the various matrices required for the WR system to work */
void dynamic_wr_set(void) {
  int i, j, k;


  /* observer matrix (actually the transpose) */
  double **rPhi = malloc_f2d(2*N, P);
  wr_observer(rPhi);


  /* compute controllable wavenumbers */
  double *K = malloc(M*sizeof(double));
  j = 1;
  for (i = 0; i < M; i++) {
    j *= -1;
    k = j*(i+1)/2; // uses integer rounding
    K[i] = (2*M_PI/LX)*k;
  } // i end


  /* compute spectral Jacobian */
  double complex **J = malloc_z2d(2*M, 2*M);
  for (i = 0; i < 2*M; i++) {
    for (j = 0; j < 2*M; j++) {
      J[i][j] = 0.0;
    } // j end
  } // i end
  for (i = 0; i < M; i++) {
    J[i][M+i] = -(I*K[i]);
    J[M+i][i] = (5.0/RE) + (4.0/7.0 - 5.0/(3.0*RE*tan(THETA)))*(I*K[i]) + (5.0/(6.0*RE*CA))*(-I*K[i]*K[i]*K[i]);
    J[M+i][M+i] = -(5.0/(2.0*RE)) - (34.0/21.0)*(I*K[i]);
  } // i end


  /* compute spectral Psi - actuator */
  double complex **Psi = malloc_z2d(2*M, M);
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      Psi[i][j] = 0.0;

      /* manually compute the DFT via dot product */
      for (k = 0; k < N; k++) {
        double xk = ITOX(k);
        Psi[i][j] += actuator(xk-Aloc[j]) * (cos(K[i]*xk) - I*sin(K[i]*xk))/N;
      } // k end

      /* q component */
      Psi[i+M][j] = Psi[i][j] * (2.0*RE/15.0);
    } // j end
  } // i end


  /* compute spectral K */
  zlqr(J, Psi, MU*LX, (1.0-MU), 2*M, M, DYNAMIC_K);

  // /* set K to zero */
  // for (i = 0; i < M; i++) {
  //   for (j = 0; j < 2*M; j++) {
  //     DYNAMIC_K[i][j] = 0.0;
  //   } // j end
  // } // i end


  /* compute spectral Phi - observer (this is actually the transpose) */
  double complex **Phi = malloc_z2d(2*M, P);
  for (i = 0; i < M; i++) {
    for (j = 0; j < P; j++) {
      Phi[i][j] = 0.0;
      Phi[i+M][j] = 0.0;

      /* manually compute the inverse DFT via dot product */
      for (k = 0; k < N; k++) {
        double xk = ITOX(k);
        Phi[i][j] += rPhi[k][j] * (cos(K[i]*xk) + I*sin(K[i]*xk));
        Phi[i+M][j] += rPhi[k+N][j] * (cos(K[i]*xk) + I*sin(K[i]*xk));
      } // k end

    } // j end
  } // i end


  /* compute L (this is actually the transpose) */
  const double DYNAMIC_MU_L = 0.9; // TODO: choose this optimally
  double complex **L = malloc_z2d(P, 2*M);
  zlqr(J, Phi, (1.0-DYNAMIC_MU_L)*LX, DYNAMIC_MU_L, 2*M, P, L);


  /* compute A */
  for (i = 0; i < 2*M; i++) {
    for (j = 0; j < 2*M; j++) {
      DYNAMIC_A[i][j] = J[i][j];
      for (k = 0; k < M; k++) {
        DYNAMIC_A[i][j] -= Psi[i][k] * DYNAMIC_K[k][j]; // NOTE opposite sign here
      } // k end
      for (k = 0; k < P; k++) {
        DYNAMIC_A[i][j] -= L[k][i] * Phi[j][k];
      } // k end
    } // j end
  } // i end


  /* compute B */
  for (i = 0; i < 2*M; i++) {
    for (j = 0; j < N; j++) {
      DYNAMIC_B[i][j] = 0.0;
      for (k = 0; k < P; k++) {
        DYNAMIC_B[i][j] += L[k][i] * rPhi[j][k];
      } // k end
    } // j end
  } // i end


  free_2d(rPhi);
  free(K);
  free_2d(J);
  free_2d(Psi);
  free_2d(Phi);
  free_2d(L);
}

/* solves the linearised system z_t = Az + B(h-1) forward by one timestep */
// TODO: this currently uses C-N but Euler looks like it would be fine too
void dynamic_wr_update_estimator(double dt, double *h) {
  /* store current estimator */
  for (int i = 0; i < 2*M; i++) {
    DYNAMIC_z0[i] = DYNAMIC_z[i];
  } // i end

  /* Euler forward step */
  for (int i = 0; i < 2*M; i++) {
    /* linearised dynamics */
    for (int j = 0; j < 2*M; j++) {
      DYNAMIC_z[i] += dt*DYNAMIC_A[i][j]*DYNAMIC_z0[j];
    } // j end

    /* linearised control */
    for (int j = 0; j < N; j++) {
      DYNAMIC_z[i] += dt*DYNAMIC_B[i][j]*(interp(ITOX(j), h) - 1.0);
    } // j end
  } // i end
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void dynamic_set(void) {
  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      DYNAMIC_K = malloc_z2d(M, M);
      DYNAMIC_B = malloc_z2d(M, N);
      DYNAMIC_A = malloc_z2d(M, M);
      DYNAMIC_z = malloc(M*sizeof(double complex));
      for (int i = 0; i < M; i++) {
        DYNAMIC_z[i] = 0.0;
      } // i end

      DYNAMIC_A1 = malloc_z2d(M, M);
      DYNAMIC_z0 = malloc(M*sizeof(double complex));
      DYNAMIC_ipiv = malloc(M*sizeof(int));

      dynamic_benney_set();
      dynamic_update = &dynamic_benney_update_estimator;
      break;
    case WR:
      // ABORT("dynamic WR not implemented yet");

      DYNAMIC_K = malloc_z2d(M, 2*M);
      DYNAMIC_B = malloc_z2d(2*M, N);
      DYNAMIC_A = malloc_z2d(2*M, 2*M);
      DYNAMIC_z = malloc(2*M*sizeof(double complex));
      for (int i = 0; i < 2*M; i++) {
        DYNAMIC_z[i] = 0.0;
      } // i end

      DYNAMIC_A1 = malloc_z2d(2*M, 2*M);
      DYNAMIC_z0 = malloc(2*M*sizeof(double complex));
      DYNAMIC_ipiv = malloc(M*sizeof(int)); // TODO: not used

      dynamic_wr_set();
      dynamic_update = &dynamic_wr_update_estimator;
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void dynamic_free(void) {
  free_2d(DYNAMIC_K);
  free_2d(DYNAMIC_B);
  free_2d(DYNAMIC_A);
  free(DYNAMIC_z);

  free_2d(DYNAMIC_A1);
  free(DYNAMIC_z0);
  free(DYNAMIC_ipiv);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void dynamic_step(double dt, double *h, double *q) {
  /* update system */
  dynamic_update(dt, h);

  /* f = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;
    for (int j = 0; j < M; j++) {
      Amag[i] += creal(DYNAMIC_K[i][j] * DYNAMIC_z[j]); // ensure this is real
    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double dynamic_estimator(double x) {
  /* account for periodicity */
  if (x < 0) {
    x += LX;
  } else if (x > LX) {
    x -= LX;
  }

  /* sum the terms from the Fourier series */
  double complex e = 0.0;
  double dk = 2 * M_PI / LX;
  int j, k;
  j = 1;
  for (int i = 0; i < M; i++) {
    j *= -1;
    k = j*(i+1)/2; // uses integer rounding
    e += DYNAMIC_z[i]*(cos(dk*k*x) + I*sin(dk*k*x));
  } // i end

  return creal(e);
}

/* [REQUIRED] outputs the internal matrices */
void dynamic_output(void) {
  output_z2d("out/A.dat", DYNAMIC_A, M, M);
  output_z2d("out/B.dat", DYNAMIC_B, M, N);
  output_z2d("out/K.dat", DYNAMIC_K, M, M);
}

#ifdef __cplusplus
}
#endif

#endif
