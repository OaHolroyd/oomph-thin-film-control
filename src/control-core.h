#ifndef CONTROL_CORE_H
#define CONTROL_CORE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "c-utils.h"
#include "control.h"

/* convert from index to location */
#define JTOX(j) (DX * (j + 0.5))
#define ITOY(i) (DY * (i + 0.5))

/* convert from location to index */
#define XTOJ(x) ((int)(x / DX - 0.5))
#define YTOI(y) ((int)(y / DY - 0.5))

/* convert from double to single index */
#define IJTOK(i, j) ((((i) + NY) % NY) * NX + ((j) + NX) % NX)

/* ========================================================================== */
/*   GLOBAL DEFINITIONS                                                       */
/* ========================================================================== */
/* film constants */
static int NX;       // number of gridcells in x
static int NY;       // number of gridcells in y
static double LX;    // domain length
static double LY;    // domain width
static double DX;    // grid spacing in x
static double DY;    // grid spacing in y
static double RE;    // Reynolds' number
static double CA;    // capillary number
static double THETA; // plate angle

/* control constants */
static int M;        // number of actuators
static int P;        // number of observers
static double W;     // width parameter
static double ALPHA; // control strength parameter
static double MU;    // control cost parameter
static double NORM;  // normalising constant
static double DEL;   // observer offset (upstream)
static rom_t RT;     // type of reduced order model

/* location arrays */
static double *Aloc; // actuator locations (stored as x1, y1, x2, y2, ...)
static double *Oloc; // observer locations (stored as x1, y1, x2, y2, ...)

/* variable arrays */
static double *Amag; // actuator magnitudes

/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* forcing matrix (N-by-M) */
void forcing_matrix(double **F);

/* Jacobian (N-by-N) */
void benney_jacobian(double **A);

/* Actuator (N-by-M) */
void benney_actuator(double **B);

/* (the transpose of the) Observer (N-by-P) */
void benney_observer(double **C);

/* Jacobian (2N-by-2N) */
void wr_jacobian(double **A);

/* Actuator (2N-by-M) */
void wr_actuator(double **B);

/* (the transpose of the) Observer (2N-by-2P) */
void wr_observer(double **C);

/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* interpolates the discretised film height h at a position x */
double interp(double x, double y, double *h) {
  /* account for periodicity */
  if (x < 0) {
    x += LX;
  } else if (x > LX) {
    x -= LX;
  }
  if (y < 0) {
    y += LY;
  } else if (y > LY) {
    y -= LY;
  }

  /* get neighbouring indices */
  double i = y / DY - 0.5;
  double _i0 = floor(i);
  double _i1 = ceil(i);
  int i0 = (int)_i0;
  int i1 = (int)_i1;

  double j = x / DX - 0.5;
  double _j0 = floor(j);
  double _j1 = ceil(j);
  int j0 = (int)_j0;
  int j1 = (int)_j1;

  /* bilinear interpolation */
  double di = i - _i0;
  double dj = j - _j0;
  return h[IJTOK(i0, j0)] * (1 - di) * (1 - dj) +
         h[IJTOK(i0, j1)] * (1 - di) * dj + h[IJTOK(i1, j0)] * di * (1 - dj) +
         h[IJTOK(i1, j1)] * di * dj;
}

/* (periodic) actuator function (only valid on [-3Lx/2, 3Lx/2]) */
double actuator(double x, double y) {
  /* account for periodicity */
  // TODO: since cos is periodic this probably isn't necessary
  if (x < 0) {
    x += LX;
  } else if (x > LX) {
    x -= LX;
  }
  if (y < 0) {
    y += LY;
  } else if (y > LY) {
    y -= LY;
  }

  // rescale x and y to [0, 2pi]
  x = 2 * M_PI * x / LX;
  y = 2 * M_PI * y / LY;

  // 2D control is the product of two 1D controls, one for x, one for y
  //  control(x, y) = exp((cos(x) - 1) / (W * W)) * exp((cos(y) - 1) / (W * W))
  return NORM * exp((cos(x) + cos(y) - 2.0) / (W * W));
}

/* sets the common control parameters and allocates common memory */
void internal_control_set(rom_t rt, int m, int p, double w, double alpha,
                          double mu, double del, double lx, double ly, int nx,
                          int ny, double re, double ca, double theta) {
  /* set constants */
  NX = nx;
  NY = ny;
  LX = lx;
  LY = ly;
  DX = LX / NX;
  DY = LY / NY;
  RE = re;
  CA = ca;
  THETA = theta;

  M = m;
  P = p;
  W = w;
  ALPHA = alpha;
  MU = mu;
  DEL = del;
  RT = rt;

  /* actuator locations/magnitudes */
  // TODO: decide on how to place the actuators
  Aloc = malloc(2 * M * sizeof(double));
  Amag = malloc(M * sizeof(double));
  int mx = round(sqrt((double)M) * LX / sqrt(LX * LY));
  int my = M / mx;
  if (mx * my != M) {
    fprintf(stderr, "ERROR: M must be divisable in proportion to Lx and Ly\n");
    exit(1);
  }
  for (int i = 0; i < mx; i++) {
    for (int j = 0; j < my; j++) {
      int k = i * my + j;
      Aloc[2 * k] = (i + 0.5) * LX / mx;
      Aloc[2 * k + 1] = (j + 0.5) * LY / my;
      Amag[k] = 0.0;
    } // j end
  } // i end

  /* observer locations */
  // TODO: decide on how to place the observers
  Oloc = malloc(2 * P * sizeof(double));
  int px = round(sqrt((double)P) * LX / sqrt(LX * LY));
  int py = P / px;
  if (px * py != P) {
    fprintf(stderr, "ERROR: P must be divisable in proportion to Lx and Ly\n");
    exit(1);
  }
  for (int i = 0; i < px; i++) {
    for (int j = 0; j < py; j++) {
      int k = i * py + j;
      Oloc[2 * k] = (i + 0.5) * LX / px - DEL;
      Oloc[2 * k + 1] = (j + 0.5) * LY / py;

      /* wrap via periodicity */
      // TODO: this will break if |DEL| > LX
      if (Oloc[2 * k] < 0) {
        Oloc[2 * k] += LX;
      } else if (Oloc[2 * k] > LX) {
        Oloc[2 * k] -= LX;
      }
    }
  }

  /* control normaliser */
  NORM = 1.0;
  double integral = 0.0;
  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      integral += actuator(JTOX(j), ITOY(i));
    }
  } // i end
  NORM = 1.0 / (DX * DY * integral);
}

/* frees the common memory */
void internal_control_free(void) {
  free(Aloc);
  free(Oloc);
  free(Amag);
}

/* outputs common arrays */
void internal_control_output(void) {
  output_d1d("output/Aloc.dat", Aloc, 2 * M);
  output_d1d("output/Oloc.dat", Oloc, 2 * P);

  double **F = malloc_f2d(NX * NY, M);
  forcing_matrix(F);
  output_d2d("output/F.dat", F, NX * NY, M);

  double **A, **B, **C;

  // /* BENNEY */
  // /* Jacobian */
  // A = malloc_f2d(NX * NY, NX * NY);
  // benney_jacobian(A);
  //
  // /* actuator matrix */
  // B = malloc_f2d(NX * NY, M);
  // benney_actuator(B);
  //
  // /* observer matrix (actually the transpose) */
  // C = malloc_f2d(NX * NY, P);
  // benney_observer(C);
  //
  // output_d2d("output/A_be.dat", A, NX * NY, NX * NY);
  // output_d2d("output/B_be.dat", B, NX * NY, M);
  // output_d2d("output/C_be.dat", C, NX * NY, P);
  //
  // free_2d(A);
  // free_2d(B);
  // free_2d(C);

  /* WR */
  /* Jacobian */
  A = malloc_f2d(2 * NX * NY, 2 * NX * NY);
  wr_jacobian(A);

  /* actuator matrix */
  B = malloc_f2d(2 * NX * NY, M);
  wr_actuator(B);

  /* observer matrix (actually the transpose) */
  C = malloc_f2d(2 * NX * NY, P);
  wr_observer(C);

  output_d2d("output/A_wr.dat", A, 2 * NX * NY, 2 * NX * NY);
  output_d2d("output/B_wr.dat", B, 2 * NX * NY, M);
  output_d2d("output/C_wr.dat", C, 2 * NX * NY, P);

  free_2d(F);
  free_2d(A);
  free_2d(B);
  free_2d(C);
}

/* ======================= */
/*  ROM Matrix Generators  */
/* ======================= */
/* forcing matrix (N-by-M) */
void forcing_matrix(double **F) {
  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      for (int k = 0; k < M; k++) {
        F[IJTOK(i, j)][k] =
            actuator(JTOX(j) - Aloc[2 * k], ITOY(i) - Aloc[2 * k + 1]);
      }
    } // j end
  } // i end
}

/* Jacobian (N-by-N) */
void benney_jacobian(double **A) {
  fprintf(stderr, "ERROR: Benney jacobian not implemented in 3D\n");
  exit(1);
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     J[i][j] = 0.0;
  //   } // j end
  // } // i end
  //
  // double c0 = -1.0/(3.0*CA) * (1.0/(DX*DX*DX*DX));
  // double c1 = 1.0/DX + (2.0/tan(THETA)/3.0 - 8.0*RE/15.0) * (1.0/(DX*DX))
  // + 1.0/(3.0*CA) * (4.0/(DX*DX*DX*DX)); double c2 = (2.0/tan(THETA)/3.0
  // - 8.0*RE/15.0) * (-2.0/(DX*DX)) - 1.0/(3.0*CA) * (6.0/(DX*DX*DX*DX));
  // double c3 = -1.0/DX + (2.0/tan(THETA)/3.0 - 8.0*RE/15.0) * (1.0/(DX*DX))
  // + 1.0/(3.0*CA) * (4.0/(DX*DX*DX*DX)); double c4 = -1.0/(3.0*CA) *
  // (1.0/(DX*DX*DX*DX)); for (int i = 2; i < N-2; i++) {
  //   J[i][i-2] = c0;
  //   J[i][i-1] = c1;
  //   J[i][i] = c2;
  //   J[i][i+1] = c3;
  //   J[i][i+2] = c4;
  // } // i end
  //
  // J[0][N-2] = c0;
  // J[0][N-1] = c1;
  // J[0][0] = c2;
  // J[0][1] = c3;
  // J[0][2] = c4;
  //
  // J[1][N-1] = c0;
  // J[1][0] = c1;
  // J[1][1] = c2;
  // J[1][1+1] = c3;
  // J[1][1+2] = c4;
  //
  // J[N-2][N-4] = c0;
  // J[N-2][N-3] = c1;
  // J[N-2][N-2] = c2;
  // J[N-2][N-1] = c3;
  // J[N-2][0] = c4;
  //
  // J[N-1][N-3] = c0;
  // J[N-1][N-2] = c1;
  // J[N-1][N-1] = c2;
  // J[N-1][0] = c3;
  // J[N-1][1] = c4;
}

/* Actuator (N-by-M) */
void benney_actuator(double **B) {
  fprintf(stderr, "ERROR: Benney actuator not implemented in 3D\n");
  exit(1);
  // /* forcing matrix */
  // double **F = malloc_f2d(N, M);
  // forcing_matrix(F);
  //
  // /* actuator matrix */
  // for (int i = 1; i < N-1; i++) {
  //   for (int j = 0; j < M; j++) {
  //     Psi[i][j] = F[i][j] + (RE/(3*DX)) * (F[i+1][j]-F[i-1][j]);
  //   } // j end
  // } // i end
  // for (int j = 0; j < M; j++) {
  //   Psi[0][j] = F[0][j] + (RE/(3*DX)) * (F[1][j]-F[N-1][j]);
  //   Psi[N-1][j] = F[N-1][j] + (RE/(3*DX)) * (F[0][j]-F[N-2][j]);
  // } // j end
  //
  // free_2d(F);
}

/* (the transpose of the) Observer (N-by-P) */
void benney_observer(double **C) {
  fprintf(stderr, "ERROR: Benney observer not implemented in 3D\n");
  exit(1);
  // /* P < N use approximate delta-functions */
  // if (P != N) {
  //   for (int i = 0; i < N; i++) {
  //     for (int j = 0; j < P; j++) {
  //       Phi[i][j] = DX*actuator(ITOX(i)-Oloc[j]);
  //     } // j end
  //   } // i end
  // }
  //
  // /* if P == N then just pass all the information */
  // else {
  //   for (int i = 0; i < N; i++) {
  //     for (int j = 0; j < N; j++) {
  //       Phi[i][j] = 0.0;
  //     } // j end
  //     Phi[i][i] = 1.0;
  //   } // i end
  // }
}

/* Jacobian (2N-by-2N) */
void wr_jacobian(double **A) {
  for (int i = 0; i < 2 * NX * NY; i++) {
    for (int j = 0; j < 2 * NX * NY; j++) {
      A[i][j] = 0.0;
    } // j end
  } // i end

  // coefficients
  double c00, c01, c02, c03, c04;
  double c10, c11, c12, c13, c14;
  double c20, c21, c22, c23, c24;
  double c30, c31, c32, c33, c34;
  double c40, c41, c42, c43, c44;

  double BETA = 1.0 / tan(THETA);

  /* ============== TOP LEFT (hh) ============== */
  c02 = 0.0;
  c11 = 0.0;
  c12 = 0.0;
  c13 = 0.0;
  c20 = 0.0;
  c21 = 0.0;
  c22 = 0.0;
  c23 = 0.0;
  c24 = 0.0;
  c31 = 0.0;
  c32 = 0.0;
  c33 = 0.0;
  c42 = 0.0;

  // 2/3 cot(theta) Dyy
  c12 += (2.0 / 3.0 * BETA) * (1.0 / (DY * DY));
  c22 += (2.0 / 3.0 * BETA) * (-2.0 / (DY * DY));
  c32 += (2.0 / 3.0 * BETA) * (1.0 / (DY * DY));

  // - 1/3/Ca Dxxyy
  c11 += -1.0 / (3.0 * CA) * (1.0 / (DX * DX * DY * DY));
  c21 += -1.0 / (3.0 * CA) * (-2.0 / (DX * DX * DY * DY));
  c31 += -1.0 / (3.0 * CA) * (1.0 / (DX * DX * DY * DY));

  c12 += -1.0 / (3.0 * CA) * (-2.0 / (DX * DX * DY * DY));
  c22 += -1.0 / (3.0 * CA) * (4.0 / (DX * DX * DY * DY));
  c32 += -1.0 / (3.0 * CA) * (-2.0 / (DX * DX * DY * DY));

  c13 += -1.0 / (3.0 * CA) * (1.0 / (DX * DX * DY * DY));
  c23 += -1.0 / (3.0 * CA) * (-2.0 / (DX * DX * DY * DY));
  c33 += -1.0 / (3.0 * CA) * (1.0 / (DX * DX * DY * DY));

  // - 1/3/Ca Dyyyy
  c02 += -1.0 / (3.0 * CA) * (1.0 / (DY * DY * DY * DY));
  c12 += -1.0 / (3.0 * CA) * (-4.0 / (DY * DY * DY * DY));
  c22 += -1.0 / (3.0 * CA) * (6.0 / (DY * DY * DY * DY));
  c32 += -1.0 / (3.0 * CA) * (-4.0 / (DY * DY * DY * DY));
  c42 += -1.0 / (3.0 * CA) * (1.0 / (DY * DY * DY * DY));

  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      A[IJTOK(i, j)][IJTOK(i - 2, j + 0)] = c02;
      A[IJTOK(i, j)][IJTOK(i - 1, j - 1)] = c11;
      A[IJTOK(i, j)][IJTOK(i - 1, j + 0)] = c12;
      A[IJTOK(i, j)][IJTOK(i - 1, j + 1)] = c13;
      A[IJTOK(i, j)][IJTOK(i + 0, j - 2)] = c20;
      A[IJTOK(i, j)][IJTOK(i + 0, j - 1)] = c21;
      A[IJTOK(i, j)][IJTOK(i + 0, j + 0)] = c22;
      A[IJTOK(i, j)][IJTOK(i + 0, j + 1)] = c23;
      A[IJTOK(i, j)][IJTOK(i + 0, j + 2)] = c24;
      A[IJTOK(i, j)][IJTOK(i + 1, j - 1)] = c31;
      A[IJTOK(i, j)][IJTOK(i + 1, j + 0)] = c32;
      A[IJTOK(i, j)][IJTOK(i + 1, j + 1)] = c33;
      A[IJTOK(i, j)][IJTOK(i + 2, j + 0)] = c42;
    } // j end
  } // i end

  /* ============== TOP RIGHT (hq) ============== */
  c02 = 0.0;
  c11 = 0.0;
  c12 = 0.0;
  c13 = 0.0;
  c20 = 0.0;
  c21 = 0.0;
  c22 = 0.0;
  c23 = 0.0;
  c24 = 0.0;
  c31 = 0.0;
  c32 = 0.0;
  c33 = 0.0;
  c42 = 0.0;

  // -Dx
  c21 += -0.5 / DX;
  c23 += 0.5 / DX;

  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      A[IJTOK(i, j)][NX * NY + IJTOK(i - 2, j + 0)] = c02;
      A[IJTOK(i, j)][NX * NY + IJTOK(i - 1, j - 1)] = c11;
      A[IJTOK(i, j)][NX * NY + IJTOK(i - 1, j + 0)] = c12;
      A[IJTOK(i, j)][NX * NY + IJTOK(i - 1, j + 1)] = c13;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 0, j - 2)] = c20;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 0, j - 1)] = c21;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 0, j + 0)] = c22;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 0, j + 1)] = c23;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 0, j + 2)] = c24;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 1, j - 1)] = c31;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 1, j + 0)] = c32;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 1, j + 1)] = c33;
      A[IJTOK(i, j)][NX * NY + IJTOK(i + 2, j + 0)] = c42;
    } // j end
  } // i end

  /* ============== BOTTOM LEFT (qh) ============== */
  c02 = 0.0;
  c11 = 0.0;
  c12 = 0.0;
  c13 = 0.0;
  c20 = 0.0;
  c21 = 0.0;
  c22 = 0.0;
  c23 = 0.0;
  c24 = 0.0;
  c31 = 0.0;
  c32 = 0.0;
  c33 = 0.0;
  c42 = 0.0;

  // 5/Re I
  c22 += 5.0 / RE;

  // + (4/7 - 5/3/RE cot(theta)) Dx
  c21 += (4.0 / 7.0 - 5.0 / 3.0 / RE * BETA) * (-0.5 / DX);
  c23 += (4.0 / 7.0 - 5.0 / 3.0 / RE * BETA) * (0.5 / DX);

  // + 5/6/Ca/Re Dxxx
  c20 += 5.0 / 6.0 / CA / RE * (-0.5 / (DX * DX * DX));
  c21 += 5.0 / 6.0 / CA / RE * (1.0 / (DX * DX * DX));
  c23 += 5.0 / 6.0 / CA / RE * (-1.0 / (DX * DX * DX));
  c24 += 5.0 / 6.0 / CA / RE * (0.5 / (DX * DX * DX));

  // + 5/6/Ca/Re Dxyy
  c11 += 5.0 / 6.0 / CA / RE * (-0.5 / (DX * DY * DY));
  c21 += 5.0 / 6.0 / CA / RE * (1.0 / (DX * DY * DY));
  c31 += 5.0 / 6.0 / CA / RE * (-0.5 / (DX * DY * DY));

  c13 += 5.0 / 6.0 / CA / RE * (0.5 / (DX * DY * DY));
  c23 += 5.0 / 6.0 / CA / RE * (-1.0 / (DX * DY * DY));
  c33 += 5.0 / 6.0 / CA / RE * (0.5 / (DX * DY * DY));

  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      A[NX * NY + IJTOK(i, j)][IJTOK(i - 2, j + 0)] = c02;
      A[NX * NY + IJTOK(i, j)][IJTOK(i - 1, j - 1)] = c11;
      A[NX * NY + IJTOK(i, j)][IJTOK(i - 1, j + 0)] = c12;
      A[NX * NY + IJTOK(i, j)][IJTOK(i - 1, j + 1)] = c13;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 0, j - 2)] = c20;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 0, j - 1)] = c21;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 0, j + 0)] = c22;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 0, j + 1)] = c23;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 0, j + 2)] = c24;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 1, j - 1)] = c31;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 1, j + 0)] = c32;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 1, j + 1)] = c33;
      A[NX * NY + IJTOK(i, j)][IJTOK(i + 2, j + 0)] = c42;
    } // j end
  } // i end

  /* ============== BOTTOM RIGHT (qq) ============== */
  c02 = 0.0;
  c11 = 0.0;
  c12 = 0.0;
  c13 = 0.0;
  c20 = 0.0;
  c21 = 0.0;
  c22 = 0.0;
  c23 = 0.0;
  c24 = 0.0;
  c31 = 0.0;
  c32 = 0.0;
  c33 = 0.0;
  c42 = 0.0;

  // -5/2/Re I
  c22 += -5.0 / 2.0 / RE;

  // -34/21 Dx
  c21 += -34.0 / 21.0 * (-0.5 / DX);
  c23 += -34.0 / 21.0 * (0.5 / DX);

  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i - 2, j + 0)] = c02;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i - 1, j - 1)] = c11;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i - 1, j + 0)] = c12;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i - 1, j + 1)] = c13;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 0, j - 2)] = c20;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 0, j - 1)] = c21;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 0, j + 0)] = c22;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 0, j + 1)] = c23;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 0, j + 2)] = c24;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 1, j - 1)] = c31;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 1, j + 0)] = c32;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 1, j + 1)] = c33;
      A[NX * NY + IJTOK(i, j)][NX * NY + IJTOK(i + 2, j + 0)] = c42;
    } // j end
  } // i end
}

/* Actuator (2N-by-M) */
void wr_actuator(double **B) {
  /* forcing matrix (TODO: see if rotating this makes it faster) */
  double **F = malloc_f2d(NX * NY, M);
  forcing_matrix(F);

  /* actuator matrix */
  for (int i = 0; i < NX * NY; i++) {
    for (int j = 0; j < M; j++) {
      B[i][j] = F[i][j];
      B[i + NX * NY][j] = (1.0 / 3.0) * F[i][j];
    } // j end
  } // i end

  free_2d(F);
}

/* (the transpose of the) Observer (2N-by-P) */
void wr_observer(double **C) {
  fprintf(stderr, "ERROR: WR observer not implemented in 3D\n");
  exit(1);

  // /* P < N use approximate delta-functions */
  // if (P != N) {
  //   for (int i = 0; i < 2*N; i++) {
  //     for (int j = 0; j < P; j++) {
  //       Phi[i][j] = 0.0;
  //     } // j end
  //   } // i end
  //
  //   /* observe interfacial height */
  //   for (int i = 0; i < N; i++) {
  //     for (int j = 0; j < P; j++) {
  //       Phi[i][j] = DX*actuator(ITOX(i)-Oloc[j]);
  //     } // j end
  //   } // i end
  //
  //   /* observe flux */ // TODO: should the q = 2/3 h assumption be here?
  //   // fprintf(stderr, "WARNING: this should use the flux approximation\n");
  //   // for (int i = 0; i < N; i++) {
  //   //   for (int j = 0; j < P; j++) {
  //   //     Phi[N+i][P+j] = DX*actuator(ITOX(i)-Oloc[j]);
  //   //   } // j end
  //   // } // i end
  // }
  //
  // /* if P == N then just pass all the information */
  // else {
  //   for (int i = 0; i < 2*N; i++) {
  //     for (int j = 0; j < P; j++) {
  //       Phi[i][j] = 0.0;
  //     } // j end
  //     Phi[i][i] = 1.0;
  //   } // i end
  // }
}

#ifdef __cplusplus
}
#endif

#endif
