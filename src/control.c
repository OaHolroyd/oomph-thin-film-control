#ifdef __cplusplus
extern "C" {
#endif

#include "control.h"

#include "control-core.h"
#include "control-dynamic.h"
#include "control-lqr.h"
#include "control-pair.h"
#include "control-static.h"

#include "c-utils.h"

static control_t CT;

/* ========================================================================== */
/*   STRATEGY-SPECIFIC FUNCTION DECLARATIONS                                  */
/* ========================================================================== */
/* sets up the specific control strategy */
void (*s_set)(void);

/* frees varaibles specific to the selected strategy */
void (*s_free)(void);

/* outputs matrices specific to the selected strategy */
void (*s_output)(void);

/* sets CM to be the control matrix (so that f = CM*(h-1)) */
void (*control_matrix)(double **CM);

/* steps the specific control system forward in time given the interfacial
   height */
void (*control_step)(double dt, double *h, double *q);

/* returns the estimator as a function of x */
double (*estimator)(double x, double y);

/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* returns the baseplate control velocity as a function of x, y */
double control(double x, double y) {
  double vc = 0.0;

  for (int i = 0; i < M; i++) {
    vc += Amag[i] * actuator(x - Aloc[2 * i], y - Aloc[2 * i + 1]);
  } // i end

  return -ALPHA * vc;
}

/* computes the incremental cost from the current timestep from the interface */
double control_cost(double *h) {
  double cost = 0.0;

  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      double x = JTOX(j);
      double y = ITOY(i);

      /* control cost */
      double a = control(x, y);
      cost += (1 - MU) * a * a;

      /* interfacial cost */
      double dh = h[IJTOK(i, j)] - 1.0;
      cost += MU * dh * dh;
    } // j end
  }

  /* integral scaling */
  cost *= DX * DY;

  return cost;
}

/* set up the control system */
// TODO: explain parameters properly
void control_set(control_t ct, rom_t rt, int m, int p, double w, double alpha,
                 double mu, double del, double lx, double ly, int nx, int ny,
                 double re, double ca, double theta) {
  fprintf(stderr, "ENTERED control_set\n");

  /* control strategy independent setup */
  internal_control_set(rt, m, p, w, alpha, mu, del, lx, ly, nx, ny, re, ca,
                       theta);

  return;

  /* set strategy specific functions */
  CT = ct;
  switch (CT) {
  case PAIR:
    s_set = &pair_set;
    s_free = &pair_free;
    control_step = &pair_step;
    estimator = &pair_estimator;
    s_output = &pair_output;
    control_matrix = &pair_matrix;
    break;
  case LQR:
    s_set = &lqr_set;
    s_free = &lqr_free;
    control_step = &lqr_step;
    estimator = &lqr_estimator;
    s_output = &lqr_output;
    control_matrix = &lqr_matrix;
    break;
  case STATIC:
    s_set = &static_set;
    s_free = &static_free;
    control_step = &static_step;
    estimator = &static_estimator;
    s_output = &static_output;
    control_matrix = &static_matrix;
    break;
  case DYNAMIC:
    s_set = &dynamic_set;
    s_free = &dynamic_free;
    control_step = &dynamic_step;
    estimator = &dynamic_estimator;
    s_output = &dynamic_output;
    control_matrix = NULL;
    break;
  default:
    ABORT("invalid control type %d", ct);
  }

  /* call strategy specific setup */
  s_set();
}

/* frees the control system */
void control_free(void) {
  internal_control_free();
  s_free();
}

/* outputs the relevant control matrices */
void control_output(void) {
  internal_control_output();
  s_output();
}

#ifdef __cplusplus
}
#endif
