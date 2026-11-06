#include "udf.h"

/* ================================================================
   UDF de viscosidad Bingham para una transicion temporal de Vartheta

   Tramo 1:       0 <= t < Tt1       -> Vartheta_0
   Rampa:         Tt1 <= t < Tt1+Tr  -> interpolacion lineal
   Tramo 2:       t >= Tt1+Tr         -> Vartheta_f

   Para obtener un escalon, usar Tr = 0.0 s.
   ================================================================ */

/* Parametros facilmente modificables (unidades SI) */
#define Tt1          1.0       /* duracion inicial [s] */
#define Tr           0.0       /* duracion de rampa [s]; 0.0 = escalon */
#define Tt2          1.0       /* duracion final [s] (informativa) */

#define VART_THETA_0 0.0       /* umbral inicial [Pa], estado 0T */
#define VART_THETA_F 680.6    /* umbral final [Pa], estado 0.06T */

#define MU_PLASTICA  0.909     /* viscosidad plastica [Pa s] */
#define GAMMA_MIN    1.0e-12  /* limite numerico [1/s] */
#define MU_MAX       1.0e6     /* limite numerico [Pa s] */

static real vartheta_transitorio(real time)
{
    real vartheta;

    if (time < Tt1)
    {
        vartheta = VART_THETA_0;
    }
    else if (Tr <= 0.0)
    {
        /* Escalon instantaneo en t = Tt1 */
        vartheta = VART_THETA_F;
    }
    else if (time < Tt1 + Tr)
    {
        vartheta = VART_THETA_0
                 + (VART_THETA_F - VART_THETA_0)
                 * (time - Tt1) / Tr;
    }
    else
    {
        vartheta = VART_THETA_F;
    }

    return vartheta;
}

DEFINE_PROPERTY(Bingham_Viscoplastico_FluidA_Transitorio, c, t)
{
    real gamma_dot = C_STRAIN_RATE_MAG(c, t);
    real vartheta = vartheta_transitorio(CURRENT_TIME);
    real mu_eff;

    if (gamma_dot < GAMMA_MIN)
        gamma_dot = GAMMA_MIN;

    mu_eff = MU_PLASTICA + vartheta / gamma_dot;

    if (mu_eff > MU_MAX)
        mu_eff = MU_MAX;

    return mu_eff;
}

