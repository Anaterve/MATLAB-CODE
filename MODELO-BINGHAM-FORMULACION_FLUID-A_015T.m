#include "udf.h"

DEFINE_PROPERTY(Bingham_Viscoplastico_FluidA_0T, c, t)
{
    real mu_eff;
    real mu_p = 0.909;     // Pa·s
    real gamma_dot = C_STRAIN_RATE_MAG(c, t);
    real vartheta;
    real expo;
	
    // Modelo exponencial para tau_y(H)
    vartheta = 1857.9;   //Pa

    // Seguridad numérica
    if (gamma_dot < 1e-12) 
        gamma_dot = 1e-12;

    // Cálculo de viscosidad efectiva
    mu_eff = mu_p + vartheta / gamma_dot;

    // Opcional: límite de seguridad máxima
    if (mu_eff > 1e6) 
        mu_eff = 1e6;

    return mu_eff;
}

