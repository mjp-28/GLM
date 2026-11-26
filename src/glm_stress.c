/******************************************************************************
 *                                                                            *
 * glm_stress.c                                                               *
 *                                                                            *
 *   Calculate layer stress                                                   *
 *                                                                            *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013-2026 : The University of Western Australia                  *
 *                                                                            *
 *  This file is part of GLM (General Lake Model)                             *
 *                                                                            *
 *  GLM is free software: you can redistribute it and/or modify               *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  GLM is distributed in the hope that it will be useful,                    *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "glm.h"
#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"
#include "glm_stress.h"

#define d_50  80e-6
#define Ks    ( 2.5*d_50 )

AED_REAL Hs = 0., T = 0., L = 0.;
static AED_REAL BottomStress = 0.;
static AED_REAL U_last = 0., F_last = 0.;   //# zone_area_stress: wind and fetch of the last calc_layer_stress

#define DEBUG_STRESS  0
#define DEBUG_STRESS_CSV  0

#if DEBUG_STRESS
#undef DEBUG_STRESS_CSV
#define DEBUG_STRESS_CSV 1
#endif
#if DEBUG_STRESS_CSV
static int seenit = FALSE;
static int stepno = 0;
static FILE *csv_dbg = NULL;
static int new_run = TRUE;
#endif

// Variant 1 is as per CAEDYM, variant 0 as per Ji (2008) and Laenen and LeTourneau (1996)
#define VARIANT 0

/******************************************************************************/
#define f_c(h) ( 0.24 / pow(log10(12*h / Ks), 2) )
/******************************************************************************/
#define f_Uorb(F, h) ((Pi * Hs) / (T * sinh((two_Pi * h) / L)))

/******************************************************************************/
static AED_REAL f_L(AED_REAL T, AED_REAL h)
{
    AED_REAL L0 = ((g * T * T) / (two_Pi));

    return L0 * pow( tanh( (two_Pi * h) / L0 ), 0.5);
}

/******************************************************************************/
static AED_REAL f_Hs(AED_REAL U_Sqr, AED_REAL F, AED_REAL h)
{
    AED_REAL Zeta = (0.53 * pow((g * h) / U_Sqr, 0.75));

#if VARIANT
    return 0.283 * (U_Sqr / g) * tanh(Zeta) * tanh(0.00565 * pow((g * F) / U_Sqr, 0.5) / tanh(Zeta));
#else
    return 0.283 * (U_Sqr / g) * tanh(Zeta) * tanh(0.0125 * pow((g * F) / U_Sqr, 0.42) / tanh(Zeta));
#endif
}

/******************************************************************************/
static AED_REAL f_T(AED_REAL U, AED_REAL U_Sqr, AED_REAL F, AED_REAL h)
{
    AED_REAL Xi = (0.833 * pow((g * h) / U_Sqr, 0.375));

#if VARIANT
    return 1.2 * two_Pi * (U / g) * tanh(Xi) * tanh(0.0379 * pow((g * F) / U_Sqr, 0.333) / tanh(Xi));
#else
    return 1.2 * two_Pi * (U / g) * tanh(Xi) * tanh(0.077 * pow((g * F) / U_Sqr, 0.25) / tanh(Xi));
#endif
}


/******************************************************************************
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
static AED_REAL wave_friction_factor(AED_REAL Uorb, AED_REAL dens)
{
    // Swart (1974) and Kleinhans & Grasmeijer (2006) Eq.9
    // also cited in http://nora.nerc.ac.uk/8360/1/POL_ID_189.pdf eqn 2.3.8

    return MIN( exp(-5.977 + 5.213 * pow((MAX(Uorb, 0.0001) * MAX(T, 1) / (2 * two_Pi * Ks)), -0.194)), 0.1);

    // and Le Roux (2001) Eq 15. as:
    //return MIN(0.00251 * exp(5.213 * pow((MAX(Uorb, 0.0001) * MAX(T, 1) / (2 * two_Pi * Ks)), -0.19)), 0.1);

    //# Le Roux (2001) ??
    //AED_REAL beta = zero
    //return MIN( (2.0*beta*9.81*2650.0*d_50)/( MAX(Uorb,0.0001)*MAX(Uorb,0.0001)*dens ), 0.1);
}

/******************************************************************************
 * calc_layer_stress                                                          *
 *                                                                            *
 *   U : Wind velocity 10m above the surface                                  *
 *   F : Fetch                                                                *
 ******************************************************************************/
void calc_layer_stress(AED_REAL U, AED_REAL F)
{
    int i;
    AED_REAL Uorb, Ucur, h, dens;
    AED_REAL U_Sqr;


    if (ice) U = 0.00001;
    U_last = U; F_last = F;

#if DEBUG_STRESS_CSV
    if ( new_run ) {
        csv_dbg = fopen("stress_dbg.csv", "w");
        new_run = FALSE;
        fprintf(csv_dbg, "stepno, U, F, Uorb(top), Uorb(bot), Tuab(top), Taub(bot)\n");
    }
#endif

    if ( U > 0.001 ) {
        U_Sqr = U * U; // Module global used by various subroutines above

        h  = Lake[surfLayer].Height / 2.;
        T  = f_T(U, U_Sqr, F, h);
        Hs = f_Hs(U_Sqr, F, h);
        L  = f_L(T, h);

        Lake[surfLayer].Umean = sqrt( (1.2/Lake[surfLayer].Density) * coef_wind_drag * U * U);

#if DEBUG_STRESS
        if ( !seenit ) {
            fprintf(stderr, "calc_layer_stress U = %f ; F  = %f ; Umean %f\n", U, F, Lake[surfLayer].Umean);
            fprintf(stderr, "calc_layer_stress T = %f ; Hs = %f ; L = %f\n", T, Hs, L);
        }
#endif

        for (i = surfLayer; i >= botmLayer; i--) {
            if (i == surfLayer) Ucur =  Lake[i].Umean;
            else                Ucur = (Lake[i].Umean = 0.);

            if ( i != botmLayer ) h = Lake[surfLayer].Height - Lake[i-1].Height;
            else                  h = Lake[surfLayer].Height;

            Uorb = MIN( f_Uorb(F, h) , 5.0);
            Lake[i].Uorb = Uorb;

            dens = Lake[i].Density;
            Lake[i].LayerStress = dens *
                     ((0.5 * wave_friction_factor(Uorb, dens) * pow(Uorb, 2)) +
                                                   (f_c(h) * pow(Ucur, 2) / 8));

#if DEBUG_STRESS
            if ( !seenit )
                fprintf(stderr, "Layer %02d : h = %f Uorb = %e Taub = %e U = %e F = %e\n",
                                                           i, h, Lake[i].Uorb, Lake[i].LayerStress, U, F);

            if (Lake[i].LayerStress < 0. || Lake[i].LayerStress > 1.00e4) {
                fprintf(stderr, "L(S-%3d) Taub out of range %e ; U = %e ; F = %e ; h %e Dens %e ; ice %s\n",
                                    surfLayer-i, Lake[i].LayerStress, U, F, h, dens, (ice)?"true":"false");
//              exit(0);
            } else {
                if ( h < 1. && i < surfLayer-3) {
                    fprintf(stderr, "step  %6d L(S-%3d) OK h %e\n", stepno, surfLayer-i, h);
//                  exit(0);
                }
            }
#endif
        }
    } else { // U is 0 thus all else should be too
        for (i = surfLayer; i >= botmLayer; i--) {
            Lake[i].Umean = 0.;
            Lake[i].LayerStress = 0.;
        }
    }
    BottomStress = Lake[botmLayer].LayerStress;
    if ( U <= 0.001 ) Hs = 0.;   //# no waves: zone_area_stress returns 0

#if DEBUG_STRESS_CSV
    fprintf(csv_dbg, "%d, %e, %e, %e, %e, %e, %e\n",
                stepno, U, F, Lake[surfLayer].Uorb, Lake[botmLayer].Uorb,
                              Lake[surfLayer].LayerStress, Lake[botmLayer].LayerStress);

//  fprintf(stderr, "step %6d L(S-%3d) Taub at the bottom %e ; U = %e ; F = %e ; h %e ice %s\n",
//                                      stepno, surfLayer, Lake[botmLayer].LayerStress, U, F, h, (ice)?"true":"false");

    seenit = TRUE;
    stepno++;
#endif
}


/******************************************************************************
 * area_at_height  (zone_taub_area helper)                                    *
 *   Plan area of the lake at height h above the bottom, from the same        *
 *   morphometry arrays GLM uses for layer areas (glm_mixu.c resize_internals) *
 ******************************************************************************/
static AED_REAL area_at_height(AED_REAL h)
{
    AED_REAL x = h * MphInc, y = AMOD(x, 1.0);
    int ij = (x - y);
    if (ij > Nmorph) { y += (ij - Nmorph); ij = Nmorph; }
    if (ij > 0) { ij--; return MphLevelArea[ij] + y * dMphLevelArea[ij]; }
    return MphLevelArea[0];
}

/******************************************************************************
 * zone_area_stress  (&wq_setup zone_taub_area; called from glm_zones.F90)    *
 *   Wave bed stress averaged over the bed area of a sediment zone (heights   *
 *   zb..zt above the bottom, cut at the water surface), in 5 cm height bands. *
 *   Each band uses the wave orbital velocity at its own depth (same Hs, T, L  *
 *   and friction factor as calc_layer_stress), with depth floored at          *
 *   max(zone_taub_hmin, 1.28 Hs) for the swash / breaking fringe.             *
 *   Returns -1 when zone_taub_area is off (the layer-mean stress is kept).    *
 ******************************************************************************/
AED_REAL zone_area_stress(AED_REAL zb, AED_REAL zt)
{
    AED_REAL surf, top, h, h2, dA, depth, dmin, Uorb, tau, sum = 0., wsum = 0., dens;
    const AED_REAL dh = 0.05;

    if ( !zone_taub_area ) return -1.;
    if ( ice || U_last <= 0.001 || Hs <= 0. ) return 0.;
    surf = Lake[surfLayer].Height;
    top = MIN(zt, surf);
    if ( top <= zb ) return 0.;
    dens = Lake[surfLayer].Density;
    dmin = MAX(zone_taub_hmin, 1.28 * Hs);
    for (h = zb; h < top; h += dh) {
        h2 = MIN(h + dh, top);
        dA = area_at_height(h2) - area_at_height(h);
        if ( dA <= 0. ) continue;
        depth = surf - 0.5 * (h + h2);
        if ( depth < dmin ) depth = dmin;
        Uorb = MIN( f_Uorb(F_last, depth), 5.0 );
        tau  = dens * 0.5 * wave_friction_factor(Uorb, dens) * pow(Uorb, 2);
        sum += tau * dA; wsum += dA;
    }
    return (wsum > 0.) ? sum / wsum : 0.;
}
