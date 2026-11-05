/******************************************************************************
 *                                                                            *
 * glm_ss_density.h                                                           *
 *                                                                            *
 * Contribution of suspended solids (SS) to water density                     *
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
#ifndef _GLM_SS_DENSITY_H_
#define _GLM_SS_DENSITY_H_

#include "glm.h"

#define MAX_SS_DENS_VARS  20

extern CLOGICAL ss_density;       //# include SS in inflow and layer density (&wq_setup)
extern CLOGICAL ss_density_diag;  //# write per-parcel inflow diagnostics (&wq_setup)

void ss_density_init(const char *wq_nml_file, char **vars, int nvars,
                                               AED_REAL *rho_override, int nrho);

//# Density excess [kg/m3] of a set of WQ values (e.g. an inflow parcel) in
//# water of density rho_w. Returns zero when ss_density is off.
AED_REAL ss_density_excess(const AED_REAL *wq, AED_REAL rho_w);

//# Density excess [kg/m3] of the SS held in lake layer k
AED_REAL ss_layer_excess(int k, AED_REAL rho_w);

//# Density of lake layer k from T, S and (when on) its SS
AED_REAL layer_density(int k);

//# Recalculate the density of every lake layer
void refresh_densities(void);

#endif
