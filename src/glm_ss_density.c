/******************************************************************************
 *                                                                            *
 * glm_ss_density.c                                                           *
 *                                                                            *
 * Contribution of suspended solids (SS) to water density.                    *
 *                                                                            *
 * GLM normally derives density from temperature and salinity only. With      *
 * ss_density = .true. in &wq_setup, the mass of the listed particulate WQ    *
 * variables (by default the AED non-cohesive groups NCS_ss1..n, in g/m3)     *
 * is added to inflow and layer density:                                      *
 *                                                                            *
 *     rho = rho(T,S) + sum_i C_i * 1e-3 * (1 - rho(T,S)/rho_s,i)             *
 *                                                                            *
 * where C_i is the dry-mass concentration [g/m3] and rho_s,i the particle    *
 * (grain) density [kg/m3]. rho_s is read from rho_ss in the                  *
 * &aed_noncohesive group of the AED namelist so that GLM and AED share one   *
 * value; ss_density_rho in &wq_setup overrides it.                           *
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
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_util.h"
#include "glm_wqual.h"
#include "glm_ss_density.h"
#include "aed_time.h"

#define NCS_GROUP    "aed_noncohesive"
#define NCS_MAX      100      //# AED dimensions its NCS parameter vectors at 100
#define NML_LINE_LEN 8192

/*----------------------------------------------------------------------------*/
CLOGICAL ss_density      = FALSE;
CLOGICAL ss_density_diag = FALSE;

static int      n_ss = 0;                          //# number of active SS vars
static int      ss_idx[MAX_SS_DENS_VARS];          //# their WQ indices
static AED_REAL ss_rho[MAX_SS_DENS_VARS];          //# their particle density

static char     diag_path[512] = "";               //# inflow diagnostics file
static FILE    *diag_f    = NULL;
static int      diag_jday = 0;
/*----------------------------------------------------------------------------*/


/******************************************************************************
 * Parse a Fortran namelist value list ("1.6e3, 2*2650, 1.5d-5") into vals,   *
 * starting at position *pos. Returns the number of values stored.            *
 ******************************************************************************/
static int parse_values(char *s, AED_REAL *vals, int *pos, int maxn)
{
    int nread = 0;
    char *tok = strtok(s, ", \t");

    while ( tok != NULL ) {
        int rep = 1;
        char *star = strchr(tok, '*'), *p;

        if ( star != NULL ) { rep = atoi(tok); tok = star + 1; }
        for (p = tok; *p; p++) if ( *p == 'd' || *p == 'D' ) *p = 'e';

        if ( *tok != 0 ) {
            AED_REAL v = strtod(tok, NULL);
            while ( rep-- > 0 && *pos < maxn ) { vals[(*pos)++] = v; nread++; }
        }
        tok = strtok(NULL, ", \t");
    }
    return nread;
}


/******************************************************************************
 * Read num_ss, settling, rho_ss, d_ss and w_ss from the &aed_noncohesive     *
 * group of the AED namelist. Written separately from the GLM namelist reader *
 * so that the rest of the AED file (stray '/', quotes in comments) cannot    *
 * stop the run. Returns 1 if the group was found, 0 if not, -1 if the file   *
 * could not be opened.                                                       *
 ******************************************************************************/
static int read_aed_ncs(const char *fname, int *num_ss, int *settling,
                          AED_REAL *rho_ss, AED_REAL *d_ss, AED_REAL *w_ss)
{
    FILE *f;
    char *line, *s, key[64] = "";
    int in_group = FALSE, found = FALSE, pos = 0, i;
    AED_REAL tmp[NCS_MAX];

    for (i = 0; i < NCS_MAX; i++) { rho_ss[i] = 1.6e3; d_ss[i] = 1e-6; w_ss[i] = 0.; }
    *num_ss = 0; *settling = 0;

    if ( (f = fopen(fname, "r")) == NULL ) return -1;
    line = malloc(NML_LINE_LEN);

    while ( fgets(line, NML_LINE_LEN, f) != NULL ) {
        char quote = 0, *eq, *vals;

        //# strip comments ('!' outside quotes) and trailing white space
        for (s = line; *s; s++) {
            if ( quote ) { if ( *s == quote ) quote = 0; }
            else if ( *s == '\'' || *s == '"' ) quote = *s;
            else if ( *s == '!' ) { *s = 0; break; }
        }
        s = line + strlen(line);
        while ( s > line && isspace((unsigned char)s[-1]) ) *--s = 0;
        s = line;
        while ( isspace((unsigned char)*s) ) s++;
        if ( *s == 0 ) continue;

        if ( !in_group ) {
            if ( *s == '&' && strncasecmp(s+1, NCS_GROUP, strlen(NCS_GROUP)) == 0 &&
                 (s[1+strlen(NCS_GROUP)] == 0 || isspace((unsigned char)s[1+strlen(NCS_GROUP)])) )
                in_group = found = TRUE;
            continue;
        }
        if ( *s == '/' || *s == '&' ) break;   //# end of group

        if ( (eq = strchr(s, '=')) != NULL ) {
            char *k = s, *e = eq, *paren;
            while ( e > k && isspace((unsigned char)e[-1]) ) e--;
            *e = 0;
            pos = 0;
            if ( (paren = strchr(k, '(')) != NULL ) {   //# e.g. rho_ss(2) = ...
                pos = atoi(paren+1) - 1;
                if ( pos < 0 ) pos = 0;
                *paren = 0;
            }
            strncpy(key, k, sizeof(key)-1); key[sizeof(key)-1] = 0;
            vals = eq + 1;
        } else
            vals = s;                           //# continuation of the previous key

        if      ( strcasecmp(key, "rho_ss") == 0 ) parse_values(vals, rho_ss, &pos, NCS_MAX);
        else if ( strcasecmp(key, "d_ss")   == 0 ) parse_values(vals, d_ss,   &pos, NCS_MAX);
        else if ( strcasecmp(key, "w_ss")   == 0 ) parse_values(vals, w_ss,   &pos, NCS_MAX);
        else if ( strcasecmp(key, "num_ss") == 0 || strcasecmp(key, "settling") == 0 ) {
            int p0 = 0;
            if ( parse_values(vals, tmp, &p0, 1) > 0 ) {
                if ( strcasecmp(key, "num_ss") == 0 ) *num_ss = (int)tmp[0];
                else                                  *settling = (int)tmp[0];
            }
        }
    }

    free(line);
    fclose(f);
    return found;
}


/******************************************************************************
 * Set up the SS density contribution. Called once from init_glm after the    *
 * WQ model has been initialised (so WQ variable indices are known).          *
 *   wq_nml_file  : AED namelist (for rho_ss etc)                             *
 *   vars, nvars  : ss_density_vars from &wq_setup (NULL/0 = all NCS groups)  *
 *   rho_override : ss_density_rho from &wq_setup (NULL/0 = use rho_ss)       *
 *   out_dir      : output directory for the ss_density_diag file             *
 ******************************************************************************/
void ss_density_init(const char *wq_nml_file, char **vars, int nvars,
                     AED_REAL *rho_override, int nrho, const char *out_dir)
{
    int num_ss, settling, found, i;
    AED_REAL rho_ss[NCS_MAX], d_ss[NCS_MAX], w_ss[NCS_MAX];
    char name[64];

    n_ss = 0;
    if ( !ss_density && !ss_density_diag ) return;

    if ( ss_density )
        fprintf(stderr, "     SS density: suspended solids included in inflow and layer density\n");
    else
        fprintf(stderr, "     SS density: OFF (density from T and S only); inflow diagnostics only\n");

    if ( !wq_calc || Num_WQ_Vars <= 0 ) {
        fprintf(stderr, "     WARNING: ss_density needs the WQ model (AED) to be running;"
                        " SS density switched OFF\n");
        ss_density = FALSE; ss_density_diag = FALSE;
        return;
    }

    found = read_aed_ncs(wq_nml_file, &num_ss, &settling, rho_ss, d_ss, w_ss);
    if ( found < 0 )
        fprintf(stderr, "     WARNING: could not open '%s' to read rho_ss\n", wq_nml_file);
    else if ( found == 0 )
        fprintf(stderr, "     WARNING: no &%s group in '%s'\n", NCS_GROUP, wq_nml_file);
    if ( num_ss > NCS_MAX ) num_ss = NCS_MAX;

    if ( nvars <= 0 ) nvars = num_ss;           //# default: every NCS group
    if ( nvars > MAX_SS_DENS_VARS ) {
        fprintf(stderr, "     WARNING: only the first %d ss_density_vars are used\n", MAX_SS_DENS_VARS);
        nvars = MAX_SS_DENS_VARS;
    }
    if ( nrho > 1 && nrho != nvars )
        fprintf(stderr, "     WARNING: %d ss_density_rho values for %d variables\n", nrho, nvars);

    for (i = 0; i < nvars; i++) {
        int idx, grp = -1;
        size_t len;
        AED_REAL rho = -1.;
        const char *src = NULL;

        if ( vars != NULL ) { strncpy(name, vars[i], sizeof(name)-1); name[sizeof(name)-1] = 0; }
        else snprintf(name, sizeof(name), "NCS_ss%d", i+1);

        len = strlen(name);
        if ( (idx = wq_var_index_c(name, &len)) < 0 ) {
            fprintf(stderr, "     WARNING: ss_density variable '%s' is not a WQ variable; ignored\n", name);
            continue;
        }

        if ( strncasecmp(name, "NCS_ss", 6) == 0 ) {
            grp = atoi(name+6) - 1;
            if ( grp < 0 || grp >= num_ss ) grp = -1;
        }

        if ( nrho > 0 ) {
            rho = rho_override[(nrho == 1) ? 0 : ((i < nrho) ? i : nrho-1)];
            src = "ss_density_rho";
            if ( grp >= 0 && fabs(rho - rho_ss[grp]) > 1e-6 )
                fprintf(stderr, "     WARNING: %s: ss_density_rho %.0f differs from AED rho_ss %.0f"
                                " (AED settling still uses %.0f)\n", name, rho, rho_ss[grp], rho_ss[grp]);
        } else if ( grp >= 0 ) {
            rho = rho_ss[grp];
            src = "aed rho_ss";
        }

        if ( rho <= 0. ) {
            fprintf(stderr, "     WARNING: no particle density for '%s' (set ss_density_rho); ignored\n", name);
            continue;
        }
        if ( rho <= 1000. )
            fprintf(stderr, "     WARNING: %s: particle density %.0f is not denser than water\n", name, rho);

        ss_idx[n_ss] = idx;
        ss_rho[n_ss] = rho;
        n_ss++;

        fprintf(stderr, "        %-10s rho_s = %6.0f kg/m3 (%s)", name, rho, src);
        if ( grp >= 0 ) {
            if ( settling == 3 ) {
                //# Stokes velocity at 15 degC for orientation only
                AED_REAL w = 9.807 * d_ss[grp] * d_ss[grp] * (rho_ss[grp] - 999.1) / (18. * 1.138e-3) * SecsPerDay;
                fprintf(stderr, ", d_ss = %.3g m, Stokes ~%.3g m/d at 15C", d_ss[grp], w);
            } else if ( settling == 1 || settling == 2 )
                fprintf(stderr, ", w_ss = %.3g m/d", w_ss[grp]);
        }
        fprintf(stderr, "\n");
    }

    if ( found > 0 && settling == 0 )
        fprintf(stderr, "     WARNING: &%s settling = 0 - particles never settle, so their"
                        " density excess never decays\n", NCS_GROUP);
    if ( mobility_off )
        fprintf(stderr, "     WARNING: mobility_off = .true. - WQ settling is disabled, so the SS"
                        " density excess never decays\n");

    if ( n_ss == 0 ) {
        fprintf(stderr, "     WARNING: no usable ss_density variables; SS density switched OFF\n");
        ss_density = FALSE; ss_density_diag = FALSE;
    }

    if ( ss_density_diag ) {
        snprintf(diag_path, sizeof(diag_path), "%s/inflow_ss_density.csv",
                                            (out_dir != NULL && *out_dir) ? out_dir : ".");
        fprintf(stderr, "     SS density: inflow diagnostics to %s\n", diag_path);
    }
}


/******************************************************************************
 * Density excess [kg/m3] of the SS in a WQ vector (e.g. an inflow parcel)    *
 ******************************************************************************/
static AED_REAL ss_excess_any(const AED_REAL *wq, AED_REAL rho_w)
{
    AED_REAL drho = zero;
    int i;

    if ( wq == NULL ) return zero;
    for (i = 0; i < n_ss; i++)
        if ( wq[ss_idx[i]] > zero )
            drho += wq[ss_idx[i]] * 1e-3 * (1.0 - rho_w / ss_rho[i]);
    return drho;
}

AED_REAL ss_density_excess(const AED_REAL *wq, AED_REAL rho_w)
{
    if ( !ss_density ) return zero;
    return ss_excess_any(wq, rho_w);
}


/******************************************************************************
 * Density excess [kg/m3] of the SS in lake layer k                           *
 ******************************************************************************/
AED_REAL ss_layer_excess(int k, AED_REAL rho_w)
{
    AED_REAL drho = zero;
    int i;

    if ( !ss_density ) return zero;
    for (i = 0; i < n_ss; i++)
        if ( _WQ_Vars(ss_idx[i], k) > zero )
            drho += _WQ_Vars(ss_idx[i], k) * 1e-3 * (1.0 - rho_w / ss_rho[i]);
    return drho;
}


/******************************************************************************
 * Density of lake layer k from T, S and (when on) its SS                     *
 ******************************************************************************/
AED_REAL layer_density(int k)
{
    AED_REAL rho = calculate_density(Lake[k].Temp, Lake[k].Salinity);
    return rho + ss_layer_excess(k, rho);
}


/******************************************************************************
 * Recalculate the density of every lake layer (e.g. after WQ settling)       *
 ******************************************************************************/
void refresh_densities()
{
    int k;

    if ( !ss_density ) return;
    for (k = botmLayer; k <= surfLayer; k++)
        Lake[k].Density = layer_density(k);
}


/******************************************************************************
 * Inflow diagnostics                                                         *
 ******************************************************************************/
void ss_density_set_day(int jday) { diag_jday = jday; }

/******************************************************************************
 * Write one row of output/inflow_ss_density.csv:                             *
 *   stage  "inflow" - the day's inflow water as it enters at the surface     *
 *          "subm"   - a submerged inflow inserted at its fixed height        *
 *          "insert" - a parcel inserted at its level of neutral buoyancy     *
 *                     (or the bottom), after plunging and entraining         *
 *   volume is the parcel volume [m3]; comparing "insert" with "inflow"       *
 *   volumes shows how much lake water was entrained on the way down.        *
 *   height_m is the middle of the receiving layer (m above the bottom).      *
 ******************************************************************************/
void ss_density_diag_inflow(int iRiver, const char *stage, AED_REAL volume,
                            AED_REAL temp, AED_REAL salt, const AED_REAL *wq,
                            int layer)
{
    AED_REAL rho_ts, drho, ss = zero, height;
    char ts[24];
    int i;

    if ( !ss_density_diag || diag_path[0] == 0 || n_ss == 0 ) return;

    if ( diag_f == NULL ) {
        if ( (diag_f = fopen(diag_path, "w")) == NULL ) {
            fprintf(stderr, "     WARNING: cannot open %s; inflow diagnostics off\n", diag_path);
            diag_path[0] = 0;
            return;
        }
        fprintf(diag_f, "time,inflow,stage,volume_m3,temp,salt,ss_gm3,"
                        "rho_ts,rho,drho_ss,height_m,layer\n");
    }

    //# drho_ss is the excess the SS carries (or would carry, when ss_density
    //# is off); rho is the density GLM actually used for this parcel.
    rho_ts = calculate_density(temp, salt);
    drho   = ss_excess_any(wq, rho_ts);
    for (i = 0; i < n_ss; i++) ss += wq[ss_idx[i]];

    if ( layer < botmLayer ) layer = botmLayer;
    height = (layer > botmLayer) ? (Lake[layer].Height + Lake[layer-1].Height) / 2.0
                                 : Lake[layer].Height / 2.0;

    write_time_string(ts, diag_jday, 0);
    ts[10] = 0;   //# date only - inflows are applied once per day

    fprintf(diag_f, "%s,%d,%s,%.6g,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.3f,%d\n",
                    ts, iRiver+1, stage, volume, temp, salt, ss,
                    rho_ts, ss_density ? rho_ts + drho : rho_ts, drho, height, layer);
    fflush(diag_f);
}
