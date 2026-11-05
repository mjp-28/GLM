# Suspended-solids density in GLM (`ss-density` branch)

This branch lets suspended sediment add to water density in GLM. Stock GLM
calculates density from temperature and salinity only, so a sediment-laden
storm inflow behaves like clear water at the same temperature. With this
branch the sediment's mass is counted. Muddy inflows plunge deeper, the
layers they reach stay denser until the sediment settles out, and the
extra density then fades as AED settles the particles.

It was written for Lake Tutira (Cyclone Gabrielle, February 2023), where a
warm, sediment-laden underflow reached the bottom of the lake.

## What it does

Density with the feature on:

    rho = rho(T, S) + sum_i  C_i * 1e-3 * (1 - rho(T, S) / rho_s,i)

- `C_i`: concentration of sediment group *i* in g/m3. By default these are
  AED's non-cohesive groups `NCS_ss1`, `NCS_ss2`, and so on.
- `rho_s,i`: particle (grain) density in kg/m3, read from `rho_ss` in
  `&aed_noncohesive` in the AED namelist. AED and GLM therefore share one
  value.

It is applied to:

- **Inflows.** The plunge point, each entrainment step on the way down
  (the parcel's sediment is diluted by entrained lake water), and the
  insertion depth. Submerged inflows are included.
- **Lake layers.** Every place GLM recalculates a layer's density: inflow
  insertion, layer merging, surface and deep mixing, stability checks,
  rain and evaporation, and re-injection. The surface mixed layer carries
  the sediment's density as it deepens.
- **After every AED step.** Densities are refreshed once AED has settled
  the particles. When a warm bottom layer loses its sediment it becomes
  lighter than the water above, and GLM's stability check mixes it upward.

With `ss_density = .false.` (the default) every density is calculated
exactly as in stock GLM.

The branch also fixes an existing bug: water-quality dilution weights during
inflow entrainment (commit "Fix WQ dilution weights during inflow
entrainment"). This fix changes results even with `ss_density` off. It can
be turned off to reproduce older runs.

## Settings

### `glm.nml` / `glm4.nml`

```
&inflow
   ...
   inflow_wq_entrain_fix = .true.   ! default .true.; .false. = old (buggy) dilution
/

&wq_setup
   wq_lib      = 'aed'
   wq_nml_file = 'aed/aed.nml'
   ...
   ss_density      = .true.                          ! default .false.
   ss_density_vars = 'NCS_ss1','NCS_ss2','NCS_ss3'   ! optional; default = all NCS_ss* groups
!  ss_density_rho  = 2650, 2650, 2650                ! optional override [kg/m3]; default = rho_ss from AED
   ss_density_diag = .true.                          ! optional; writes inflow_ss_density.csv
/
```

- AED must be running. `ss_density` needs the sediment to be AED state
  variables. Without AED, GLM prints a warning and switches the feature off.
- `&wq_setup` must still contain the usual `wq_lib` and `wq_nml_file`. GLM
  turns water quality on as soon as `&wq_setup` contains any setting, so
  adding only `ss_density` would switch AED on with the default file name.
- `ss_density_rho` is for sensitivity runs. If it differs from AED's
  `rho_ss`, GLM warns you, because AED's Stokes settling still uses
  `rho_ss`.

### `aed.nml` (no new entries)

The sediment groups are AED's existing non-cohesive module:

```
&aed_noncohesive
   num_ss   = 3
   settling = 3                        ! 3 = Stokes; 1 = constant w_ss. NOT 0 (= never settles)
   rho_ss   = 2650, 2650, 2650         ! grain density, kg/m3 (AED default 1600 is a floc value)
   d_ss     = 1.25e-6, 1.5e-5, 1.0e-4  ! particle diameter, m - controls how fast the excess decays
   ...
/
```

The inflow files need `NCS_ss1`, `NCS_ss2`, ... columns in g/m3, listed in
`inflow_vars`. Keep `mobility_off = .false.` (the default), otherwise
nothing settles.

At start-up GLM prints what it is using, for example:

```
     SS density: suspended solids included in inflow and layer density
        NCS_ss1    rho_s =   2650 kg/m3 (aed rho_ss), d_ss = 1.25e-06 m, Stokes ~0.107 m/d at 15C
        NCS_ss2    rho_s =   2650 kg/m3 (aed rho_ss), d_ss = 1.5e-05 m, Stokes ~15.4 m/d at 15C
```

It warns about `settling = 0`, `mobility_off`, unknown variable names and
missing densities.

## Diagnostics

These are only produced when `ss_density` is on.

| Output | Content |
|---|---|
| `output.nc`: `ss_rho_excess` | Per-layer density added by sediment, kg/m3. `dens` already includes it. |
| Point CSVs: `ssrho` column | The same quantity at your `csv_point_at` depths. Add `'ssrho'` to `csv_point_vars` and count it in `csv_point_nvars`. |
| `<out_dir>/inflow_ss_density.csv` | One row per inflow event when `ss_density_diag = .true.` (see below). |

`inflow_ss_density.csv` columns:

    time,inflow,stage,volume_m3,temp,salt,ss_gm3,rho_ts,rho,drho_ss,height_m,layer

- `stage`:
  - `inflow`: the day's inflow as it enters at the surface;
  - `insert`: the parcel at the layer where it was inserted, after
    plunging and entraining lake water (compare its `volume_m3` with
    the `inflow` row for the entrained volume);
  - `subm`: a submerged inflow.
- `rho_ts` excludes sediment and `rho` includes it.
- `height_m` is the middle of the receiving layer, in metres above the
  bottom.

## Building on Windows

GLM is normally built from AED's **GLM-AED** bundle, which contains GLM and
all the libraries it needs:

1. Install the tools listed in `README.Windows` (chocolatey: `make`,
   `mingw`, `git`, `cmake`).
2. In Git Bash:
   ```
   git clone --recurse-submodules https://github.com/AquaticEcoDynamics/GLM-AED.git
   cd GLM-AED/glm-source/GLM
   git remote add fork https://github.com/mjp-28/GLM.git
   git fetch fork
   git checkout -b ss-density fork/ss-density
   cd ..
   ./build_glm.sh
   ```
3. The new source file `src/glm_ss_density.c` is listed in the `Makefile`.

The Visual Studio project files in `win/` also list the new file. However,
they were already missing several existing GLM sources (oxygenation,
bubbler, heat exchange, restart) before this branch, so use the MinGW /
`build_glm.sh` route above.

## Testing checklist

1. **Switch off.** Run an existing model with `ss_density = .false.`. Results
   should match a stock GLM build, apart from the entrainment fix. To check
   that too, also set `inflow_wq_entrain_fix = .false.`; results should then
   be identical.
2. **On, with no sediment.** Set `ss_density = .true.` with all `NCS_ss*`
   inflow columns at 0. Results should match run 1, and `drho_ss` in
   `inflow_ss_density.csv` should be 0.
3. **On, with real sediment.** For an event such as Gabrielle:
   - `insert` rows should sit lower (smaller `height_m`) than in run 2;
   - `ss_rho_excess` near the bottom should rise and then decay at a rate
     set by `d_ss`;
   - the warm bottom water should mix upward as it decays.

## Known simplifications

- Nothing settles while a parcel is still plunging; AED only settles
  particles once they are in a lake layer. Coarse material therefore
  counts fully during the plunge. A separate box model of stream-mouth
  deposition is planned to handle this.
- Inflows are inserted once per day using the day's mean concentration.
- `resize_internals()` merges the surface layer's T and S but not its WQ.
  This is existing GLM behaviour and is left unchanged.
- Re-injected water (heat exchange, oxygenation recirculation) and
  bubbler plumes keep T/S-only densities.
