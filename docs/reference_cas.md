# Référence des clés du fichier de cas

Document généré à partir de la liste des clés que le logiciel vérifie (`microrans/fv2d/validate.py`) : il contient exactement les clés reconnues. Ne pas le modifier à la main ; le régénérer avec
`python -m microrans.fv2d.validate > docs/reference_cas.md`.

Toute autre clé est signalée (« clé inconnue, ignorée — vouliez-vous dire … ? »). « incompressible seulement » / « compressible seulement » : clé lue par un seul des deux solveurs (`[physics] compressible = true`) ; « maillage … » : types de maillage qui l'utilisent. Les valeurs par défaut de `[solver]` sont lues dans le code (`Settings`, `CompressibleSettings`).

Exemples complets : `microrans examples` ; tutoriel : `docs/tutoriel.md`.

## `[mesh]` — maillage

| Clé | Signification | S'applique à |
|---|---|---|
| `type` | rectangle \| blocks \| ogrid \| unstructured \| hybrid \| file \| box (défaut unstructured ; box : pavé 3D) |  |
| `preset` | maillage prédéfini (cavity, channel, cylinder-ogrid…) : remplace type |  |
| `cut_axis` | garde la moitié y > 0 (axisymétrique autour d'un corps) |  |
| `x0` | x min | maillage rectangle, box |
| `x1` | x max | maillage rectangle, box |
| `y0` | y min | maillage rectangle, box |
| `y1` | y max | maillage rectangle, box |
| `z0` | z min | maillage box |
| `z1` | z max | maillage box |
| `nx` | nombre de mailles selon x | maillage rectangle, box |
| `ny` | nombre de mailles selon y | maillage rectangle, box |
| `nz` | nombre de mailles selon z | maillage box |
| `grading` | resserrement [gx, gy] (box : [gx, gy, gz] ; rapport dernière / première maille) | maillage rectangle, box |
| `patch_types` | types des frontières {nom = "wall" \| "patch" \| …} | maillage rectangle, file, box |
| `periodic` | paires de frontières périodiques [["a", "b"], …] | maillage rectangle, blocks, box |
| `vertices` | sommets [[x, y], …] | maillage blocks |
| `n_around` | mailles autour du corps (défaut 128) | maillage ogrid |
| `n_radial` | mailles dans la direction radiale (défaut 64) | maillage ogrid |
| `farfield_radius` | rayon du champ lointain (défaut 20) | maillage ogrid |
| `first_height` | hauteur de la 1re maille à la paroi (défaut 1e-3) | maillage ogrid |
| `center` | centre du maillage en O (défaut : centre du corps) | maillage ogrid |
| `h_max` | taille maximale des triangles (défaut 1) | maillage unstructured, hybrid |
| `h_surface` | taille des mailles sur les corps (défaut 0.05) | maillage unstructured, hybrid |
| `growth` | croissance de la taille avec la distance au corps (défaut 0.2) | maillage unstructured, hybrid |
| `max_iter` | itérations du mailleur DistMesh (défaut 300) | maillage unstructured, hybrid |
| `path` | fichier .msh (Gmsh) ou .su2 | maillage file |

### `[mesh.names]` — noms des frontières du rectangle (du pavé : + back, front) (maillage rectangle, box)

| Clé | Signification | S'applique à |
|---|---|---|
| `left` | nom du côté x = x0 |  |
| `right` | nom du côté x = x1 |  |
| `bottom` | nom du côté y = y0 |  |
| `top` | nom du côté y = y1 |  |
| `back` | box : nom du côté z = z0 |  |
| `front` | box : nom du côté z = z1 |  |

### `[mesh.extrude]` — 3D : extrusion du maillage 2D selon z (quadrilatères → hexaèdres, triangles → prismes)

| Clé | Signification | S'applique à |
|---|---|---|
| `z0` | z de départ (défaut 0) |  |
| `z1` | z d'arrivée (défaut 1) |  |
| `nz` | nombre de couches (défaut 1) |  |
| `grading` | resserrement selon z (rapport dernière / première couche) |  |
| `patch_types` | types des frontières {nom = "wall" \| "symmetry" \| …} (défaut : type 2D ; patch pour back et front) |  |
| `periodic` | paires périodiques supplémentaires [["back", "front"]] |  |

#### `[mesh.extrude.names]` — noms des faces d'extrémité

| Clé | Signification | S'applique à |
|---|---|---|
| `back` | nom de la face z = z0 (défaut back) |  |
| `front` | nom de la face z = z1 (défaut front) |  |

### `[[mesh.blocks]]` — blocs (maillage blocks)

| Clé | Signification | S'applique à |
|---|---|---|
| `vertices` | 4 indices de sommets (sens trigonométrique) |  |
| `cells` | [nx, ny] |  |
| `grading` | [gx, gy] |  |

### `[[mesh.edges]]` — arêtes courbes (maillage blocks)

| Clé | Signification | S'applique à |
|---|---|---|
| `type` | arc \| polyline \| spline |  |
| `vertices` | [i, j] : arête courbe entre deux sommets |  |
| `point` | arc : point intermédiaire [x, y] |  |
| `points` | polyline, spline : points intermédiaires |  |

### `[mesh.patches]` — frontières nommées (maillage blocks)

#### `[mesh.patches.<nom>]` — une section par nom (frontière, scalaire…)

| Clé | Signification | S'applique à |
|---|---|---|
| `type` | wall \| patch \| symmetry \| empty |  |
| `faces` | [[i, j], …] : paires de sommets |  |

### `[[mesh.refinements]]` — zones de raffinement (maillage unstructured, hybrid)

| Clé | Signification | S'applique à |
|---|---|---|
| `h` | taille des mailles dans la zone |  |
| `growth` | croissance hors de la zone (défaut 0.2) |  |

#### `[mesh.refinements.shape]` — zone (même syntaxe que les corps)

| Clé | Signification | S'applique à |
|---|---|---|
| `type` | circle \| rectangle \| ellipse \| polygon \| naca \| spline \| file |  |
| `name` | nom du corps = nom de la frontière (condition aux limites) |  |
| `patch_type` | type de frontière (défaut wall pour un corps) |  |
| `center` | centre [x, y] (circle, ellipse ; défaut [0, 0]) |  |
| `radius` | rayon (circle) |  |
| `x0` | rectangle : x min |  |
| `y0` | rectangle : y min |  |
| `x1` | rectangle : x max |  |
| `y1` | rectangle : y max |  |
| `a` | ellipse : demi-axe selon x |  |
| `b` | ellipse : demi-axe selon y |  |
| `points` | polygon, spline : [[x, y], …] |  |
| `code` | naca : 4 chiffres (défaut "0012") |  |
| `chord` | naca : corde (défaut 1) |  |
| `n` | naca, spline : nombre de points du contour |  |
| `closed_te` | naca : bord de fuite fermé (défaut true) |  |
| `path` | file : fichier de contour (.dat Selig/Lednicer, .csv, .dxf) |  |
| `sharp_angle` | file : angle (°) au-delà duquel un sommet est un coin (défaut 60) |  |
| `angle` | rotation (°, sens trigonométrique) |  |
| `incidence` | incidence (°) : rotation de −incidence (profil cabré) |  |
| `rotation_center` | centre de rotation [x, y] (défaut [0, 0]) |  |
| `scale` | facteur d'échelle (défaut 1) |  |
| `translate` | translation [dx, dy] |  |

##### `[mesh.refinements.shape.names]` — noms des frontières du rectangle

| Clé | Signification | S'applique à |
|---|---|---|
| `left` | nom du côté x = x0 |  |
| `right` | nom du côté x = x1 |  |
| `bottom` | nom du côté y = y0 |  |
| `top` | nom du côté y = y1 |  |

### `[mesh.layers]` — couches de quadrilatères aux parois (maillage hybrid)

| Clé | Signification | S'applique à |
|---|---|---|
| `n` | nombre de couches (défaut 10) |  |
| `first_height` | hauteur de la 1re couche (défaut 1e-3) |  |
| `ratio` | rapport de croissance (défaut 1.2) |  |

## `[domain]` — domaine extérieur (défaut : rectangle [-10, 30] × [-10, 10]) (maillage unstructured, hybrid)

| Clé | Signification | S'applique à |
|---|---|---|
| `type` | circle \| rectangle \| ellipse \| polygon \| naca \| spline \| file |  |
| `name` | nom du corps = nom de la frontière (condition aux limites) |  |
| `patch_type` | type de frontière (défaut wall pour un corps) |  |
| `center` | centre [x, y] (circle, ellipse ; défaut [0, 0]) |  |
| `radius` | rayon (circle) |  |
| `x0` | rectangle : x min |  |
| `y0` | rectangle : y min |  |
| `x1` | rectangle : x max |  |
| `y1` | rectangle : y max |  |
| `a` | ellipse : demi-axe selon x |  |
| `b` | ellipse : demi-axe selon y |  |
| `points` | polygon, spline : [[x, y], …] |  |
| `code` | naca : 4 chiffres (défaut "0012") |  |
| `chord` | naca : corde (défaut 1) |  |
| `n` | naca, spline : nombre de points du contour |  |
| `closed_te` | naca : bord de fuite fermé (défaut true) |  |
| `path` | file : fichier de contour (.dat Selig/Lednicer, .csv, .dxf) |  |
| `sharp_angle` | file : angle (°) au-delà duquel un sommet est un coin (défaut 60) |  |
| `angle` | rotation (°, sens trigonométrique) |  |
| `incidence` | incidence (°) : rotation de −incidence (profil cabré) |  |
| `rotation_center` | centre de rotation [x, y] (défaut [0, 0]) |  |
| `scale` | facteur d'échelle (défaut 1) |  |
| `translate` | translation [dx, dy] |  |

### `[domain.names]` — noms des frontières du rectangle

| Clé | Signification | S'applique à |
|---|---|---|
| `left` | nom du côté x = x0 |  |
| `right` | nom du côté x = x1 |  |
| `bottom` | nom du côté y = y0 |  |
| `top` | nom du côté y = y1 |  |

## `[[bodies]]` — corps (obstacles) (maillage ogrid, unstructured, hybrid)

| Clé | Signification | S'applique à |
|---|---|---|
| `type` | circle \| rectangle \| ellipse \| polygon \| naca \| spline \| file |  |
| `name` | nom du corps = nom de la frontière (condition aux limites) |  |
| `patch_type` | type de frontière (défaut wall pour un corps) |  |
| `center` | centre [x, y] (circle, ellipse ; défaut [0, 0]) |  |
| `radius` | rayon (circle) |  |
| `x0` | rectangle : x min |  |
| `y0` | rectangle : y min |  |
| `x1` | rectangle : x max |  |
| `y1` | rectangle : y max |  |
| `a` | ellipse : demi-axe selon x |  |
| `b` | ellipse : demi-axe selon y |  |
| `points` | polygon, spline : [[x, y], …] |  |
| `code` | naca : 4 chiffres (défaut "0012") |  |
| `chord` | naca : corde (défaut 1) |  |
| `n` | naca, spline : nombre de points du contour |  |
| `closed_te` | naca : bord de fuite fermé (défaut true) |  |
| `path` | file : fichier de contour (.dat Selig/Lednicer, .csv, .dxf) |  |
| `sharp_angle` | file : angle (°) au-delà duquel un sommet est un coin (défaut 60) |  |
| `angle` | rotation (°, sens trigonométrique) |  |
| `incidence` | incidence (°) : rotation de −incidence (profil cabré) |  |
| `rotation_center` | centre de rotation [x, y] (défaut [0, 0]) |  |
| `scale` | facteur d'échelle (défaut 1) |  |
| `translate` | translation [dx, dy] |  |

### `[bodies.names]` — noms des frontières du rectangle

| Clé | Signification | S'applique à |
|---|---|---|
| `left` | nom du côté x = x0 |  |
| `right` | nom du côté x = x1 |  |
| `bottom` | nom du côté y = y0 |  |
| `top` | nom du côté y = y1 |  |

## `[physics]` — physique

| Clé | Signification | S'applique à |
|---|---|---|
| `compressible` | true : solveur compressible (grandeurs SI, section [flow]) |  |
| `nu` | viscosité cinématique ν (m²/s) | incompressible seulement |
| `reynolds` | nombre de Reynolds (ν = U_ref L_ref / Re) — au lieu de nu | incompressible seulement |
| `reference_velocity` | vitesse de référence U_ref (défaut 1) | incompressible seulement |
| `reference_length` | longueur de référence L_ref (défaut 1) |  |
| `reference_area` | axisymétrique : aire de référence de Cd (défaut π L²/4) | incompressible seulement |
| `model` | laminar \| sa \| ke \| kw \| sst \| sst_gamma (défaut laminar) | incompressible seulement |
| `body_force` | force volumique [fx, fy] (3D : [fx, fy, fz] ; conduite périodique) | incompressible seulement |
| `angle_of_attack` | incidence (°) de l'écoulement amont |  |
| `axisymmetric` | true : axisymétrique (x = axe, y = rayon) |  |
| `swirl` | true : rotation propre (axisymétrique) | incompressible seulement |

### `[physics.model_options]` — options du modèle (sa : ft2 ; sst_gamma : kato_launder) (incompressible seulement)

Contenu libre, sous la forme `nom = valeur`.

### `[physics.viscosity]` — fluide non newtonien (incompressible seulement)

| Clé | Signification | S'applique à |
|---|---|---|
| `model` | newtonian \| power_law \| carreau \| cross \| herschel_bulkley \| bingham \| casson |  |
| `K` | indice de consistance |  |
| `n` | indice de comportement |  |
| `nu0` | viscosité à cisaillement nul |  |
| `nu_inf` | viscosité à cisaillement infini |  |
| `lambda` | carreau : temps caractéristique |  |
| `a` | carreau : exposant (défaut 2) |  |
| `m` | cross : temps caractéristique |  |
| `tau_y` | contrainte seuil |  |
| `nu_min` | borne inférieure de ν |  |
| `nu_max` | borne supérieure de ν |  |
| `relax` | sous-relaxation de ν (défaut 0.5) |  |
| `gamma_ref` | taux de cisaillement de référence (défaut U_ref / L_ref) |  |

## `[flow]` — écoulement amont (compressible) (compressible seulement)

| Clé | Signification | S'applique à |
|---|---|---|
| `mach` | nombre de Mach amont |  |
| `velocity` | vitesse [u, v] (m/s) |  |
| `pressure` | pression (Pa) |  |
| `temperature` | température (K) |  |
| `density` | masse volumique (kg/m³) |  |
| `gamma` | γ (défaut 1.4) |  |
| `gas_constant` | R (J/kg/K, défaut 287.058) |  |
| `prandtl` | Pr (défaut 0.72) |  |
| `viscosity` | inviscid \| constant \| sutherland |  |
| `mu` | viscosité (Pa·s) |  |
| `reynolds` | Re = ρ U L_ref / μ (au lieu de mu) |  |
| `T_ref` | Sutherland : température de référence (défaut 273.15 K) |  |
| `sutherland_S` | Sutherland : constante S (défaut 110.4 K) |  |

## `[initial]` — état initial

| Clé | Signification | S'applique à |
|---|---|---|
| `U` | vitesse initiale [ux, uy] (3D : [ux, uy, uz] ; valeurs ou formules en x, y, z) |  |
| `perturbation` | amplitude d'un tourbillon initial (déclenche le lâcher) | incompressible seulement |
| `perturbation_center` | centre du tourbillon initial (défaut [1.5, 0]) | incompressible seulement |
| `restart` | fichier checkpoint.npz de reprise |  |
| `restart_mode` | exact (défaut) \| fields |  |
| `restart_shift_U` | décalage de vitesse à la reprise | incompressible seulement |
| `rho` | masse volumique (formule possible) | compressible seulement |
| `density` | = rho | compressible seulement |
| `p` | pression (formule possible) | compressible seulement |
| `pressure` | = p | compressible seulement |
| `T` | température (formule possible) | compressible seulement |
| `temperature` | = T | compressible seulement |
| `mach` | Mach initial (direction amont) | compressible seulement |

## `[turbulence]` — turbulence amont (incompressible seulement)

| Clé | Signification | S'applique à |
|---|---|---|
| `intensity` | intensité turbulente amont (défaut 0.001) |  |
| `viscosity_ratio` | ν_t/ν amont (défaut 0.1) |  |

## `[energy]` — thermique (incompressible seulement)

| Clé | Signification | S'applique à |
|---|---|---|
| `enabled` | false : thermique désactivée |  |
| `Pr` | nombre de Prandtl |  |
| `Pr_t` | Prandtl turbulent (défaut 0.85) |  |
| `beta` | dilatation thermique (Boussinesq) |  |
| `T_ref` | température de référence |  |
| `gravity` | gravité [gx, gy] (3D : [gx, gy, gz]) |  |
| `T0` | température initiale |  |
| `delta_T` | écart de température de référence (nombre de Nusselt) |  |
| `source` | source de chaleur (valeur ou formule) |  |

## `[scalars]` — scalaires passifs [scalars.<nom>] (incompressible seulement)

### `[scalars.<nom>]` — une section par nom (frontière, scalaire…)

| Clé | Signification | S'applique à |
|---|---|---|
| `diffusivity` | diffusivité D |  |
| `schmidt` | Sc = ν / D (au lieu de diffusivity) |  |
| `Sc_t` | Schmidt turbulent (défaut 0.7) |  |
| `source` | source (formule) |  |
| `initial` | valeur initiale |  |
| `scheme` | upwind \| linearUpwind \| linearUpwindLimited (défaut) |  |

## `[[porous]]` — zones poreuses (incompressible seulement)

| Clé | Signification | S'applique à |
|---|---|---|
| `name` | nom de la zone |  |
| `region` | rectangle \| circle \| expression |  |
| `x0` | rectangle : x min |  |
| `x1` | rectangle : x max |  |
| `y0` | rectangle : y min |  |
| `y1` | rectangle : y max |  |
| `center` | circle : centre |  |
| `radius` | circle : rayon |  |
| `expression` | condition en x, y |  |
| `darcy` | coefficient de Darcy (1/m²) |  |
| `permeability` | perméabilité K = 1/darcy (m²) |  |
| `forchheimer` | coefficient de Forchheimer (1/m) |  |
| `angle` | orientation des axes principaux (°) |  |

## `[[actuator_disk]]` — disques actuateurs (incompressible seulement)

| Clé | Signification | S'applique à |
|---|---|---|
| `name` | nom |  |
| `x0` | début (x) |  |
| `x1` | fin (x) |  |
| `radius` | rayon |  |
| `hub_radius` | rayon du moyeu |  |
| `thrust` | poussée |  |
| `thrust_coefficient` | coefficient de poussée |  |
| `mode` | turbine \| propeller |  |
| `torque` | couple |  |
| `distribution` | uniform \| optimal |  |
| `center` | centre |  |

## `[boundary]` — conditions aux limites [boundary.<frontière>]

### `[boundary.<nom>]` — une section par nom (frontière, scalaire…)

| Clé | Signification | S'applique à |
|---|---|---|
| `type` | type de condition (incompressible : wall, inlet, outlet, symmetry, farfield, axis, pressure_inlet ; compressible : farfield, inlet, outlet, supersonic_inlet, supersonic_outlet, slip_wall, symmetry, wall) |  |
| `U` | vitesse [ux, uy] (3D : [ux, uy, uz] ; valeurs ou formules en x, y, z) |  |
| `p` | pression |  |
| `p0` | pression totale |  |
| `T` | température |  |
| `q` | flux de chaleur pariétal | incompressible seulement |
| `flow_rate` | débit (au lieu de U) | incompressible seulement |
| `profile` | profil du débit : uniform \| parabolic | incompressible seulement |
| `U_theta` | vitesse de rotation (swirl) | incompressible seulement |
| `omega` | paroi tournante (rad/s, swirl) | incompressible seulement |
| `mach` | Mach | compressible seulement |
| `pressure` | = p | compressible seulement |
| `temperature` | = T | compressible seulement |
| `density` | masse volumique | compressible seulement |
| `velocity` | vitesse [u, v] | compressible seulement |
| `angle` | direction de l'écoulement (°) | compressible seulement |
| `direction` | direction [dx, dy] | compressible seulement |
| `total_pressure` | = p0 | compressible seulement |
| `T0` | température totale | compressible seulement |
| `total_temperature` | = T0 | compressible seulement |
| `vortex` | correction de tourbillon : point [x, y] | compressible seulement |

#### `[boundary.<nom>.scalars]` — valeurs des scalaires (incompressible seulement)

Contenu libre, sous la forme `nom = valeur`.

#### `[boundary.<nom>.scalar_flux]` — flux des scalaires (incompressible seulement)

Contenu libre, sous la forme `nom = valeur`.

## `[solver]` — réglages numériques

| Clé | Signification | Défaut (incompressible) | Défaut (compressible) | S'applique à |
|---|---|---|---|---|
| `adjust_dt` | pas de temps adaptatif (défaut false) | false |  | incompressible seulement |
| `algorithm` | SIMPLE \| SIMPLEC (défaut) \| coupled | `SIMPLEC` |  | incompressible seulement |
| `backend` | cpu (défaut) \| cuda \| rocm \| intel | `cpu` |  | incompressible seulement |
| `cfl` | nombre CFL |  | 0.8 | compressible seulement |
| `cfl_adapt` | implicite : réduction du CFL si le résidu croît |  | 1.2 | compressible seulement |
| `cfl_cuts` | implicite : divisions du plafond de CFL (défaut 3) |  | 3 | compressible seulement |
| `cfl_growth` | implicite : croissance du CFL |  | 1.1 | compressible seulement |
| `cfl_max` | implicite : CFL maximal |  | 1000 | compressible seulement |
| `cn_theta` | Crank-Nicolson : θ (défaut 0.5) | 0.5 |  | incompressible seulement |
| `convection_T` | upwind \| linearUpwind (défaut) \| linearUpwindLimited | `linearUpwind` |  | incompressible seulement |
| `convection_U` | upwind \| linearUpwind (défaut) \| linearUpwindLimited (gradient limité, plus robuste sur maillage déformé) | `linearUpwind` |  | incompressible seulement |
| `convection_turb` | upwind (défaut) \| linearUpwind \| linearUpwindLimited | `upwind` |  | incompressible seulement |
| `ddt_phi_coeff` | correction ddtCorr (défaut 1) | 1 |  | incompressible seulement |
| `dt` | pas de temps (instationnaire) |  |  |  |
| `entropy_fix` | correction d'entropie de Harten (défaut 0.1) |  | 0.1 | compressible seulement |
| `first_order_iter` | itérations d'ordre 1 au départ |  | 0 | compressible seulement |
| `flux` | roe (défaut) \| hllc |  | `roe` | compressible seulement |
| `fmg_levels` | démarrage multigrille : niveaux grossiers (défaut 0) |  |  | incompressible seulement |
| `fmg_tol` | démarrage multigrille : tolérance des niveaux grossiers |  |  | incompressible seulement |
| `implicit_jacobian` | roe \| rusanov |  | `roe` | compressible seulement |
| `limiter` | venkatakrishnan (défaut) \| barth_jespersen \| none |  | `venkatakrishnan` | compressible seulement |
| `limiter_freeze` | limiteur gelé après N itérations (défaut 0 : jamais) |  | 0 | compressible seulement |
| `linear_iter` | gmres : itérations |  | 20 | compressible seulement |
| `linear_solver` | gmres (défaut) \| sgs |  | `gmres` | compressible seulement |
| `linear_sweeps` | implicite : balayages SGS |  | 2 | compressible seulement |
| `linear_tol` | gmres : tolérance |  | 0.05 | compressible seulement |
| `log_every` | affichage tous les N itérations / pas (défaut 100) |  | 100 |  |
| `max_co` | Courant visé si adjust_dt (défaut 1) | 1 |  | incompressible seulement |
| `max_dt` | pas de temps maximal si adjust_dt | aucun |  | incompressible seulement |
| `max_iter` | itérations maximales (stationnaire) | 3000 | 5000 |  |
| `mode` | steady (défaut) \| transient |  |  |  |
| `monitor_tol` | arrêt quand les efforts varient de moins de monitor_tol | auto | auto |  |
| `monitor_window` | fenêtre (itérations) de monitor_tol | 100 | 200 |  |
| `n_corr` | corrections de pression (défaut 2) | 2 |  | incompressible seulement |
| `n_nonorth` | corrections non orthogonales (défaut 1) | 1 |  | incompressible seulement |
| `n_outer` | boucles externes PIMPLE (défaut 2) | 2 |  | incompressible seulement |
| `nonorth_limit` | limiteur de la correction non orthogonale (défaut 0.5) | 0.5 |  | incompressible seulement |
| `numba` | noyaux Numba (défaut false) | false |  | incompressible seulement |
| `order` | 1 \| 2 (défaut) |  | 2 | compressible seulement |
| `pseudo_cfl` | CFL de pseudo-temps (stationnaire pseudo-transitoire) | auto |  | incompressible seulement |
| `pseudo_dt` | pas de pseudo-temps (stationnaire pseudo-transitoire) | auto |  | incompressible seulement |
| `relax_T` | sous-relaxation de la température (défaut 0.9) | 0.9 |  | incompressible seulement |
| `relax_U` | sous-relaxation de la vitesse (défaut auto) | auto |  | incompressible seulement |
| `relax_p` | sous-relaxation de la pression (défaut 1) | 1 |  | incompressible seulement |
| `relax_scalar` | sous-relaxation des scalaires (défaut 1) | 1 |  | incompressible seulement |
| `relax_turb` | sous-relaxation de la turbulence (défaut 0.8) | 0.8 |  | incompressible seulement |
| `solver_U` | auto \| direct \| amg \| bicgstab \| pyamg | `auto` |  | incompressible seulement |
| `solver_p` | auto \| direct \| amg \| bicgstab \| pyamg | `auto` |  | incompressible seulement |
| `solver_turb` | auto \| direct \| amg \| bicgstab \| pyamg | `auto` |  | incompressible seulement |
| `steady_scheme` | rk3 (défaut) \| rk5 \| implicit |  | `rk3` | compressible seulement |
| `t_end` | temps final (instationnaire) |  |  |  |
| `threads` | fils de calcul Numba (défaut 1) | 1 |  | incompressible seulement |
| `time_scheme` | auto \| euler \| backward \| crankNicolson \| rk1 \| rk2 \| rk3 \| rk4 \| ab2 | `auto` |  | incompressible seulement |
| `tol` | critère de convergence sur les résidus | 1e-05 | 1e-06 |  |
| `turbulence_every_outer` | turbulence à chaque boucle externe (défaut false) | false |  | incompressible seulement |
| `venkat_k` | seuil du limiteur de Venkatakrishnan (défaut 0.05) |  | 0.05 | compressible seulement |
| `viscous_factor` | coefficient du pas de temps visqueux |  | 2 | compressible seulement |
| `wall_treatment` | resolved (défaut) \| wall_function | `resolved` |  | incompressible seulement |

## `[output]` — sorties

| Clé | Signification | S'applique à |
|---|---|---|
| `directory` | dossier des résultats |  |
| `probes` | sondes [[x, y], …] (3D : [[x, y, z], …]) |  |
| `forces` | frontières où calculer les efforts (défaut : parois) |  |
| `moment_center` | centre des moments (défaut [0, 0] ; 3D : [x, y, z], moment autour de l'axe z) |  |
| `checkpoint` | écrit checkpoint.npz (défaut true) |  |
| `checkpoint_minutes` | sauvegarde périodique (min, défaut 5) |  |
| `vtk` | écrit fields.vtk (défaut true) |  |
| `plots` | figures (défaut true) |  |
| `average_from` | instationnaire : moyennes à partir de t | incompressible seulement |
| `animate` | instationnaire : grandeur animée (vorticity, U_mag, p…) | incompressible seulement |
| `animate_every` | une image tous les N pas | incompressible seulement |
| `animate_fps` | images par seconde (défaut 15) | incompressible seulement |
| `vtk_every` | instationnaire : fields_N.vtk tous les N pas | incompressible seulement |
| `nusselt` | bulk : Nusselt local par température de mélange | incompressible seulement |

### `[[output.lines]]` — profils sur des segments

| Clé | Signification | S'applique à |
|---|---|---|
| `name` | nom |  |
| `start` | [x, y] (3D : [x, y, z]) |  |
| `end` | [x, y] (3D : [x, y, z]) |  |
| `n` | nombre de points |  |

## `[sweep]` — balayage

| Clé | Signification | S'applique à |
|---|---|---|
| `parameter` | clé balayée (ex. physics.reynolds) |  |
| `values` | liste de valeurs |  |
| `range` | [début, fin, pas] |  |
| `continuation` | part du point précédent (défaut true) |  |
| `jobs` | calculs en parallèle (défaut 1) |  |
