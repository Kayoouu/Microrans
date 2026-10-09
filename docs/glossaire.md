# Glossaire

Les termes de la mécanique des fluides numérique (CFD) employés par microrans, avec
l'endroit où on les rencontre dans le logiciel. Clés du fichier de cas :
[`reference_cas.md`](reference_cas.md) ; premier calcul : [`tutoriel.md`](tutoriel.md) ;
problèmes : [`depannage.md`](depannage.md).

## Écoulement

**Viscosité cinématique ν** (m²/s) — viscosité dynamique divisée par la masse volumique
(eau ≈ 10⁻⁶, air ≈ 1.5·10⁻⁵). Le solveur incompressible travaille en grandeurs divisées
par ρ (ρ = 1). Clé `[physics] nu`.

**Nombre de Reynolds Re = U L / ν** — rapport des effets d'inertie aux effets visqueux, avec
U et L une vitesse et une longueur de référence (`reference_velocity`,
`reference_length`). Petit (< ~1 000 pour un obstacle) : écoulement laminaire ; grand :
turbulent. Donner `[physics] reynolds` au lieu de `nu` : ν = U L / Re.

**Vitesse de référence U_ref** — une seule valeur, utilisée pour ν = U L / Re et pour tous
les coefficients (½ U² A). Si `reference_velocity` n'est pas donnée : vitesse d'entrée
(moyenne sur la frontière d'entrée, de débit ou de champ lointain ; la plus grande s'il y
en a plusieurs, avec un avertissement), sinon vitesse de la paroi mobile (couvercle d'une
cavité), sinon 1. Jamais la vitesse initiale. Valeur et origine affichées au début du
calcul et écrites dans `summary.json`.

**Laminaire / turbulent** — laminaire : écoulement régulier, calculé sans modèle ;
turbulent : fluctuations chaotiques, représentées en moyenne par un modèle de turbulence
(`[physics] model`).

**Incompressible / compressible** — incompressible : masse volumique constante, valable
tant que la vitesse reste petite devant celle du son (Mach < ~0.3). Compressible :
`[physics] compressible = true`, solveur en densité, grandeurs SI, section `[flow]`.

**Nombre de Mach M = U / c** — vitesse rapportée à celle du son c. Au-delà de 1
(supersonique) apparaissent des chocs.

**Stationnaire / instationnaire** — stationnaire : on cherche l'état qui ne varie plus
(`[solver] mode = "steady"`, itérations) ; instationnaire : on suit l'évolution dans le temps
(`mode = "transient"`, pas de temps `dt` jusqu'à `t_end`), par exemple un sillage oscillant.

**RANS / URANS** — équations de Navier-Stokes moyennées (Reynolds-Averaged Navier-Stokes) :
la turbulence est remplacée par une viscosité turbulente ν_t donnée par un modèle. URANS :
la même chose en instationnaire (grandes structures suivies dans le temps). microrans ne fait
ni LES ni DNS.

**Axisymétrique** — écoulement de révolution autour de l'axe x (tuyau, jet, sphère) : on
calcule un demi-plan (x = axe, y = rayon ≥ 0). `[physics] axisymmetric = true`, frontière
d'axe `type = "axis"`. **Rotation propre (swirl)** : vitesse azimutale u_θ en plus
(`swirl = true`).

**Fluide non newtonien** — viscosité qui dépend du cisaillement (sang, peintures, boues).
Lois disponibles : `[physics.viscosity] model` (loi puissance, Carreau, Cross,
Herschel-Bulkley, Bingham, Casson) ; laminaire uniquement.

**Scalaire passif** — grandeur transportée par l'écoulement sans l'influencer
(concentration, colorant, âge du fluide). `[scalars.<nom>]`.

**Boussinesq** — approximation de la convection naturelle : la masse volumique ne varie que
dans le terme de gravité, proportionnellement à l'écart de température (`[energy] beta`,
`gravity`, `T_ref`).

**Milieu poreux (Darcy-Forchheimer)** — zone (filtre, grillage, végétation) qui freine
l'écoulement par une force proportionnelle à la vitesse (Darcy) et à son carré
(Forchheimer). `[[porous]]`.

**Disque actuateur** — hélice ou éolienne représentée par une force répartie sur un disque,
sans dessiner les pales. `[[actuator_disk]]`.

## Maillage

**Maillage, cellule** — découpage du domaine en petites cellules (quadrilatères ou
triangles ; en 3D, hexaèdres ou prismes) ; le solveur calcule une valeur moyenne par
cellule (volumes finis). Plus il y a
de cellules, plus c'est précis… et lent.

**Frontière (patch)** — partie nommée du bord du domaine (« inlet », « wall »,
« cylinder »…) sur laquelle on impose une condition aux limites. Le nom d'un corps
(`[[bodies]] name`) est celui de sa frontière.

**Structuré / non structuré / hybride** — structuré : quadrilatères rangés en lignes et
colonnes (`rectangle`, `blocks`, maillage en O) ; non structuré : triangles
(`unstructured`) ; hybride : couches de quadrilatères le long des parois et triangles
ailleurs (`hybrid`).

**Maillage en O (O-grid)** — maillage structuré qui entoure un corps unique par des anneaux
successifs jusqu'au champ lointain (`type = "ogrid"`) : cylindre, profil d'aile.

**Couches de paroi, première maille** — mailles très fines et aplaties le long d'une paroi,
pour représenter la couche limite (`[mesh.layers]`, `first_height`).

**Qualité du maillage** — non-orthogonalité (angle entre la droite reliant deux centres de
cellules et la normale à leur face commune) et asymétrie (skewness). Plus elles sont
grandes, moins le calcul est précis et robuste (le seuil d'alerte usuel de
non-orthogonalité est ~70°). Affichées par `microrans mesh` (`quality.json`).

**Pavé, extrusion (3D)** — deux façons d'obtenir un maillage 3D : `[mesh] type = "box"`,
pavé d'hexaèdres (`x0` … `z1`, `nx`, `ny`, `nz`) ; `[mesh.extrude]`, n'importe quel
maillage 2D répété en couches selon z (`z0`, `z1`, `nz`). Les faces d'extrémité d'une
extrusion s'appellent `back` (z = z0) et `front` (z = z1) et ont besoin d'une condition aux
limites (ou d'être périodiques). Interface : type « Pavé 3D », ou case « Extruder le
maillage 2D en 3D » (tutoriel, § 7).

**Hexaèdre, prisme** — cellules 3D : l'hexaèdre est une « brique » à six faces
quadrilatères (pavé, quadrilatère extrudé) ; le prisme est un triangle extrudé (deux faces
triangulaires, trois quadrilatères).

**Frontières périodiques** — deux frontières traitées comme une seule : ce qui sort par
l'une entre par l'autre. Pour un écoulement établi (conduite, canal : `periodic =
[["inlet", "outlet"]]`) ou invariant selon z (`[mesh.extrude] periodic = [["back",
"front"]]`). Elles ne reçoivent pas de condition aux limites. Maillages `rectangle`,
`blocks`, `box` et faces z d'une extrusion seulement.

**Plan de coupe** — en 3D, les figures montrent les champs dans un plan x, y ou
z = constante : interface, page « 7. Résultats », « Plan de coupe (3D) » ; ligne de
commande, `[output] slice_axis` et `slice_value` (défaut : plan z médian). Le champ complet
est dans `fields.vtk` (ParaView).

## Parois et efforts

**Couche limite** — mince zone près d'une paroi où la vitesse passe de 0 (adhérence) à la
vitesse de l'écoulement.

**y⁺** — distance à la paroi en unités de paroi, y⁺ = u_τ y / ν (u_τ : vitesse de
frottement). Calcul résolu jusqu'à la paroi : il faut y⁺ ≈ 1 pour la première maille ;
lois de paroi : y⁺ ≈ 30 à 300. Valeurs obtenues : colonne `yplus` de `wall_<paroi>.csv`,
« y⁺ max » du résumé.

**Loi de paroi** — formule qui donne le frottement sans résoudre la couche limite, ce qui
permet des mailles plus épaisses à la paroi (`[solver] wall_treatment = "wall_function"`,
loi de Spalding). SA, k-ω, SST seulement.

**Cd, Cl, Cm** — coefficients de traînée (dans le sens de l'écoulement amont), de portance
(perpendiculaire) et de moment : effort divisé par ½ U² L (par unité de profondeur ; en
axisymétrique, effort total divisé par ½ U² A_ref). Résumé du calcul et `summary.json`.

**Cp, Cf** — coefficient de pression (p / ½ U²) et de frottement (τ_w / ½ U²) le long d'une
paroi : `wall_<paroi>.csv`, page « Résultats », tracés pariétaux.

**Polaire** — courbes Cl(α), Cd(α) d'un profil en fonction de l'incidence α
(`microrans polar`, ou balayage de `physics.angle_of_attack`).

**Nombre de Strouhal St = f L / U** — fréquence adimensionnée du lâcher de tourbillons
(≈ 0.16 derrière un cylindre à Re = 100). Calculé en instationnaire.

**Nombre de Nusselt** — flux de chaleur à une paroi rapporté au flux de conduction pure
(`[energy] delta_T`) ; **nombre de Prandtl** Pr = ν / α (α : diffusivité thermique).

## Modèles de turbulence

**Spalart-Allmaras (SA)** — une équation ; conçu pour l'aérodynamique externe (profils,
plaques). `model = "sa"`.

**k-ε** — deux équations (énergie turbulente k, dissipation ε) ; version Launder-Sharma
« bas Reynolds » : y⁺ ≈ 1 obligatoire, pas de lois de paroi. `model = "ke"`.

**k-ω, k-ω SST** — deux équations (k et ω = ε/k) ; SST (Menter) combine k-ω près des parois
et k-ε loin : choix courant en industrie. `model = "kw"`, `"sst"`.

**Transition (SST-γ)** — SST complété d'une équation d'intermittence γ qui prédit le passage
de la couche limite du laminaire au turbulent (Menter 2015). `model = "sst_gamma"`.

## Méthodes numériques

**Volumes finis** — méthode qui écrit la conservation (masse, quantité de mouvement…) sur
chaque cellule : ce qui entre par les faces moins ce qui sort.

**Schéma de convection** — façon d'estimer une grandeur sur une face à partir des cellules
voisines. `upwind` : valeur de la cellule amont (robuste, précision d'ordre 1, diffuse) ;
`linearUpwind` : corrigée par le gradient (ordre 2, défaut pour la vitesse) ;
`linearUpwindLimited` : gradient limité pour ne pas créer de nouveaux extrema (plus
robuste sur maillage déformé). `[solver] convection_U`.

**SIMPLE, SIMPLEC, couplé** — algorithmes qui résolvent ensemble vitesse et pression en
stationnaire : SIMPLE et SIMPLEC (défaut) les résolvent l'un après l'autre et corrigent ;
couplé (`algorithm = "coupled"`) les résout en un seul système : beaucoup moins
d'itérations, nombre presque indépendant de la taille du maillage. Sur les cas mesurés
(README, § 5.5), il est 1,4 à 27 fois plus rapide, sauf pour la transition SST-γ et la
convection naturelle où il est un peu plus lent ; mémoire mesurée de 0,3 Go
(16 000 cellules) à 1,65 Go (160 000 cellules).

**Sous-relaxation** — on n'applique qu'une fraction de chaque correction (0 < relax ≤ 1)
pour stabiliser le calcul stationnaire ; plus petit = plus stable mais plus lent.
`relax_U`, `relax_p`, `relax_turb`.

**Résidus, convergence** — les résidus mesurent combien les équations sont encore loin
d'être satisfaites ; le calcul stationnaire s'arrête quand ils passent sous `[solver] tol`
(« convergé »), sinon après `max_iter` itérations (« NON CONVERGÉ » : résultats à
vérifier). Courbes : `convergence.png`, `history.csv`.

**Convergence en maillage, GCI** — un calcul « convergé » (résidus) garde une erreur due au
maillage : on refait le cas sur trois maillages de plus en plus fins. L'indice de
convergence de maillage (GCI, procédure de Celik et al. 2008, `microrans.gci`) en tire
l'ordre apparent du schéma (≈ 2 ici), une valeur extrapolée à maillage infini et une bande
d'incertitude du maillage fin. Il ne couvre ni la taille du domaine ni le modèle de
turbulence. Exemples chiffrés : README § 6, « Incertitude de maillage ».

**Divergence** — les valeurs explosent (vitesses infinies ou démesurées) : le calcul est
arrêté avec des pistes. Voir [`depannage.md`](depannage.md).

**Nombre de Courant (CFL)** — distance parcourue par le fluide pendant un pas de temps,
rapportée à la taille d'une cellule. Les schémas explicites exigent un nombre de Courant
petit (`max_co`, `cfl`).

**Pseudo-transitoire** — calcul stationnaire mené comme un instationnaire avec un pas de
pseudo-temps, souvent plus robuste au démarrage (`pseudo_cfl`, `pseudo_dt`).

**Démarrage multigrille (FMG)** — on calcule d'abord sur des maillages plus grossiers, puis
on interpole sur le maillage fin : départ plus proche de la solution (`fmg_levels`).

**Solveur linéaire (direct, AMG)** — résout les systèmes d'équations de chaque itération :
direct (factorisation, petits maillages), multigrille algébrique (AMG, gros maillages).
`solver_p`, `solver_U` (`auto` choisit).

**Reprise (checkpoint)** — `checkpoint.npz`, écrit pendant et à la fin du calcul, permet de
poursuivre (« Continuer », `--continue`) ou de repartir sur un autre maillage
(`--restart`, champs interpolés).

## Fichiers

**TOML** — format texte du fichier de cas (`[section]`, `clé = valeur`). L'interface
l'affiche dans l'onglet « Fichier de cas (TOML) ».

**VTK** — `fields.vtk` : champs du calcul, à ouvrir dans ParaView pour des visualisations
avancées. Binaire par défaut (valeurs exactes) ; `[output] vtk_format = "ascii"` pour un
fichier texte lisible dans un éditeur.

**summary.json** — toutes les valeurs du résumé de fin de calcul (efforts, convergence,
sondes…), lisibles par un script.
