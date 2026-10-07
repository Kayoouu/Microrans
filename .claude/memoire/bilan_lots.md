# Bilan des lots : estimé / réel, rétrospectives, indicateur du jalon A

À compléter à la fin de chaque lot (skill `lot`, étape « Fin »). Sert à :
1. comparer estimé et réel, pour que les délais donnés à l'utilisateur reposent sur des
   chiffres (jusqu'ici : estimations seulement) ;
2. suivre l'indicateur du jalon A (résultats faux silencieux trouvés par audit) ;
3. garder les changements de méthode issus des rétrospectives (et vérifier qu'ils servent).

## Lots

Réel = durée de session mesurée par l'heure de la relance et des commits (git log), pas le
temps de calcul de l'utilisateur. « Défauts après coup » : défauts introduits par ce lot
et trouvés plus tard (à remplir rétroactivement).

| Date | Lot | Estimé | Réel | Tests ajoutés | Défauts après coup | Rétro : ce qui a coûté → changement |
|---|---|---|---|---|---|---|
| 2026-10-04 | F1 (9 points) | 1 lot | ≈ 50 min (relance 03:53 → 04:41) | 11 | — | C15 bloqué sur une décision de l'utilisateur → poser la question tôt, continuer le reste |
| 2026-10-04 | C15 | 1 lot | ? (début non noté) | 3 | — | — |
| 2026-10-05 | F1b (5 points) | 1 à 2 lots | ≈ 45 min (03:53 → 04:35) | 6 | — | `git stash` à la main pour prouver chaque test, comparaisons avant / après réécrites à chaque fois, suite lancée pendant l'édition de la doc (course connue) → outils tools/dev (echoue_avant, ab, suite) |
| 2026-10-05 | Outillage (demande de l'utilisateur) | non estimé | ≈ 35 min (≈ 23:05 → CI verte 23:37) | 4 | — | 3 défauts des outils trouvés en les essayant (--ref avalé, erreurs déclarées identiques, « git push » cité pris pour une poussée) → essayer chaque outil sur un cas connu avant de s'y fier |
| 2026-10-06 | F2 (12 points) | 1 lot (F2 à F4 : 3 à 5) | 34 min (03:53 → CI verte 04:27 ; 2 commits, 2 suites, 2 constructions) | 12 | — | cas de reproduction mal écrits deux fois (compressible sans `compressible = true` ni état amont ; 3D sans faces d'extrusion) → partir d'un exemple existant modifié plutôt qu'écrire un cas à la main. Outils : echoue_avant 12 preuves en 1 appel, ab egalite 5 exemples en < 1 min : utiles |
| 2026-10-07 | F3 (7 points) | 1 à 2 lots | 55 min (03:53 → commit 04:48 ; 3 suites dont 2 en trop, 10 min chacune) | 8 (+1 outil) | — | suite lancée trop tôt et en-tête CSV figé dans un ancien test → skill § 3.1 (suite sur le code fini, grep de l'ancien format). Outils : echoue_avant a signalé 2 preuves faibles, corrigées ; ab.py amélioré (colonnes ajoutées) |

## Indicateur du jalon A : résultats faux silencieux trouvés par audit

Définition : le calcul ou l'affichage se termine normalement, sans message, et le résultat
est faux ou trompeur. Les plantages et les messages obscurs ne comptent pas. Jalon A
atteint quand un audit complet en trouve 0. Les audits n'ont pas eu la même profondeur :
la tendance est indicative.

| Audit | Nombre | Points |
|---|---|---|
| 1 (avant D4) | 8 | C1, C4, C5, C6, C7, C8, C9, C10 |
| 2 | 3 | C13, U14, U17 |
| 2 approfondi | 2 (+1 limite) | C15, C21 ; limite : C17 (reprise annoncée exacte, ne l'était pas) |

Trouvés hors audit, en travaillant : 2026-10-07 (F3, en mesurant P2) figure de convergence
de la cavité illisible depuis toujours (résidu nul tracé à 1e-300 : axe sur 330 décades),
interface et ligne de commande. Affichage trompeur sans message : l'audit 3 doit regarder
les figures produites, pas seulement les chiffres.

## Campagnes de non-régression (tools/audit) : dernier relevé

Relancer le lundi (ou après 7 jours) ; écart avec la ligne ci-dessous = régression à
expliquer avant le lot. Sorties en dehors du dépôt (scratchpad), seulement le verdict ici.

| Campagne | Date | Résultat |
|---|---|---|
| c1 exemples en ligne de commande | 2026-10-04 (F1) | 25 / 25 |
| c6 interface = ligne de commande | 2026-10-04 (C15) | 22 exemples identiques au bit près (hors chemins), Sod ignoré |
| c7 reprise exacte | 2026-10-05 (F1b) | 18 / 18 au bit près, Sod non découpable |
| c8 2D = 3D une couche | 2026-10-03 (audit 2 approfondi) | convergé : écart ≤ 2·10⁻⁶ ; perturbation 3D corrigée depuis (C16) |
| c3b / c3c clés numériques (texte, −1) | 2026-10-06 (F2) | 92 / 92 et 24 / 24 refusés clairement avant le calcul |

## Changements de méthode (rétrospectives) et effet constaté

- 2026-10-05 : outils `tools/dev` (suite sur copie, echoue_avant, ab egalite / temps, ci,
  avant_push), hook avant poussée, hook de fin de tour (travail non poussé), état de la CI
  au démarrage, skill `lot`. Effet constaté au lot F2 (2026-10-06) : preuve des 12 tests en
  un appel (avant : un `git stash` par test) ; comparaisons avant / après sur 5 exemples
  sans script à écrire ; état de la CI lu au démarrage. Pas encore d'erreur évitée par le
  hook avant poussée.
- 2026-10-06 (rétro F2) : heure de fin écrite de tête (04:48, réel 04:12) → la procédure
  (skill `lot`, § 5) impose de lire les heures (git log, date -u). Effet constaté au lot F3 (2026-10-07) : heures lues (date -u) à chaque étape.
- 2026-10-07 (rétro F3) : suite complète lancée sur un arbre provisoire (à refaire) puis
  échec d'un ancien test qui figeait l'en-tête de wall_*.csv (vu après 10 min) : ≈ 20 min
  de suites en trop → skill `lot` § 3.1 : suite une fois le code fini ; format de sortie
  changé → `grep` de l'ancien texte dans tests/ avant. Preuves faibles évitées : import des
  modules nouveaux en fin de test (§ 2.2).
