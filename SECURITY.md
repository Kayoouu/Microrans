# Sécurité

## Signaler une faille

Ne pas ouvrir d'issue publique. Utiliser **Security → Report a vulnerability** (signalement
privé GitHub) sur ce dépôt. Réponse visée sous 14 jours.

## Versions suivies

Seule la dernière version publiée (onglet Releases) et la branche principale reçoivent des
correctifs.

## Ce que fait et ne fait pas microrans

- **Aucun accès réseau, aucune télémétrie.** Le programme ne se connecte à rien ; seul le
  menu Aide → Documentation ouvre la page GitHub dans le navigateur, sur clic.
- **Fichiers de cas (`.toml`) venus d'ailleurs** : les formules en x, y (profils, valeurs
  initiales, sources) sont évaluées par un analyseur à liste blanche (`microrans/safe_expr.py`
  : nombres, x, y, opérateurs, fonctions mathématiques) — pas par `eval`, qui se contourne.
  Un cas ne peut donc pas exécuter de code. En revanche un cas choisit **où** les résultats
  sont écrits (`[output] directory`) et quels fichiers de maillage sont lus : comme pour tout
  fichier reçu, regarder ces chemins avant de lancer un cas inconnu.
- **Fichiers de reprise (`checkpoint.npz`)** : lus sans `pickle` (`allow_pickle=False`) ;
  un fichier piégé ne peut pas exécuter de code.
- **Exécutables** : construits par GitHub Actions à partir de ce code (journal public des
  constructions), non signés numériquement. Ne les télécharger que depuis l'onglet Releases
  ou les artefacts des workflows de ce dépôt.
- **Intégration continue** : jeton GitHub en lecture seule (écriture uniquement pour
  attacher les archives à une release), actions tierces épinglées par empreinte de commit,
  aucun secret utilisé.
