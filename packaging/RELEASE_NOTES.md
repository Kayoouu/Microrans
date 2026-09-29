## Télécharger et lancer (aucune installation nécessaire)

**1. Choisir le fichier de votre ordinateur** (plus bas, rubrique « Assets ») :

| Votre ordinateur | Fichier à télécharger |
|---|---|
| Windows | `microrans-Windows.zip` |
| Linux | `microrans-Linux.zip` |

Pas de version Mac à télécharger : sur Mac, installer la version Python (voir le README,
section « Installation Python »).

**2. Décompresser** l'archive (clic droit → « Extraire tout… » sous Windows). Le dossier
`microrans` fait environ 350 Mo : c'est normal (tout ce qu'il faut est dedans).

**3. Lancer `microrans-gui`** (`microrans-gui.exe` sous Windows) par un double-clic.
(`microrans` sans « -gui » est la version ligne de commande pour utilisateurs avancés ;
double-cliqué sans rien taper, il ouvre lui aussi l'interface.)

Le programme n'est pas signé numériquement (cela demande un certificat payant) :
- **Windows** affiche « Windows a protégé votre ordinateur » → cliquer
  **« Informations complémentaires »** puis **« Exécuter quand même »** ;
- **Linux** : `./microrans-gui` dans un terminal ouvert dans le dossier.

## Premier essai (2 minutes)

1. Page **Accueil** : double-cliquer sur l'exemple `cavite_re100`.
2. Page **2. Maillage** : cliquer **« Générer le maillage »**.
3. Page **6. Calcul** : cliquer **« Lancer le calcul »** (quelques secondes ; les courbes
   de convergence s'affichent en direct).
4. Page **7. Résultats** : choisir une grandeur (vitesse, pression…) puis
   **« Tracer le champ »**. Les fichiers sont écrits dans le dossier `microrans_resultats`
   de votre dossier personnel.

Autres exemples : cylindre (écoulement autour d'un obstacle, lâcher de tourbillons),
plaque plane et profil d'aile turbulents, polaire d'un profil, sphère et tuyaux
(axisymétrique), convection naturelle (air chauffé), canal 1D.

## Contenu de cette version

- Simulation d'écoulements 1D, 2D plans et **axisymétriques** (tuyaux, jets, corps de
  révolution), laminaires ou turbulents (modèles Spalart-Allmaras, k-ε, k-ω, k-ω SST),
  stationnaires ou instationnaires, avec thermique et lois de paroi ;
- **polaires** Cl(α), Cd(α), Cm(α) et balayage de n'importe quel paramètre ;
- **sauvegarde automatique et reprise** : « Continuer le calcul précédent », ou démarrer un
  maillage fin depuis un calcul grossier ;
- mailleur intégré (rectangles, maillages autour d'objets, triangles, hybride) et
  import/export Gmsh, SU2, VTK, OpenFOAM ;
- interface graphique et ligne de commande (`microrans` dans le même dossier).

Méthodes, résultats de validation mesurés et **limites connues** : voir le README du dépôt.
