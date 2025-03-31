# Simulation de Particules avec Barnes-Hut

Ce projet simule les interactions gravitationnelles entre particules en utilisant l’algorithme de Barnes-Hut. Voici les instructions pour compiler et exécuter la simulation.

## Compilation

1. Créer le répertoire de compilation :
```bash
    mkdir build
```
2. Compilation avec scripts utilitaires

Au lieu d'exécuter manuellement les commandes de compilation, nous fournissons quatre scripts utilitaires permettant de compiler facilement chaque version :

* MPI + OpenMP (scripts/compile_mpi_omp.sh)

* MPI (scripts/compile_mpi.sh)

* OpenMP (scripts/compile_omp.sh)

* Séquentiel (scripts/compile_sequential.sh)

Chaque script accepte un paramètre optionnel "sdl" pour activer la visualisation SDL (nécessite les paquets libsdl2-dev et libsdl2-ttf-dev sur les systèmes Debian).

## Exécution de la simulation

Le binaire de la simulation se trouve dans le répertoire build. Vous pouvez exécuter la version séquentielle de la simulation avec :
```bash
./build/particleSimulation <nombre_de_particules> <nombre_d_itérations> [sdl]
```
De la même manière, nous fournissons des scripts d’exécution pour simplifier le lancement des versions OpenMP, MPI et MPI+OpenMP :

* MPI + OpenMP

Utilisation : `./run_mpi_omp.sh <nombre_de_processus> <nombre_de_threads> <nombre_de_particules> <nombre_d_itérations> <version_mpi> <build_dir> [sdl]`

* MPI

Utilisation : `./run_mpi.sh <nombre_de_processus> <nombre_de_threads> <nombre_de_particules> <nombre_d_itérations> <version_mpi> <build_dir> [sdl]`

* OpenMP

Utilisation : `./run_omp.sh <nombre_de_threads> <nombre_de_particules> <nombre_d_itérations> <build_dir>`

## Profilage

Pour générer un profil visuel avec gprof (nécessite la compilation avec -pg) :
```bash
gprof ./build/particleSimulation | gprof2dot | dot -Tpng -o output.png
```
Cela créera un fichier output.png contenant les informations de profilage.
### Utilisation de MAQAO

Pour le profilage MPI :
```bash
OMP_NUM_THREADS=<nombre_de_threads> maqao oneview -R1 --number-processes=<nombre_de_processus> --mpi-command="mpirun -n <nombre_de_processus>" --number-processes-per-node=<process_par_noeud> xp=ov_mpi --replace -- ./build/particleSimulationMPI <nombre_de_particules> <nombre_d_itérations>
```
Pour le profilage OpenMP :
```bash
OMP_NUM_THREADS=<nombre_de_threads> maqao oneview -R1 xp=ov_omp --replace -- ./build/particleSimulation <nombre_de_particules> <nombre_d_itérations>
```
Pour le profilage séquentiel :
```bash
maqao oneview -R1 xp=ov_orig --replace -- ./build/particleSimulation <nombre_de_particules> <nombre_d_itérations>
```