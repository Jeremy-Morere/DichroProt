# DichroProt
Dichroprot generates theoretical near-UV circular dichroism spectra for proteins. It splits the protein into photo-active residues and starts a QM-MM calculation for each of them using amber16 and gaussian09.
The program couples the gaussian outputs using two different methods (a reference will be added after release). 

## Prerequisite
The current version of DichroProt has been developed and tested on _Centos Linux 7_:

* MKL library version 2019.4
* Intel version 2019.4 ( Fortran compiler)
* Amber version 16
* Gaussian version 09

Paths to mkl libraries, Fortran compiler, Gaussian and Amber must be add in the installation file.

To compile ***DichrProt*** use the install bash file  
`bash install`.

The DichroProt repository contains all the programs and executables necessary to obtain a near-UV CD spectrum.
+ dichroprot_prog: Coupled Gaussian 09 calculations; produces a CD spectrum with or without coupling.
+ geometry_prog: Analyses the protein structure from a .pdb; generates a geometry.dat file which is used to obtain the CD spectrum.
+ other_prog: Other small post-processing programs.
+ executable: Contains all the executables needed.
+ test_snapshot: Test repository (Lysozyme (PDB:2YVB) is given as an example).

The executable must be added to the PATH. This can be done by adding to the .bashrc file  
`export PATH=/path_to_executable:$PATH`

## Tutorial

### Snapshot Tutorial
> [!IMPORTANT]
> For a specific conformation, a topology file (.prmtop) and a geometry file (.pdb) are required.

1) Enter your working repository containing structure and topology files (.prmtop and .pdb).
2) Run `create_input` and adapt the GlobalInput file to your protein:
- PDB : protein.pdb
- PRMTOP : protein.prmtop
- Level of theory : B3LYP/6-311+G**
- Number of state : 5
- Full contribution : n

> [!TIP]
> Nearly all transitions in the near-UV are covered for 5 calculated states.

> [!TIP]
> If the full contribution is set to yes, the computational effort is multiplied by the number of states, for a very small correction of the CD spectra.You can get a first impression without the full contribution on.

3) Run `prep_qmmm GlobalInput`.
4) Adapt submit.sh to your architecture and gau_job.tpl if you need to add specific keywords for the Gaussian computation.
5) Run `launch_qmmm GlobalInput`.
6) Run `prep_dichroprot GlobalInput`, all output files will be collected in the DichroProt repository.
7) Go to DichroProt and customize the DichroProtInput:
- perturbation: Coupled Gaussian output with a perturbation approach
- Matrix: Coupled Gaussian output with an effectif Hamiltonian
- Superposition: Each coupling is perfermed, the resulting spectrum is the superposition of the residual spectra.
  
> [!TIP]
> The computational cost of this step is negligible, all calculations can be done at the same time.
9) Run 'dichroprot DichroProtInput'.

### Scan Tutorial

> [!IMPORTANT]
> For a particular conformation, a topology file (.prmtop) and a trajectory file (.dcd or .nc) are required.

1) Enter your working repository containing the trajectory and topology files (.prmtop and .dcd/.nc).
2) Run `create_input` and adapt the GlobalInput for QMMM calculations and the DichroProt coupling following the advice in steps 2) and 7) of the Snapshot tutorial.
3) Run `prep_scan GlobalInput`.
4) Adapt submith.sh to your architecture and gau_job.tpl if you need to add specific keywords for the Gaussian computation.
5) Run `launch_scan GlobalInput` (it will start every qmmm calculation)

>[!WARNING]
> This step is very expensive. If you are running the program on a shared cluster, we recommend that you do it during a low traffic period (e.g. at night), or start the calculation in batch following the snapshot tutorial.

6) When all Qmmm calculations are finished, check if any crashed and restart them following the snapshot tutorial.
7) Run `launch_dichroprot GlobalInput`.
8) Spectra of all snapshots will be collected into the spectrum repository and the average spectra will be calculated.
