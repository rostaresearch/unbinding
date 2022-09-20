Ligand Unbinding
--------------------
![cover](banner.jpeg)

This is an iterative ligand-protein unbinding repository
for the method detailed in:
[*J. Chem. Theory Comput.* 2022, 18, 4, 2543–2555](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00924)

*Although this method has been successfully applied to a
list of cases by the developers, it has not been thoroughly
tested. Feel free to report problems and unexpected behaviour
in the issues.*

### Dependencies

##### Python
Check `requirements.txt` 

##### VMD
The optional output for visualisation works with version 1.9.1.  
pbctools package needed (works with 2.6 version)  

---

### Walk-through

##### A. Initial unbias MD simulation
1. Generate the structure from pdb, homology model, etc  

2. Set up your system with your favorite method (for example charmm-gui or amber)

3. Run *x* ns (for example 20 ns, 10000000 steps) of free MD simulations  
In this project we are using NAMD for the MD simulations  
Folder name (for future reference): traj_0  

##### B. Distances Selection Process / Unbinding initial trajectory

4. Unbinding simulation iteratively
The script analyses the trajectory obtained from the previous step and calculate distances
to be included or excluded in the CVs.
    ```bash
    python main.py
    ``` 
Options:

| Flag             | Argument  |           Default |
|------------------|:----------|------------------:|
| --step           | enumerate |                   |
| --stride      -s | int       |                 5 |
| --lig         -l | string    |               LIG |
| --cutoff      -c | float(Å)  |               3.5 |
| --maxdist     -m | float(Å)  |               9.0 |
| --trajectory  -t | int       | "last trajectory" |
| --cumulative     |           |                   |
| -ns              | int(ns)   |                10 |
| --writeDCD       | Boolean   |             False |
| --processonly -p | Boolean   |             False |
| --nosave         | Boolean   |             False |

### step
Define the subprocess of the unbinding step, only for debugging. 

### stride (default=5)
Define how many frame to skip when analysing the DCD trajectory

### lig (default=LIG)
Define the ligand resname, correspond in the resname present in the psf/prmtop/pdb file

### cutoff (default=3.5)
Initial cutoff for identifying neighbours in Å. The script analyses the last trajectory and
if in more than 50% of the frames the distance appears, that distance will be included in
the new set of CVs

### maxdist (default=9.0)
Cutoff for excluding contacts in Å.

### trajectory
Which trajectory to analyze

### cumulative 
Rerun for all the trajectories till the one defined in the *trajectory* argument

### ns (default=10)
Length of the next simulation 

### writeDCD
Write a new DCD strided

### processonly  (default=False)
Do not write VMD and NAMD input. Other outputs will be written

### nosave  (default=False)
Do not save the checkpoint. For debug only


## Output generated:
Main folder
- **unbinding.out**
    Provides all the distances used as sum for the CVs used in the trajectories, and a summary
    of the distances excluded for the next trajectory
    
- **distances_traked.csv**
    Summary file tracking all the distances included and excluded

- **.checkpoint**
    Binary file with all the information stored
    
New **traj_(last)** folder:
- input file for NAMD simulation
- **sum_(last).col** Colvar file for NAMD sumulation  
 
5. Run the simulation
  
6. Repeat from 4-5 until the ligand is completely outside  

##### C. String Setup

The option `-step string` attempts to create colvar files
for a subsequent finite temperature string calculation.
*This function is experimental.*