# openmm-scripts
Helpful python scripts written to be able to easily submit an MD trajectory submitted to either a GPU or CPU on a cluster or local pc. 
Andrew Powers (apowers4@uoregon.edu OR apowers@flatironinstitute.org)


## Repo Setup
---
- `md_complex.py`: This is for running a target + peptide MD trajectory. This can be done for some specified nubmer of steps or a production run of 100ns.
- `md_monomer.py`: This is for running a monomer MD trajectory. This can be done for some specified nubmer of steps or a production run of 100ns.

### Flags
---
```
p.add_argument("--file", type=str, help="Path to Complex (target + cyclic peptide).")
p.add_argument("--output-dir", type=str, default="out/", help="Directory to store output files. Will generate dir, if not yet created.")
p.add_argument("--temperature", type=float, default=300.0, help="Temperature in Kelvin")
p.add_argument("--pressure", type=float, default=1.0, help="Atmospheric Pressure in atm")
p.add_argument("--timestep", type=float, default=2.0, help="Time step in femtoseconds")
p.add_argument("--total-steps", type=int, default=30_000, help="Total number of steps to take in Simulation")
p.add_argument("--postfix-filename", type=str, default="out", help="Specify this to make your output .pdb and .dcd unique")
p.add_argument("--add-hydrogens", action="store_true", help="Add hydrogens (Good to use if you pass --pdbfix)")
p.add_argument("--production", action="store_true", help="Specify if you want to run a production 100ns run, if not only PS and total-steps are run.")
p.add_argument("--debug", action="store_true", help="Specify if you want output to STDOUT during Production run for testing.")
p.add_argument("--cyclic", action="store_true", help="Specify if you want the peptide to be cyclic.")
p.add_argument("--pdbfix", action="store_true", help="Perform PDB cleanup first.")
```


### Version
---
`0.0.1`
