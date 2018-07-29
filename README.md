## VASP-toolkit

​	VASP-toolkit is a set of tools for analyzing VASP (The Vienna Ab initio simulation package) calculation result files (EIGENVAL, DOSCAR etc.). Every functional module (xxxAnalyst.py) in company with internal tool module (internalTools.py) has independent function.

#### Usage

1. Band structure (EIGENVAL file) analysis

   _bandAnalyst.py_ and _internalTools.py_ are needed for band structure analysis.

   `import bandAnalyst`

   `bandAnalyst.parseEIGENVAL(EIGENVAL, rl)`

   where _EIGENVAL_ is the EIGENVAL file name (include path), _rl_ is reciprocal lattice vectors of optimized structure which can get from OUTCAR file.

2. Density of States (DOSCAR file) analysis

   _dosAnalyst.py_ and _internalTools.py_ are needed for DOS analysis.

   `import dosAnalyst`

   `dosAnalyst.parseDOSCAR(DOSCAR, atoms)`

   where _DOSCAR_ is the DOSCAR file name (include path), _atoms_ is the atom No. collection whose LDOS (Local DOS) you want to analyse.

#### Dependencies

​	Pure Python and No dependencies.

#### Author

​	E-mail: zjwwho@gmail.com
	
​	Github: https://github.com/zjwwho/VASP-toolkit

​	Please contact me with any questions.