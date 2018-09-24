#BACH-SixthSense

###QUICK START
To score a set of .pdb structures:
  *make a text file containing the relative (with respect to the executable) or absolute path of each .pdb file you want to score. Put just one path per line
  *Type
```<FOLDERPATH>/BSS.x -COMPUTE_ENE -PDBLIST <listname>``` where <listname> is the relative (with respect to the executable) or absolute path of the list file you created
  - in the standard output the scoring of the structures will be printed ("ENERGY" rows). 
    Warnings about the structure and .pdb format will be printed as well ("WARNING" rows). 
    An output file output.bss with only the "ENERGY" rows will be created


###OPTIONS
  mandatory
  -COMPUTE\_ENE   computes the score of the structures in the structure list file
     <xor>
   COMPUTE\_PAR   learns the parameters on the structures in the structure list file. Do not
                 use unless you really know what you are doing
  -PDBLIST       specifies the relative or absolute path of the structure list file

  optional
  -FILE\_PAR      specifies an alternative relative or absolute path of the parameter file.
                 Default: ./BSS.par
  -o             specifies an alternative relative or absolute path of the output file.
                 Default: ./output.bss (take care not overwriting an existing output file!)
  -q             disables all the standard output messages, except "FATAL" errors
  -qwarn         disables all the standard output warnings ("WARNING" rows), but not "ENERGY"
                 rows or "FATAL" errors
  -NO\_SOLV       does not compute the solvation term. The code will run faster, yet the
                 score will be coarser
By browsing the code you can find that other options are implemented. We are still testing them.


###EXAMPLES WALKTHROUGH
1)  A CASP (monomer) decoy set and a CAPRI (dimer) decoy sets are reported as an example.
    You can find the already generated structures list files in the lists/ folder.
    To run the program over these decoy sets, just go in the executable folder and type

  ./BSS.x -COMPUTE\_ENE -PDBLIST lists/CASP.list -o my\_CASP.bss;
  ./BSS.x -COMPUTE\_ENE -PDBLIST lists/CAPRI.list -o my\_CAPRI.bss

    then compare with results/CASP.bss and results/CAPRI.bss respectively.


2)  You can build a rank by sorting the energy scores and see at what place the native
    is ranked. If the native has the lowest score of all, it will be ranked as first.
    The native is called "T0628.pdb" in the CASP decoy set and "native.pdb" in the
    CAPRI decoy set. Try and type
  
    for CASP:
  awk 'NR>1' my\_CASP.bss | sort -nk4 | awk 'index($2,"T0628.pdb") {n=NR} END {print "RANK (absolute, normalized, # of structures): ",n,n/NR,NR}'

    for CAPRI:
  awk 'NR>1' my\_CAPRI.bss | sort -nk4 | awk 'index($2,"native.pdb") {n=NR} END {print "RANK (absolute, normalized, # of structures): ",n,n/NR,NR}'

 
3)  You will notice that not all the .pdb structures in the CAPRI decoy set have the same
    number of residues. To overcome this, a CAPRI_R15_T36/cut_decoys/ directory is included.
    There you can find all the decoys structures (now called "cut_decoy_####") paired with their 
    corresponding natives ("cut_nat_decoy_####"). Each couple has the same number of residues,
    and can then be safely compared.
   
    The files lists/CAPRI_cut_d.list and lists/CAPRI_cut_n.list already contain the lists of the
    cut decoys and the corresponding cut natives, respectively. Try and type

  ./BSS.x -COMPUTE\_ENE -PDBLIST lists/CAPRI\_cut\_n.list -o my\_CAPRI\_cut\_n.bss;
  ./BSS.x -COMPUTE\_ENE -PDBLIST lists/CAPRI\_cut\_d.list -o my\_CAPRI\_cut\_d.bss

    To get a score, you will have to compare the score of each structure with the score of the
    corresponding native, and count how many times the score of the decoy is lower than the one
    of the native. If it happens n times, the native will be at place n+1 (so that if the native
    has always the lowest score, the rank will be 1).

    Try and type
  paste  my\_CAPRI\_cut\_n.bss my\_CAPRI\_cut\_d.bss | awk 'NR>1' | awk 'BEGIN{n=0} $8<\$4 {n++} END {print "RANK (absolute, normalized, \# of structures): ",n+1,(n+1)/(NR+1),NR+1}'


###FURTHER INSTRUCTIONS
  Moving the executable to another folder
    If you want to change the location of the executable, just remember to move ATOMIC_PARAMETERS_BSS
    and BSS.par files in the same folder. You can use -FILE_PAR option to specify an alternative
    folder for the .par parameter file, but you cannot do the same for ATOMIC\_PARAMETERS\_BSS.

  BACH-SixthSense reading method
    The program will read and score ONLY residues with the correct number of heavy atoms. Incomplete
    residues will be discarded, and a "WARNING" line will be generated before the corresponding
    structure "ENERGY" line in the standard output. You can disable these warning by means of
    -q or -qwarn options.

  -READ_MODE and -PREFIX options (beta!)
    By putting -READ_MODE instead of -COMPUTE_ENE (or -COMPUTE_PAR), the program will only read the
    .pdb structures and will output a .pdb file for each .pdb file read. The new .pdb files will
    reflect exactly what the program keeps into account: there won't be any incomplete residue,
    and the residue index will be strictly consecutive.
    WARNING: -READ_MODE works only if the structures to be modified are in the same folder of the
    executable. The output .pdb structures will be named by default "mod_<sctructurename>", where
    <structurename> is the name of the input .pdb files. If you want to change the prefix, use the
    -PREFIX option. In this way, by typing

  ./BSS.x -READ\_MODE -PDBLIST list -PREFIX ../folder/modif\_    (example)

    to the structure "name.pdb" mentioned in the structure list file "list" will correspond a structure
    "../folder/modif_name.pdb".
