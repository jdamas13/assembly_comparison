# assembly_comparison
Identifies differences between two assemblies that were aligned using [Mashmap2](https://github.com/marbl/MashMap).

Reported results are number of changes in query genome when cmpared to reference.   

E.g., for comparison of intermediate assemblies from an the assembly pipeline, reference genome should be more fragmented assembly (i.e., contig assembly) and query genome the more contigous assembly (i.e., scaffold assembly). Then, reported results will be number of changes introduced by the scaffolding step. 

## Usage:
perl assembly_comparison.pl \<mashmap.out\> \<reference lengths\> \<resolution\>

## Input:
- \<mashmap.out\>

Output file from Mashmap2 alignment.
   
- \<reference lenghts\>

Lenghts of reference chromosomes, scaffolds or contigs.
Format: ID Length

- \<resolution\>

Minimum length of blocks to be included, in base pairs. 
Usually same as resolution used for Mashmap2 alignment (s parameter).

## Output:
- \<mashmap.out\>.diff

List of end to end and inconsistent adjacencies between reference and query assemblies.

- \<mashmap.out\>.diff.broken

List of breaks introduced in reference where ends are kept free in query.

- \<mashmap.out\>.diff.summary

Summary of number of differences found between reference and query assemblies.
Reports number of end to end joins, number of inconsistent joins, no of scaffolds involved in free-end breaks and number of breaks where ends are kept free.

![image](https://github.com/jdamas13/assembly_comparison/assets/30689196/88bbd50e-9457-42d5-97aa-1acf4e9c11da)

