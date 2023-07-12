Description:

This function takes multi-FASTA files, a list of motifs as well as the GenBank file of all the sequences of interest. Using these motifs, the sequences of each multi-FASTA file are scanned for any occurrence of the motifs and grouped. Additionally, users can provide a value to represent the minimum number of motif occurrences required for a sequence to be included in that motif’s group. In summary, if a motif is found in a sequence, that sequence is grouped with that motif. Otherwise, sequences containing none of the motifs from the input list are noted and added to a group containing all such members. If a sequence is found to have two or more motifs of the input list, this overlap is noted and all individuals sharing the same overlapped motifs are grouped. 
The recorded positions are then used to generate output.

- summary of motifs (count, distribution, gene location)
- bar graphs displaying positional information  
- Location of motifs both in genes and non-coding regions using GenBank files



Functions:

fastaFormatFileInput
- Opens FASTA file and stores them in headerAndSequence

genBankFileLink
- this function links .gb files with input FASTA sequences using the sequence itself as identifier


unique
- useful to change set() to list[] while keeping variable order




How to Run:

1) In the section for REQUIRED user input variables (line 165) all the variables required for motifSurvey function are listed and commented. 

2) In the section for OPTIONAL user input variables (line 141) edit these variables if needed. Otherwise the default values will be used. 

3) Once all variables are entered, simply run.



Overview of Structure:



1) fastaFormatFileInput(path) is the first function called. Users provide a path to their FASTA file(s) and this function separates all the contained sequences into the headerAndSequence dictionary. This function will determine if the path is a file or directory and act accordingly



2) genBankFileLink is the second function requiring user input. Users input the path to their GenBank files and this function automatically pairs file with their partner sequence if they are present.



4) Following this, headerAndSequence dictionary is searched and used to generate all the files of this analysis
