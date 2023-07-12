Description:

In the context of DNA, the sliding window method is an algorithm that uses a chosen window size (k) to iteratively discover all motifs of length k in each sequence. Once the motif in a window is examined, the window “slides” one position right allowing a new window with a different nucleotide at the end to be examined. Within motifSurvey, this algorithm allows the entirety of each input sequence to be scanned for motifs matching the chosen window size. Each motif matching window size (k) is selected iteratively and compared to the rest of the genome to see if it is repeated. Following each iteration, the window slides one nucleotide to allow a new motif to be observed. In the event an exact match is found downstream, the motif and positions of occurrence are recorded for further analysis. These matches can be filtered so that only motifs matching a user input pattern will be recorded. This user input filter pattern allows for regex notation in the event users require complex patterns to define their motifs of interest matching window size. Once all sequences are searched, the motifs found in each sequence are compared and those that occur in more than one of the input sequences are recorded. Using the recorded motifs, a series of output files are generated.


- summary of motifs found (count, distribution, gene location)
- bar graphs displaying positional information  
- an alignment of the recorded motifs using Clustal Omega. 
- Location of motifs both in genes and non-coding regions using GenBank files

Main Functions:


inputSeqCapsule

- Used to create object to store information about each sequence


fastaFormatFileInput

- Opens FASTA file and stores them in headerAndSequence



sequenceInformation

- Adds additional information to headerAndSequence



genBankFileLink

- this function links .gb files with input FASTA sequences using the sequence itself as identifier


genBankLocationObj

- Iterates through .gb files of each object to find where positions occur



seqSearch

- key function that calls the seqSearch function in each object to search itself for user input parameters. Generates all output files of motifSurvey




rangeLadder


- Used to make note of areas where substrings appear, bins are decided by user.





How to Run:

1) In the section for REQUIRED user input variable (line 1538) all the variables required for motifSurvey function are listed and commented. 

2) In the OPTIONAL user varaibles section (line 1515) please edit any paramaters that may be helpful. Otherwise, the default values will be used.


3) Once all variables are entered, simply run.


Overview of Structure:


1) fastaFormatFileInput(path) is the first function called. Users provide a path to their FASTA files and this function separates all the contained sequences into the headerAndSequence dictionary. 

2) sequenceInformation() adds additional information such as GC content to the value of all headerAndSequence entries.


3) masterSequenceDictionary is created. This uses the class creation function and an output from headerAndSequence (inputSeqMasterList) to give each entry in headerAndSequence dict its own object. In this object each sequence can store its own variables. This dictionary storing the object of each entry is termed 'masterSequenceDictionary'



4) genBankFileLink is the second function requiring user input. Users input the path to their GenBank files and this function automatically pairs file with their partner sequence if they are present.



5) seqSearch is called. This function takes masterSequenceDictionary, user input window (k), max sequence length, custom file prefix, bin size and an output directory from the user. Following this, the seqSerach function of each object is called. This is the last function required for file output.