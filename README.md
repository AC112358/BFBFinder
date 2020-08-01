# BFBFinder
Breakage Fusion Bridge detection algorithms

The BFB Finder project provides algorithms for detecting evidence for Breakage Fusion Bridge (BFB) genome rearrangement events given genomic segment copy numbers. In particular, given (exact or noisy) segment count vectors, it allows (a) deciding if the input may be produced by BFB, (b) output one or the subset of all BFB count vectors of the input (noisy counts only), and (c) output one or the set of all BFB architectures supporting the input.

The implemented algorithms were described in the following papers:

1. Algorithms for BFB detection in tumor genomes. Shay Zakov, Marcus Kinsella, and Vineet Bafna. Proceedings of the National Academy of Sciences Apr 2013, 110 (14) 5546-5551; DOI: 10.1073/pnas.1220977110 

2. Reconstructing breakage fusion bridge architectures using noisy copy numbers. Shay Zakov and Vineet Bafna. Journal of Computational Biology. Jun 2015.577-594. http://doi.org/10.1089/cmb.2014.0166

If you publish results obtained by using BFB Finder, we would appreciate if your publication contains citations of these paper. 


Installation:

1. Open a terminal window and set current directory to where you would like to install BFB Finder, e.g. “/home/my_username/projects”.
2. Run “git clone https://github.com/shay-zakov/BFBFinder”.
3. Change directory to “BFBFinder”.
4. Run “gradle build”
4. To check if build was successful, run “java -jar ./build/libs/BFBFinder-0.1.0.jar -s -a [4,4,4]”.
Expected output should be 

ABCCBAABCCBA

ABCCBBCCBAAA

ABCCCCBAABBA

Total number of strings: 3


Running:

Run “java -jar <path to BFB Finder jar>/BFBFinder-0.1.0.jar <argument list>”

Main arguments:

• Last argument should be a count vector of the format "[c1,c2,...,ck]" (c_i's are integers, no spaces).

• Modes of run are either "Decision", "Vector" (finding neighboring BFB vectors), or "String". If "-v" or "-s" arguments are supplied, the mode is set to "Vector" or "String" respectively, otherwise it is "Decision".

• If an "-a" argument is supplied, all solutions are returned, and otherwise a single solution is returned ("Decision" mode ignores this argument).

• To set error bounds, use "-e=<error model class name>" and "-w=<min weight>". Error model class names are NoErrorModel (default), PoissonErrorModel, MaxRelativeErrorModel, and CanberraErrorModel. The min weight value should be set between 0 (don't do that!) and 1. In short, a value of 1 allows no errors, while the closer the value is to 0 more errors are allowed. Good chances there are still some bugs in this mechanism.

Some examples of possible arguments:

    1. Valid BFB counts: "[3,3,4,6,34]"
    
    2. Invalid BFB counts: "[2,3,4,6,34]"
    
    3. Getting a single BFB string corresponding to the given counts: "-s [3,3,4,6,34]"
    
    4. Getting all BFB strings corresponding to the given counts: "-s -a [3,3,4,6,34]" 
    
    5. Getting all BFB counts with weight at least 0.858 with respect to the given counts, according to the Poisson error model: "-v -a -e=PoissonErrorModel -w=0.858 [4,3,4,6,34]"
    
    6. Getting all BFB strings corresponding to counts with weight at least 0.858 with respect to the given counts, according to the Poisson error model: "-s -a -e=PoissonErrorModel -w=0.858 [4,3,4,6,34]"
    


