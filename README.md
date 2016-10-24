## roProfile ##

Generation of pan-genome profile files using [Roary] (https://sanger-pathogens.github.io/Roary/) output.

It will use [Roary's] (https://sanger-pathogens.github.io/Roary/) gene presence and absence file to generate a profile for all genes in the pan-genome. If a gene is absent, it's indicated with the allele number 0. 
roProfile also allows for the profile for the core (genes present in all samples) to be obtained and to transpose the gene presence absence rtab file from [Roary] (https://sanger-pathogens.github.io/Roary/) to be used as profile, indicating the pan-genome gene presence and absence with binarys. 
roProfile will produce the profile files requested and the fasta sequence files, with all alleles, for all loci, indicating the ones belonging to the accessory and core genomes. These sequence files can be analyzed with [chewBBACA's SchemaEvaluator] (https://github.com/mickaelsilva/chewBBACA/tree/master/SchemaEvaluator).

## Updates:

24/10/2016 - version 1.3.0 - Now it's possible to obtain a roary gene_presence_absence.csv file without the loci that were removed from the pan-genome profile. A log file containing the loci removed is now implemented.

21/10/2016 - version 1.2.0 - The loci with a size variation greater than mode+-(mode*threshold) are now removed. the threshold value is set to 0.2 by default, but can be altered. 

20/10/2016 - version 1.1.0 - Started process to remove problematic loci. The Loci with multiple alleles are now being removed from the profile. The removal of the loci with allele size too variable and generation of a log file for the removed loci are planned to be implemented next.

19/10/2016 - version 1.0.2 - fixed a bug where the full gene sequence wasn't being retrieved by one nucleotide

## Usage
    roProfile.py [-h] [-r ROARY] [-d GFFDIR] [-p] [-t] [-f] [--version]

    Generation of pan-genome profile files using Roary output. By default, it will
    generate a profile for the full pan-genome, with Locus Not Fund represented as 0

    optional arguments:
        -h, --help            show this help message and exit
        -r ROARY, --roary ROARY
                        Path to directory containing all output files from
                        Roary (https:/sanger-pathogens.github.io/Roary)
        -d GFFDIR, --gffdir GFFDIR
                        Path to directory containing all gff files used in the
                        Roary analysis.
        -c, --core      Generate profile file for the core-genome only, with
                        genes present in all isolates.
        -t, --transpose       transpose the gene presence absence rtab file from
                        roary to be used as profile
        -f, --frequency       Generate pan-genome frequency plot
        -th [THRESHOLD], --threshold [THRESHOLD]
                        Threshold for the allele size (default=0.2).
        -g, --genefile        Obtain a roary's gene presence and absence csv file
                        without the removed loci.
        --version             Display version, and exit.

## Dependencies

- Python (2.7.x)
- [Biopython] (http://biopython.org/) (1.66 or similar)
- [matplotlib] (http://matplotlib.org/index.html) (1.5.2 or similar)
- [mpld3] (http://mpld3.github.io/) (0.2 or similar)
- [pandas] (http://pandas.pydata.org/) (0.15.0 or similar) (if the option to transpose is used)

## Installation

roProfile is a standalone python script and does not require any installation. Simply clone the git repository:

    git clone https://github.com/cimendes/roProfile

## Input
roProfile  requires two input paths, the path for [Roary's] (https://sanger-pathogens.github.io/Roary/) output directory and the path to the directory containing the GFF files used in the Roary analysis.

## Output
roProfile outputs the requested profile files and the fasta sequence files, with all alleles for all loci, indicating the ones belonging to the accessory and core genomes, in a seperate directory.

## License
roProfile is freely available under a GPLv3 license.

## Contact
Catarina Mendes (cimendes@medicina.ulisboa.pt)

