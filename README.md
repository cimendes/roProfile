## roProfile ##

Generation of wgMLST allelic profile using [Roary's] (https://sanger-pathogens.github.io/Roary/) output.

This script will use [Roary's] (https://sanger-pathogens.github.io/Roary/) gene presence and absence file to generate a profile for all loci in the pan-genome, with loci absence indicated with the number 0. 
roProfile provides aditional options, allowing for the profile for the core (genes present in all samples) to be obtained, to transpose the gene presence absence rtab file from [Roary] (https://sanger-pathogens.github.io/Roary/) to be used as profile, and to obtain a frequency plot for the pan-genome.
roProfile will also save the fasta sequence files in a seperate directory, with all alleles for all loci, indicating the ones belonging to the accessory and core genomes.

## Updates:

3/11/2016 - version 1.4.0 - Fixed bug where the sequences in the reversed strand weren't being retrieved propperly. 

24/10/2016 - version 1.3.0 - Now it's possible to obtain a roary gene_presence_absence.csv file without the loci that were removed from the pan-genome profile. A log file containing the loci removed is now implemented.

21/10/2016 - version 1.2.0 - The loci with a size variation greater than mode+-(mode*threshold) are now removed. the threshold value is set to 0.2 by default, but can be altered. 

20/10/2016 - version 1.1.0 - Started process to remove problematic loci. The Loci with multiple alleles are now being removed from the profile. The removal of the loci with allele size too variable and generation of a log file for the removed loci are planned to be implemented next.

19/10/2016 - version 1.0.2 - fixed a bug where the full gene sequence wasn't being retrieved by one nucleotide

## Usage
    roProfile.py [-h] [-r ROARY] [-d GFFDIR] [-p] [-t] [-f] [--version]

    Generation of pan-genome profile files using Roary output. By default, it will
    generate a profile for the full pan-genome, with Locus Not Found represented as 0

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
roProfile  requires two input paths, the path for [Roary's] (https://sanger-pathogens.github.io/Roary/) output directory and the path to the directory containing the GFF files used in the Roary analysis, like the ones obtained with [prokka] (https://github.com/tseemann/prokka). 
From [Roary's] (https://sanger-pathogens.github.io/Roary/) output directory the script will use the gene presence and absence file, that can be seen bellow, then use that information to parse all the GFF files. Other files in the Roary output directory may be needed for aditional options. 

## Output
roProfile outputs the a wgMLST profile file, with missing data represented as '0'. The fasta sequence files, with all alleles for all loci, indicating the ones belonging to the accessory and core genomes, in a saved in a seperate directory named "sequences" within the working directory. These sequence files for the loci can then be analyzed with [chewBBACA's SchemaEvaluator] (https://github.com/mickaelsilva/chewBBACA/tree/master/SchemaEvaluator).
All loci (Roary's gene groups) with a multiple alleles (paralog genes that weren't split by Roary) are removed from the wgMLST profile, so I recomend not using this script if the option to not split paralogous genes in Roary was used (-s). The loci with a size variation greater than mode+-(mode*threshold) are also removed. This threshold is by default set to 0.2, but can be altered by the user with the '-threshold' flag.

#### The '-core' flag:
This option will produce a second profile file containing the cgMLST for the dataset, meaning only the genes present in all isolates after loci size and multiple allele correction. This file can be used in software to generate and analyze Minimum Spanning Trees (MST), like [PHYLOViZ Online](https://online.phyloviz.net/index)

#### The '-transpose' flag:
This option will use the existing gene_presence_absence.Rtab file in [Roary's] (https://sanger-pathogens.github.io/Roary/) output directory and transpose it so it can be used as an presence and absence profile file. 

#### The '-frequency' flag:
This option will generate an interactive pan-genome frequency plot in a html, allowing a quick overview of the pan-genome.
![pan-genome frequency](https://cloud.githubusercontent.com/assets/15690332/19932061/97ed2354-a106-11e6-8ddd-1f693d31e6b4.png)

#### The '-genefile' flag:
This option will generate a new [Roary's] (https://sanger-pathogens.github.io/Roary/) gene presence and absence csv file without the genes with a large size variation or with unsplit paralog genes (considered multiple alleles). 
![gene_presence_absence.csv output] (http://sanger-pathogens.github.io/Roary/images/gene_presence_and_absence.png)

## Coming soon
- Remove the loci with large size variation and multiple alleles from the presence and absence profile ('-transpose' flag)
- Remove the loci with large size variation and multiple alleles from the pan-genome frequency plot ('-frequency' flag)

## License
roProfile is freely available under a GPLv3 license.

## Acknowledgements
-Andrew Page, for providing valuable assistance in understanding Roary's behaviour

## Contact
Catarina Mendes (cimendes@medicina.ulisboa.pt)

