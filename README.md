# Look4LTRs
Long Terminal Repeat Retrotransposon detection tool capable of finding recently nested LTR RTs de novo.

## Requirements

GNU g++ 11.1.0 or later

## How to Compile

```
mkdir bin
cd bin
(If your default compiler meets the version requiremet) 
cmake ..
(Or if you would like to specify a different compiler that meets the requirement)
cmake .. -DCMAKE_CXX_COMPILER=your_compiler_name_for_example_g++-7
(Or if this fails, try using this to set a different compiler. Replace the paths with your own.)
cmake .. -DCMAKE_CXX_COMPILER=$HOME/C++/GCC/bin/g++ -DCMAKE_C_COMPILER=$HOME/C++/GCC/bin/gcc -DCMAKE_PREFIX_PATH=$HOME/C++/GCC
make
```

## INPUT
Look4LTRs accepts FASTA format files, as well as multi-FASTA format. It is suggested that, at minimum, an entire genome is given to Look4LTRs, to enhance its self-supervised capabilities. Multiple genomes may be passed into Look4LTRs, but take caution with the memory requirements, and lower the number of threads (if any) if too much memory is being utilized. Additionally, Look4LTRs accepts training genomes that will not be predicted upon but are used to enhance the prediction of other genomes.

## OUTPUT
In the output directory, Look4LTRs will create three new directories. Look4LTRs will not overwrite these directories or files within UNLESS the files are generated from the same FASTA files. The directories created are as follows:

  * Bed - Files containing the start and ends of LTR RTs in BED format
  * Rtr - Files containing useful information about the LTR RTs such as the poly purine trail location, start and ends of the LTRs, target site duplication locations, etc.
  * Cpx - Files containing regions that Look4LTRs notes to be complicated; i.e., many elements of the same family may be residing in these areas and so any prediction made here is unsure.

### Bed
The typical BED format file consists of three columns.  These BED files only report the start and ends of LTR RTs, not solo LTRs.

  1. *chrom* is the chromosome identifier; in the case of multi-FASTA format files, Look4LTRs will have multiple chromosomes in the same BED file.
  2. *Start* is the start of the LTR RT. 
  3. *End* is the end of the LTR RT, exclusive.
  
### Rtr
RTR stands for RetroTransposon Relationship format file. It is Look4LTRs specific format that outputs information about detected LTR RTs and their relationship to other LTR RTs found by Look4LTRs. RTR format files also contain the locations of solo LTRs. RTR format has 16 columns.
  
  1. *chrom* is the chromosome identifier.
  2. *ID* is a unique identifier for the LTR RT per chromosome.
  3. *LeftStart* is the start of the LTR RT (start of the left LTR) or start of a solo LTR.
  4. *LeftEnd* is the end of the left LTR or end of a solo LTR.
  5. *RightStart* is the start of the right LTR or NA if a solo LTR.
  6. *RightEnd* is the end of the LTR RT (end of the right LTR) or NA if a solo LTR.
  7. *NestIn* is a curly-brace enclosed list of IDs of LTR RTs or solo LTRs that are nested in the LTR RT. If a solo LTR, this is NA.
  8. *NestOut* is an ID of the LTR RT that the current element is nested inside of. If there is no outer nest, then this is NA.
  9. *RC* is a boolean for Reverse Complement. If the LTR RT is reverse complemented, this is 1. If not, this is 0. If a solo LTR, this is NA.
  10. *PPTStart* is the start of the poly purine tract, or NA if a solo LTR.
  11. *PPTEnd* is the end of the poly purine tract, or NA if a solo LTR.
  12. *TSDStart* is the start of the target site duplication. This is NA if a TSD could not be found.
  13. *TSDEnd* is the end of the target site duplication. This is NA if a TSD could not be found.
  14. *CaseType* reflects which matching case this LTR RT came from; for example, a recently nested LTR RT would be marked with RecentlyNested.
  15. *GraphGroup* Look4LTRs utilizes a graph when matching LTRs. Each graph represents nearby elements of the same family. This column is an ID to the graph the LTR RT/solo LTR belongs to.
  16. *LTRIdentity* is the identity between the two LTRs of an LTR RT, NA if a solo LTR.

### Cpx
CPX stands for ComPleX format file. It is a Look4LTRs specific format that outputs complex regions of many same-family elements that Look4LTRs was unsure about. There are a variable number of columns.

  1. *chrom* is the chromosome identifier.
  2. *GraphID* is the ID of the graph that found this complex region.
  3. *Start* is the start of an element in this complex region
  4. *End* is the end of an element in this complex region
  5. *...* repeated *Start* and *End* columns for each element in the region.

## Parameters
Look4LTRs is activated from the command line. The following table describes the parameters as well as if they are required.

| Parameter | Description | Required? |
|-----------------|-----------------|-----------------|
| -f/--fasta | Fasta file directory for training and predicting. If you wish to train multiple genomes, you can pass multiple directories here. | Yes |
| -o/--output | Output directory. | Yes |
| -t/--train | Fasta file directory for training only. Can be given a variable number of arguments.
| -p/--parallel | Number of threads to use. If not given, defaults to 1 | No |
| -h/--help | Prints a help message and stops execution of the program | No |

## Usage

1. Predicting on a genome:
    ```bash
    ./look4ltrs --fasta /###/###/Sorghum_bicolor/Fasta/ --out /###/###/outputdir/ --parallel 8

2. Predicting on multiple genomes from different directories is this:
    ```bash
    ./look4ltrs --fasta /###/###/Phaseolus_vulgaris/Fasta/ /###/###/Vigna_radiata/Fasta/ /###/###/Vigna_angularis/Fasta/ --out /###/###/outputdir/ --parallel 8

3. Training on a genome and predicting on another:
    ```bash
    ./look4ltrs --fasta /###/###/Phaseolus_vulgaris/Fasta/ --train /###/###/Vigna_radiata/Fasta --out /###/###/outputdir/ --parallel 8

## License
Academic use: The software is provided as-is under the GNU GPLv3. Any restrictions to use for-profit or non-academics: License needed.
