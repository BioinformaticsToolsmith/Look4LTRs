# Look4LTRs
Long Terminal Repeat Retrotransposon detection tool capable of finding recently nested LTR RTs de novo.

Copyright (C) 2023 Anthony B. Garza and Hani Z. Girgis, PhD 

Academic use: Affero General Public License version 1.

Any restrictions to use for profit or non-academics: Alternative commercial license is required.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.

The src directory holds the source code.
The CMakeLists is the cmake file.
The script directory holds useful scripts for studying the results of Look4LTRs. These scripts were made with Python 3.10.

## Requirements

GNU g++ 11.1.0
cmake 3.10.3

If you do not have the required g++ and cmake, we suggest creating a conda environment.
Please follow these steps to install the required g++ and cmake, assuming an enviornment called "myenv".
You may change the environment name.

```
conda create -n myenv
conda activate myenv
conda install -c conda-forge cmake=3.10.3
conda install -c conda-forge gcc=11.1.0
```

## How to Compile

```
mkdir bin
cd bin
(If your default compiler meets the version requirement) 
cmake ..
(Or if you would like to specify a different compiler that meets the requirement)
cmake .. -DCMAKE_CXX_COMPILER=your_compiler_name_for_example_g++-7
(Or if this fails, try using this to set a different compiler. Replace the paths with your own.)
cmake .. -DCMAKE_CXX_COMPILER=$HOME/C++/GCC/bin/g++ -DCMAKE_C_COMPILER=$HOME/C++/GCC/bin/gcc -DCMAKE_PREFIX_PATH=$HOME/C++/GCC
make look4ltrs
```

## INPUT
Look4LTRs accepts FASTA format files, as well as multi-FASTA format. It is suggested that, at minimum, an entire genome is given to Look4LTRs, to enhance its self-supervised capabilities. Multiple genomes may be passed into Look4LTRs, but take caution with the memory requirements, and lower the number of threads (if any) if too much memory is being utilized. Additionally, Look4LTRs accepts training genomes that will not be predicted upon but are used to enhance the prediction of other genomes.

The FASTA format files given to Look4LTRs MUST have a .fa extension. No other extension is allowed!

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
  15. *GraphGroup* is the ID of a graph. Look4LTRs utilizes a graph when matching LTRs. Each graph represents nearby elements of the same family. This column is an ID to the graph the LTR RT/solo LTR belongs to.
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
| -t/--train | Fasta file directory for training only. Can be given a variable number of arguments.
| -o/--output | Output directory. | Yes |
| -c/--config | Config file that contains a machine learning model's parameters. Used to replace the model in the detector module and downstream parameters. | No |
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
    
When passing more than one fasta directory through --fasta, --train, or combined, give the full path to these directories. Otherwise, Look4LTRs can not properly build symbolic links.
    
## Scripts

| Script | Description | Usage | Example |
|-----------------|-----------------|-----------------|-----------------|
| findRecentNest.py | will find the IDs of recently nested LTR RTs and print out to terminal a tree of these nests. | Pass in a path to a directory of RTR files or a path to a single RTR file. If given a directory, each file, before the extension, must end with \_chr#, where # is the chromosome identifier. | python3.10 findRecentNest.py /###/###/outputdir/Rtr/Glycine_max_chr1.rtr |
| findSameGraphNest.p | will find the IDs of LTR RTs nested into a same-graph LTR RT (same family possibly). Prints to a tree like findRecentNest.py. | Pass in a path to a directory of RTR files or a path to a single RTR file. If given a directory, each file, before the extension, must end with \_chr#, where # is the chromosome identifier. | python3.10 findSameGraphNest.py /###/###/outputdir/Rtr/Glycine_max_chr1.rtr |
| findRT.py | will return the line belonging to an LTR RT from an RTR file given its ID. Use this in conjunction with the above scripts instead of searching by hand. | Pass in a path to a single RTR file. | python3.10 findRT.py /###/###/outputdir/Rtr/Glycine_max_chr1.rtr |

findSameGraphNest.py 
findRT.py will return the line belonging to an LTR RT from an RTR file given its ID. Use this in conjunction with the above scripts instead of searching by hand.

## Training
Look4LTRs is trained on *Arabidopsis thaliana*, *Oryza sativa japonica*, *Glycine max*, and *Sorghum bicolor* using elements delineated by RepeatMasker and Repbase (2018). Non-model organisms may not be well-represented as a result. To accomodate this, the pipeline for training Look4LTRs has been provided.

We advise caution with retraining Look4LTRs as it may result in unexpected results. In Look4LTRs, an SGD classifier (from scikit-learn) was trained on the aforementioned genomes. The SGD classifier itself was chosen after consideration of other models such as linear regressors and random forests. As such, the SGD classifier may not result in the optimal results on non-model genomes.

A few steps are required to set up the training of Look4LTRs.

1. **Compile required executables**. Two executables need to be compiled for the training pipeline to work. Navigate to the bin directory and run the following commands.
    ```
    cd bin
    make generateTrainingData
    make generateGraphData
    ```

2. **Set up the input genomes**. The training pipeline requires the FASTA files of the genomes to train on.
    - a. Each genome the user wishes to train on must be split into their own directories.
    - b. No multi-fasta format files allowed. Each FASTA file may only have one sequence. Each chromosome should be its own file. This simplifies the mapping of training data to chromosomes.
    - c. The FASTA files must have the '.fa' extension.

3. **Set up the LTR-retrotransposon locations**. The user must provide the locations of LTR-retrotransposons within their provided genomes in BED format.
    - a. Just like the genomes' FASTA files, the BED files for the genomes must be split into separate directories.
    - b. Each BED file must correspond to exactly one chromosome.
    - c. The name of each BED file must be exactly the same as the corresponding FASTA file's. The only difference should be the extension.
    - d. The BED files must have the following columns:
      - *chrom* is the chromosome identifier
      - *start* is the start position of the LTR-retrotransposon
      - *end* is the end position of the LTR-retrotransposon
      - *left_start* is the start position of the upstream LTR (typically the 5' LTR)
      - *left_end* is the end position of the upstream LTR
      - *right_start* is the start position of the downstream LTR (typically the 3' LTR)
      - *right_end* is the end position of the downstream LTR. 
    - e. The BED files must include the header line of column names


4. **Download Python**. Python 3.8 and higher is recommended. A conda environment can be used.

5. **Download the required packages**. In the Training folder, there is a requirements.txt file. Please install the packages within this file for Python.
    ```
    cd ../Training
    pip install -r requirements.txt
    ```

The pipeline can be run by calling **trainModel.py** located in the Training folder. This script takes the following parameters.


| Parameter | Description | Required? |
|-----------------|-----------------|-----------------|
| -fd/--fastadirs | Takes multiple arguments. The paths to each genome's FASTA file directory. Separate each path with a space | Yes |
| -bd/--beddirs | Takes multiple arguments. The paths to each genome's BED file directory. Separate each path with a space. | Yes |
| -o/--output | Output directory. If it doesn't exist, it will be created (assuming the base path exists). WARNING: If it does exist, everything in the folder will be deleted beforehand. | Yes |

An example usage is the following:

    python3 trainModel.py -fd /###/Genome1/Fasta/ /###/Genome2/Fasta/ -bd /###/Genome1/Bed/ /###/Genome2/Bed -o /###/Output/
    

The result of this pipeline is located in the provided output folder in a file called **config.txt**. This file can then be passed to the **look4ltrs** executable with the -c/--config parameter detailed above.
