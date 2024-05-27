# signor2Hipathia
In this project, we aim to integrate pathways from the Signor database to make them compatible with the Hipathia method. 
Here is an overview of the project and instructions for running the main R script.

## Introduction

The main R script in this project parses pathways from the Signor database [signor.uniroma2.it] and prepares them for analysis using the Hipathia method. 
It retrieves pathway data from the Signor API, processes it, and generates output files in a format suitable for further analysis.

## Usage

To use the main R script, follow these steps:

1. Ensure that you have R installed on your system.
   - Note: This parser was tested and developed for R versions 4.1.2 and 4.3.1. If it does not function under other versions that you use, please report it to us.
2. Clone this repository to your local machine.
   ```
   git clone https://github.com/kinzaR/signor2Hipathia.git
   ```

3. Navigate to the directory containing the main R script (`main.R` and `signor2Hipathia.sh`).
4. You have three options:
- Use the `signor2Hipathia.sh` file with your installed R version as the first argument from your shell:
  ```
  example: ./signor2Hipathia.sh 4.1.2 --help (for help)
  ```
  This option is recommended for users with multiple installed R versions.
- Use the `main.R` script:
  ```
  example: ./main.R --help
  ```
- Alternatively, open it using R Studio and adjust the arguments as needed.

5. Set any desired options in the script, such as the species variable or output folder. Most options have default values to optimize processing time and ease of use.
   available options are :
   ```
	-s SPE, --spe=SPE
		Species variable. Allowed choices: 'hsa', 'mmu', 'rno'. (default: hsa)

	-v, --verbose
		Enable verbose mode.

	-r, --readySifs
		Read from already created sif and att files.

	-p PATHWAYS_LIST, --pathways_list=PATHWAYS_LIST
		Vector of the IDs of the pathways to parse fromSignor. By default, all available pathways are loaded. Example: 'SIGNOR-AML,SIGNOR-LBC,SIGNOR-PDAP'.

	--score=SCORE
		The minimum significance score allowed, the range is from 0.1 to 1. Signor set 0.1 as the minimum score as 0 stands for no evidence of interaction.[ for more information : https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://signor.uniroma2.it/documentation/SIGNOR_3_score_Documentation_final.docx&ved=2ahUKEwim0JGkkI2GAxX1YPEDHRT8BKIQFnoECBoQAQ&usg=AOvVaw2y_b2VjYMFJgoA3BilRe95]
              By default 0.1

	-o OUTPUT_FOLDER, --output_folder=OUTPUT_FOLDER
		Output folder

	-h, --help
		Show this help message and exit
   ```
6. Example to run the script:
   ```
    ./signor2Hipathia.sh 4.3.1 -p SIGNOR-AML -v -o tmp/my_first_report OR ./main.R -p SIGNOR-AML -v -o tmp/my_first_report
   ```
   This command selects only the "SIGNOR-AML: Acute Myeloid Leukemia" pathway for parsing [https://signor.uniroma2.it/pathway_browser.php?pathway_list=SIGNOR-AML]. as a result from the commandLine, you will see:
    ![image](https://github.com/kinzaR/signor2Hipathia/assets/12510444/c41df08d-a5ba-4f0d-ab16-1d1a7871863a)

   You can also view the web report:
   ![image](https://github.com/kinzaR/signor2Hipathia/assets/12510444/b1be3af5-c7d9-4478-98ac-3b575ce5a87b)

   For more information about all molecules/nodes, a tooltip has been added containing a direct link to its record in the Signor database:
   ![image](https://github.com/kinzaR/signor2Hipathia/assets/12510444/3b3c5d75-b862-48d7-9eb8-ff9c148f1353)


## Dependencies

The main R script relies on the following R packages:

- `magrittr`
- `optparse`
- `hipathia`
- Additional functions are sourced from `utils.R` and `configs.R`.

## Output

The output of the main R script includes:

- RDS file: the meta-graph-info object for Hipathia algorithm for further mechanistic modeling analysis.
- SIF and ATT files for each parsed pathway.
- Phenotype and stimulus annotations files.
- JSON file containing pathway information (ids and names).
- A folder of pathway-viewer after parsing it, to visualize it localy.
- `log.txt` file: Contains all warnings and information about all chosen pathways, indicating if they were parsed or not, and reasons if they were not parsed.
- `not_parsed.tsv` file: Contains the IDs of pathways that were not parsed.
- 
## Testing and Validation

This parser is currently undergoing testing and validation. A generated report using all parsed 84 pathways from Signor and the Breast Cancer (BRCA) dataset can be found [here](http://hipathia.babelomics.org/signor_tests/pathway-viewer/). Please note that this report is still under review by the team.

![image](https://github.com/kinzaR/signor2Hipathia/assets/12510444/442cece7-6748-414e-a7d2-09bc5eba1746)


Reproducibility: To work with the Breast Cancer dataset from The Cancer Genome Atlas repository, you can download the expression matrix and the experimental design from these links:

- Expression matrix: [brca_genes_vals_bn.txt](https://github.com/kinzaR/signor2Hipathia/blob/main/files/data/brca_genes_vals_bn.txt)
- Experimental design: [brca_normal-basal_ed.txt](https://github.com/kinzaR/signor2Hipathia/blob/main/files/data/brca_normal-basal_ed.txt)

## Conclusion
The signor2Hipathia project aims to facilitate the integration of Signor pathways into the Hipathia method. By following the instructions provided in this README, users can parse Signor pathways and prepare them for further analysis using Hipathia.
For more information about the project or assistance with running the script, please refer to the documentation or contact us at kinza.rian@juntadeandalucia.es.

## Notes:

- "Protein families" and "fusion proteins" were treated the same, with the separator `,`.
  Complexes were separated by `/`.
- Small molecules, chemicals, and miRNA were kept with NA in gene lists.
- `,or,` were replaced by a simple `,`, and `,and` by `/`.
- Annotations from UniProt to Entrez: Some proteins (with UniProt id) has no Entrez ID, these were kept puting NA in the genesList.
