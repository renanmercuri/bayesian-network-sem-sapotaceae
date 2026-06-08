# Bayesian Network and SEM Analysis of Phenotypic Traits in Sapotaceae

This repository contains the data, R scripts, and output files used to reproduce the analyses reported in the manuscript:

**Phenotypic dependency networks of fruit traits in mamey sapote and star apple germplasm**

The study investigates phenotypic dependency networks among fruit, peel, and seed traits in two Sapotaceae species:

- **Mamey sapote** (*Pouteria sapota* (Jacq.) H.E. Moore & Stearn)
- **Star apple** (*Chrysophyllum cainito* L.)

The analyses combine Bayesian network structure learning, leave one out jackknife stability assessment, inspection of Markov equivalence classes, biological interpretation, and structural equation modeling.

## Repository structure

```text
bayesian-network-sem-sapotaceae/
├── data/
│   ├── sapote.txt
│   └── caimito.txt
├── output/
│   ├── arcs_mamey_iamb_alpha005.csv
│   ├── arcs_mamey_iamb_alpha020.csv
│   ├── arcs_mamey_tabu_bge.csv
│   ├── arcs_star_iamb_alpha005.csv
│   ├── arcs_star_iamb_alpha020.csv
│   ├── arcs_star_tabu_bge.csv
│   ├── descriptive_mamey_sapote.csv
│   ├── descriptive_star_apple.csv
│   ├── final_sem_model_fit_summary.csv
│   ├── final_sem_path_diagrams.pdf
│   ├── jk_*_connection_percent.csv
│   ├── jk_*_direction_percent.csv
│   ├── jk_*_raw_arc_percent.csv
│   ├── path_coefficients_wald_final_models.csv
│   ├── path_coefficients_wald_mamey_final.csv
│   ├── path_coefficients_wald_star_final.csv
│   └── session_info.txt
├── scripts/
│   └── supplementary_R_script_sapotaceae_network_sem.R
├── LICENSE
└── README.md
```

## Data

The data were obtained from the systematic characterization of CATIE's Sapotaceae germplasm collection reported by Gazel Filho (1995). The analyses use phenotypic measurements from:

- 63 mamey sapote trees
- 49 star apple trees

Each record corresponds to one tree and summarizes measurements obtained from fruit samples. The analyses were conducted separately for each species.

The eight quantitative traits analyzed were:

| Code | Trait | Unit |
|---|---|---|
| FRW | Fruit weight | g |
| FRL | Fruit length | mm |
| FRD | Fruit diameter | mm |
| PET | Peel thickness | mm |
| PEW | Peel weight | g |
| SEL | Seed length | mm |
| SED | Seed diameter | mm |
| SEW | Seed weight | g |

## Statistical workflow

The main analysis is implemented in:

```text
scripts/supplementary_R_script_sapotaceae_network_sem.R
```

The script reproduces the complete workflow used in the manuscript:

1. Reads and prepares the phenotypic datasets.
2. Computes descriptive statistics for each species.
3. Assesses multivariate normality.
4. Applies natural logarithmic transformations to selected star apple traits.
5. Learns Bayesian network structures using:
   - IAMB with alpha = 0.05
   - IAMB with alpha = 0.20
   - Tabu search using the Bayesian Gaussian equivalent score
6. Performs leave one out jackknife stability assessment.
7. Exports connection occurrence, direction recovery, and raw arc matrices.
8. Inspects completed partially directed acyclic graphs.
9. Fits the final retained structural equation models.
10. Exports path coefficients, standard errors, Wald tests, model fit summaries, and session information.

## Bayesian network analysis

The Bayesian network analysis was performed using the `bnlearn` package. The IAMB algorithm was used as a constraint based learning procedure, and Tabu search was used as a score based learning procedure with the Bayesian Gaussian equivalent score.

The IAMB results were useful for identifying recurrent conditional dependence relationships among traits. Tabu search provided fully directed graphs that were interpreted together with jackknife stability, Markov equivalence class information, and biological plausibility.

## Jackknife stability assessment

Leave one out jackknife resampling was used to evaluate the stability of the learned networks.

Three types of jackknife output are provided:

- `connection_percent`: occurrence of a connection between two traits, irrespective of direction.
- `direction_percent`: recovery of a specific directed path.
- `raw_arc_percent`: raw adjacency matrix output from `bnlearn`.

Stable components were defined as those recovered in at least 80% of the leave one out resamples.

## Final structural equation models

The final structural equation models were specified from the joint evaluation of learned network structures, jackknife stability evidence, Markov equivalence class information, and biological plausibility.

### Mamey sapote final retained SEM

The final model includes the following directed paths:

```text
FRW -> PET
FRW -> PEW
FRW -> FRD
PET -> PEW
FRL -> FRD
PEW -> FRL
FRL -> SEL
SEL -> SEW
SEW -> SED
```

### Star apple final retained SEM

The final model includes the following directed paths:

```text
FRW -> PET
FRW -> PEW
FRW -> FRD
FRW -> FRL
PET -> PEW
FRL -> PEW
PEW -> FRD
FRD -> SEW
SEW -> SED
SED -> SEL
FRW -> SED
```

## Main output files

The main output files used in the manuscript are:

```text
output/final_sem_model_fit_summary.csv
output/path_coefficients_wald_final_models.csv
output/path_coefficients_wald_mamey_final.csv
output/path_coefficients_wald_star_final.csv
output/final_sem_path_diagrams.pdf
output/session_info.txt
```

The jackknife matrices are also included in the `output/` folder and allow the reproduction of the tables reporting stable connections and stable directions.

## Software requirements

The analyses were performed in R. The main packages used were:

```text
bnlearn
MVN
sem
psych
semPlot
```

The exact R version, package versions, operating system, and session details are reported in:

```text
output/session_info.txt
```

## How to reproduce the analysis

After cloning or downloading this repository, open R or RStudio and run:

```r
source("scripts/supplementary_R_script_sapotaceae_network_sem.R")
```

The script assumes that the data files are available in the `data/` folder and writes the output files to the `output/` folder.

## Citation

If you use this repository, please cite the corresponding manuscript and the archived Zenodo version associated with this repository.

Manuscript title:

```text
Phenotypic dependency networks of fruit traits in mamey sapote and star apple germplasm
```

Zenodo DOI:

```text
To be updated after the final release is archived.
```

## License

This repository is released under the MIT License. See the `LICENSE` file for details.

## Contact

For questions about the data, scripts, or reproducibility of the analyses, please contact:

**Renan Mercuri Pinto**  
GitHub: [@renanmercuri](https://github.com/renanmercuri)
