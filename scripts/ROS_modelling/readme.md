# Evaluating the impact of ROS on *P. putida* metabolism
This codes allow to assess the impact of ROS-mediated inhibition of enzymes in the metabolism of *P. putida*. 

Codes were programmed as Jupyter Notebooks for a better manipulation and step-by-step evaluation. 

## What does each script do? 
- `metab_enrich_analysis.ipynb`: Assesses whether each metabolic subsystem belongsto primary or secondary metabolism according to enrichment analyses. 
- `Pputida_ROS_modelling.ipynb`: Performs single knock-out and subsystems analyses in experimentally evidenced susceptible reactions (those in [Pputida_evidenced.txt](../../data/Pputida_evidenced.txt)).
- `Pputida_ROS_modelling_extended.ipynb`: Performs single knock-out and subsystems analyses in potentially susceptible reactions (those with at least one ROS-susceptible structural feature, can be filtered in [ROS_summary_full.tsv](../../data/ROS_summary_full.tsv)).
- `ROS_metab.ipynb`: Performs Flux Balance Analysis (FBA), Flux Variability Analysis (FVA) and flux sampling, comparing an updated version of model `iJN1480` and a ROS-constrained version of it (constrained using experimentally evidenced susceptible reactions). Requires modules from `sampling_utils.py`. 
- `ROS_metab_extended.ipynb`: Jupyter Notebook used for FBA, FVA and flux sampling analyses using the ROS-constrained model restricted with structural feature’s information.

## Dependencies
All scripts depend on the python module COnstrained-Based metabolic modeling in Python: [`COBRApy`](https://github.com/opencobra/cobrapy).
