# Evaluating the impact of ROS on *P. putida* metabolism
This codes allow to assess the impact of ROS-mediated inhibition of enzymes in the metabolism of *P. putida*. 

Codes used for this aim were programmed in Jupyter Notebooks for a better manipulation and step-by-step evaluation. 

## What does each script do? 
- `metab_enrich_analysis.ipynb`: Assesses whether each metabolic subsystem belongsto primary or secondary metabolism according to enrichment analyses. 
- `Pputida_ROS_modelling.ipynb`: Performs single knock-out and subsystems analyses in experimentally evidenced susceptible reactions (those in [Pputida_evidenced.txt](../../data/Pputida_evidenced.txt).
- `Pputida_ROS_modelling_extended.ipynb`: Performs single knock-out and subsystems analyses in potentially susceptible reactions (those with at least one ROS-susceptible structural feature, can be found in [ROS_summary_full.tsv](../../data/ROS_summary_full.tsv).
