# Structure and metabolic potential of the prokaryotic communities from the hydrothermal system of Paleochori Bay, Milos, Greece.
Sven Le Moine Bauer, Guang-Sin Lu, Steven Goulaouic, Valentine Puzenat, Anders Schouw, Thibaut Barreyre, Vera Pawlowski-Glahn, Juan Jose Egozcue, Jean-Emmanuel Martelat, Javier Escartin, Jan P. Amend, Paraskevi Nomikou, Othonas Vlasopoulos, Paraskevi Polymenakou, Steffen Leth Jørgensen


| ![](Picture_bubles.jpg) | 
|:--:| 
| *Intense hydrothermal degassing in Paleochori Bay. Picture: Anders Schouw* |


This repository contains the integrality of the scripts and files needed to reproduce the data analysis presented in the aforementioned article. Note that the authors are in no case Unix/R professionals, and the code can certainly be written in a more idiomatic way. Do not hesitate to reach out for further help. 

The following links will bring you to:
- [The processing of the sequences and picking of OTUs](Pipeline%20explanations.md)
- [The decontamination protocol for all OTUs](Decontamination_pipeline.md)
- The statistical analysis: 
  - [The CoDA analysis (Figures 2A, 2B, 5, Supplementary material 5)](CoDA_analysis.md)
  - [The qPCR analysis (Figure 2D, 3)](qPCR_analysis.md)
  - [The barplots (Figure 4, Supplementary material 7)](Barplots.md)
  - [The Shannon diversity analysis (Figure 2C)](Shannon_analysis.md)

The following files are also given:
- [Metadata.csv](Metadata.csv): The context information to all samples taken.
- [Otutab.sorted.tsv](Otutab.sorted.tsv): The OTU table produced by the sequence processing pipeline, prior to decontamination.
- [assignments.txt](assignments.txt): The taxonomy information assigned by CREST4, prior to decontamination.
- [otutab_decontam.csv](otutab_decontam.csv): The OTU table post decontamination, used in the statistical analysis.
- [tax_decontam.csv](tax_decontam.csv): The taxonomy information post decontamination, used in the statistical analysis.
- [Otus_decontam.fasta](Otus_decontam.fasta): Centroids of the OTUs after decontamination. Not used in any script presented here.

==THIS DEPOSITORY IS NOT COMPLETE==
