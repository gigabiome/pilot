# Glossary

| Term | Meaning |
|---|---|
| **AMR** | Antimicrobial resistance. |
| **AMRFinderPlus** | NCBI tool that finds known resistance genes and point mutations in a genome. Provides the T7 reference calls. |
| **AST** | Antimicrobial susceptibility testing: the lab measurement of resistance. |
| **ATU** | Area of technical uncertainty: a EUCAST MIC range where R/S calls are unreliable. |
| **Bakta** | Genome annotation tool (gene coordinates and names). |
| **Breakpoint** | The MIC threshold that separates S from R for a species × drug. Set by EUCAST. |
| **Broth microdilution** | Reference AST method: growth is tested across two-fold antibiotic dilutions. |
| **Canonical k-mer** | The lexicographically smaller of a k-mer and its reverse complement, so the feature does not depend on strand. |
| **CARD** | Comprehensive Antibiotic Resistance Database. |
| **Censored MIC** | A value only known as a bound (`<=0.0625`, `>64`) because it lies outside the tested range. |
| **Clone cluster** | A group of near-identical genomes (here: k-mer Jaccard distance below a threshold). It is the unit kept together by the splits. |
| **Confounder** | A feature correlated with resistance without causing it, e.g. a plasmid backbone gene or a lineage marker. |
| **Coverage** (sequencing) | Average number of reads covering each base. |
| **Coverage** (harness) | Fraction of required test pairs a submission actually predicted. |
| **Decoy** | In the mock data, a planted non-causal plasmid backbone gene. |
| **ECE** | Expected calibration error: how far predicted probabilities are from observed frequencies (0 = perfect). |
| **EA (essential agreement)** | Fraction of predictions within ±1 two-fold dilution of the measured MIC. The primary metric for T1–T4. |
| **ESKAPE-E** | *Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, Acinetobacter baumannii, Pseudomonas aeruginosa, Enterobacter* spp., *Escherichia coli*. |
| **EUCAST** | European Committee on Antimicrobial Susceptibility Testing: publishes breakpoints and methods. |
| **F1** | Harmonic mean of precision and sensitivity (R = positive class). |
| **GNN** | Graph neural network. Here it operates on gene-adjacency graphs. |
| **Harness** | The scoring software: reads submissions, compares them to hidden labels, outputs scores. |
| **Intrinsic resistance** | Resistance shared by a whole species, not caused by an acquired gene or mutation. |
| **Jaccard distance** | 1 − (shared k-mers / all k-mers) between two genomes. |
| **k-mer** | A substring of length k of a DNA sequence. |
| **Leakage** | Test information reaching training, which inflates scores (e.g. clones on both sides of a split). |
| **Major error (ME)** | Predicted R, actually S. |
| **Manifest** | JSON file listing the isolates of each partition, plus provenance and audit. |
| **Mash** | Fast genome-distance tool based on k-mer sketches. |
| **MIC** | Minimum inhibitory concentration: the lowest antibiotic concentration that stops growth (mg/L). |
| **MIxS** | Minimum Information about any (x) Sequence, the metadata standard. |
| **MLST / ST** | Multi-locus sequence typing / sequence type: a standard strain classification. |
| **Pan-genome** | All gene families found across a set of genomes. The presence/absence matrix is the pan-genome model's input. |
| **Panaroo** | Pan-genome tool that clusters genes into families and outputs the presence/absence matrix. |
| **Plasmid** | Mobile DNA molecule; often carries resistance genes between strains and species. |
| **Point mutation / SNP** | A single-base change, e.g. in gyrA, causing fluoroquinolone resistance. |
| **PopPUNK** | Tool for defining bacterial strain clusters from genomes. |
| **R / S / I** | Resistant / susceptible / susceptible, increased exposure (EUCAST categories). |
| **Ridge regression** | Linear regression with an L2 penalty. The k-mer baseline. |
| **StratifiedGroupKFold** | A split that keeps groups (clone clusters) intact while balancing strata (species × site). |
| **Submission** | Folder of prediction files a model hands to the harness. |
| **Top-k recall** | Fraction of an isolate's known determinants that appear in the model's k highest-ranked features. |
| **Unicycler** | Genome assembler used in WP2. |
| **Unitig** | A non-branching path in a de Bruijn graph: a compact alternative to raw k-mers. |
| **VME (very major error)** | Predicted S, actually R. The clinically dangerous error. |
