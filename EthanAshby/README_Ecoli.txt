README 
E. coli timing analysis and the sigmoidal model
Ethan Ashby
Email: ecadeashby@gmail.com, Phone #: (925) 586-4191

Hello! This Github folder contains all the data and code needed to produce a
writeup (E.coli_final_analysis.pdf) which details the analysis of the timing of E. coli 
gene expression in response to cell starvation. Gene expression timing analysis is
conducted through sigmoidal models fit through the ImpulseDE2 package. The properties of
these sigmoidal models are first assessed through simulations, and then are applied to the
data at hand.

This README will give you a tour of this folder and will walk you through how to reproduce
these analyses.

A brief tour of the 'EthanAshby' folder:
The main file in this folder is the 'E.coli_final_analysis.Rmd' file, which includes the
text and code required to generate the technical report: 'E.coli_final_analysis.pdf'. In addition,
there is a 'Next Steps' file that contains possible directions for future research. 
Included in the folder are several types of data files (detailed below). There are a few
other folders included in this master folder. The 'Extra' folder contains general miscellany (data,
code, etc.) generated during the summer of work, none of which are critical to the project
but may include some useful bits and pieces for future work. 'Images' contains screenshots 
that are linked into the final technical report. 'GO Clusters' and 'Go Stricter Clusters' 
contain gene ontologies for DEGs identified by the intersections of 2+ and 3 DE tools 
respectively; these DEGs were run through a fuzzy clustering pipeline to identify genes 
with similar expression profiles. The main goal of this report was not a clustering 
analysis, so these data should be considered supplemental. More detail on clustering 
can be obtained from Annie Cohen's folder.

Data files:
• meanExptoStationaryRpoSLevels.csv: this is the first data file in the Rmarkdown. This file
contains mean RpoS and mean culture density for 3 trials of 3 WT and ∆RpoS strains each over 150
Min time course. These data are used to make the RpoS kinetics figures on pg 4.
• LB_Time_Course_GCA_000005845.2_ASM584v2_genomic_counts.tsv: this is the second data file
in the Rmarkdown. This file contains read count data for genomic features (unparsed ident-
ifier info in the first column) over different experimental conditions and times. The file
stationary_phase_data_sample_list_revised.xlsx provides a key for the acronyms and the 
experimental conditions and times. For example, JH01_A01 denotes the first replicate of 
the WT strain at time 0. In the Data Pre-Processing section, the genomic features were 
parsed and the read counts were filtered for CDS.
	***A major consideration here: we only focus on CDS sequences because we were particularly
	interested in the timing of expressions of genes. However, other genomic features such
	as IGRs, ncRNAs, etc that are highly differentially expressed could also be of interest.
	We suggest this as a potential route of further research.
•Sensitivity.xlsx: These data are used in the results section under the Exploratory and Wilcox
R chunk. This file (provided by Prof. Daniel Stoebel) contains all the E. coli gene names, 
bnums, direction of regulation, and the sensitivity to RpoS. This allows us to match DEGs 
identified by our pipeline to sensitive, insensitive, and linear genes.
•genesofinterest.txt: These data are used under the 'Overlay Trajectories and Relative mRNA plot'
R chunk. This is simply a list of genes, their bnums, and their sensitivities. These genes
were hand-picked by Prof. Stoebel as genes he was particularly interested in.
•uniprot-yourlist_M201907266746803381A1F0E0DB47453E0216320D5EB8C6I.tab: These data are used
in the Exploratory Analysis R chunk and is located within the 'Extra' folder. This file 
contains predicted functions, GO terms, metal associations, etc. This file was generated 
by submitting the DEG list to uniprot online GUI and selecting interesting features.
	***Mining UNIPROT Db and finding ways to integrate Rstudio with this Db more seamlessly
	could be a direction of further research.
	
General Steps Taken in Rmarkdown Document:
The following section goes through the Rmarkdown document chunk-by-chunk, and explains
the process and thought behind each section.

#Data processing, exploratory plots

1) Rpos kinetics: RpoS data from 'meanExptoStationaryRpoSLevels.csv' is visualized thru
ggplot2, demonstrating that RpoS expression is induced upon cell starvation time course.
2) Data-pre-processing (5 unnamed chunks): important packages and functions are loaded.
Read count data is read into R and the gene identifiers are parsed and filtered for CDS
sequences.
3) DE tools: Runs DESeq2, maSigPro, and ImpulseDE2 on the read count data. Read count data
is pre-normalized (using the DESeq2 RLE method) prior to passage through the maSigPro pipeline.
Venn diagrams are plotted to visualize the overlap between the tools.
4) DEGs: the intersection of 2+ (2 or more) tools is used to generate a DEG list.
5) Visualizing t params: generates density plots of t (onset time) parameters for differentially expressed genes. Gives a good picture of when all differentially expressed genes are 'turning on'.
6) Plotting Other and Transient DEGs: not called in the pdf. Just overlays plots all the transient DEG
and explores the genes that are ID'ed as DE by ImpulseDE2 but aren't labeled as 'transient' or 'monotonic'.

###Simulation work

7) Normal Dispersion Sim: general idea is to understand the stability/reproducibility of the sigmoid parameters. Intuition is to select a gene and it's fitted sigmoid profile, extract the estimated dispersion, then generate 100 simulated sigmoid profiles from the initial one by adding in the correct amount of Negative Binomial noise (w/ the dispersion estimate as the argument). This chunk produces density plots of the simulated parameters with the initial parameters colored, as well as a plot of the initial trajectory (in black) and all the simulated trajectories (colored by spearman's correlation).
8) More Noise: same as #7, except the dispersion-based noise is increased by a factor of 100.
9) Less Noise: Not included in the pdf. Same as #7, except the dispersion-based noise is decreased by a factor of 100.
10) More TPs: same as #7, except the initial sigmoid function is used to calculate mean read counts are additional time points, simulating an experiment with more TPs. 
11) Outlier Sim: same as #7, except with outlier data points. The simulated trajectories of 10 genes remained unchanged, while six groups of 15 genes had their read counts for one replicate at each of the six time points increased by a factor of 10. To illustrate, 15 simulated genes had one replicate's read counts increased by a factor of 10 at the first time point. 15 simulated genes had one replicate's read counts increased by a factor of 10 at the second time point.
12) Several Genes: This is the code chunk used to generate the screenshots of genes and their simulated trajectories for many different genes. Dispersion estimates are generated from ImpulseDE2 (unchanged from ImpulseDE2 estimates).

###Applying model to data at hand

13) Timing Plots and Wilcox: Plots onset times (t params) faceted by sensitivity to RpoS. Also conducts Wilcox rank sum tests on onset times for the sensitivity groups (both with/without sensitive genes that turn on abnormally late). Generates list of 8 sensitive genes that turn on late.
14) Overlay Trajectories and Relative mRNA plot: not included in pdf. Overlays the sigmoid trajectories of sensitive and insensitive DEGs, faceted by sensitivity type and alphaed by Prof Stoebel's 'interest' in them. Also plotted are IQRs for the sigmoid functions at each time point, colored by sensitivity group.
15) Overlay Sigmoids: Overlays the sigmoid trajectories of sensitive and insensitive DEGs, colored by sensitivity type. Also plots the median sensitive and insensitive sigmoid profile for sensitive and insensitive DEGs.
16) Exploratory Analysis: not printed in pdf. Assorted code from things I tried that didn't end up being all that fruitful. I tried some clustering based off model params, visualizing scatterplots of model params colored by different functional info (from the uniprot dataset), and functional (GO) enrichment
17) Overlay Sigmoids RpoS x-axis: not printed in pdf. We were interested in rerunning the timing analysis by fitting the models to RpoS concentrations instead of time. Unfortunately, simply supplying RpoS concentrations instead of time to ImpulseDE2 didn't work (the program never finished running). So I used a spline to predict RpoS concentrations from input times, and then just had the RpoS concentrations on all the plot axes. Didn't really change much, but coud be a direction of future research.


