
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VarSubnetDetection

This supports the detection and analysis of altered subnetworks for
epilepsy mutations using the hierarchical hotnet algorithm \[Reyna,
Leiserson, and Raphael
(2018)\]\[<https://github.com/raphael-group/hierarchical-hotnet.git>\].
The two WES epilepsy datasets CoGIE and EpiPGX were analyzed in the
study (May et al. 2018) and the loss-of-function and Non-synonmyous
mutations were considered in this study.

# Study design

<img src="figures/rare_variant_wfl.png" title="workflow"
style="width:40.0%" />

# Implementation details

The datasets were received in-house. The analysis was carried out in two
parts.

-   In the first part, subnetworks and genes were detected for the
    above-mentioned datasets using the HPC cluster Iris at the
    University of Luxembourg (Varrette et al. 2014).

-   In the second part, over-representation of detected genes and
    networks were further analyzed via the clusterProfiler (Yu et
    al. 2012) functions in R.

For network construction, gene interaction information were considered
from the following public databases,

-   STRING (Szklarczyk et al. 2019)
-   humanbase (Krishnan et al. 2016)
-   Reactome (Fabregat et al. 2018)

# Contact

For questions and comments, please contact <sreyoc2999@gmail.com>

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-fabregat2018reactome" class="csl-entry">

Fabregat, Antonio, Florian Korninger, Guilherme Viteri, Konstantinos
Sidiropoulos, Pablo Marin-Garcia, Peipei Ping, Guanming Wu, Lincoln
Stein, Peter D’Eustachio, and Henning Hermjakob. 2018. “Reactome Graph
Database: Efficient Access to Complex Pathway Data.” *PLoS Computational
Biology* 14 (1): e1005968.

</div>

<div id="ref-krishnan2016genome" class="csl-entry">

Krishnan, Arjun, Ran Zhang, Victoria Yao, Chandra L Theesfeld, Aaron K
Wong, Alicja Tadych, Natalia Volfovsky, Alan Packer, Alex Lash, and Olga
G Troyanskaya. 2016. “Genome-Wide Prediction and Functional
Characterization of the Genetic Basis of Autism Spectrum Disorder.”
*Nature Neuroscience* 19 (11): 1454–62.

</div>

<div id="ref-may2018rare" class="csl-entry">

May, Patrick, Simon Girard, Merle Harrer, Dheeraj R Bobbili, Julian
Schubert, Stefan Wolking, Felicitas Becker, et al. 2018. “Rare Coding
Variants in Genes Encoding GABAA Receptors in Genetic Generalised
Epilepsies: An Exome-Based Case-Control Study.” *The Lancet Neurology*
17 (8): 699–708.

</div>

<div id="ref-reyna2018hierarchical" class="csl-entry">

Reyna, Matthew A, Mark DM Leiserson, and Benjamin J Raphael. 2018.
“Hierarchical HotNet: Identifying Hierarchies of Altered Subnetworks.”
*Bioinformatics* 34 (17): i972–80.

</div>

<div id="ref-szklarczyk2019string" class="csl-entry">

Szklarczyk, Damian, Annika L Gable, David Lyon, Alexander Junge, Stefan
Wyder, Jaime Huerta-Cepas, Milan Simonovic, et al. 2019. “STRING V11:
Protein–Protein Association Networks with Increased Coverage, Supporting
Functional Discovery in Genome-Wide Experimental Datasets.” *Nucleic
Acids Research* 47 (D1): D607–13.

</div>

<div id="ref-VBCG_HPCS14" class="csl-entry">

Varrette, S., P. Bouvry, H. Cartiaux, and F. Georgatos. 2014.
“Management of an Academic HPC Cluster: The UL Experience.” In *Proc. Of
the 2014 Intl. Conf. On High Performance Computing & Simulation (HPCS
2014)*, 959–67. Bologna, Italy: IEEE.

</div>

<div id="ref-clusterProfiler" class="csl-entry">

Yu, Guangchuang, Li-Gen Wang, Yanyan Han, and Qing-Yu He. 2012.
“clusterProfiler: An r Package for Comparing Biological Themes Among
Gene Clusters.” *OMICS: A Journal of Integrative Biology* 16 (5):
284–87. <https://doi.org/10.1089/omi.2011.0118>.

</div>

</div>
