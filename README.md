# Genotype and Phenotype Network (GPN)

A novel method for multiple phenotype association studies by constructing a bipartite signed network, linking phenotypes and genotypes into a Genotype and Phenotype Network (GPN), which is a new insight to investigate correlation among phenotypes. The GPN can be constructed by both quantitative and qualitative traits, especially phenotypes have extremely unbalanced case-control ratios.

This package "GPN" includes the following parts:

- **GPN construction**: Construction of the Genotype and Phenotype Network with the individual-level phenotype and genotype data;
- **Community detection method**: Novel community detection method to partition K phenotypes into disjoint network modules based on the similarity matrix from GPN;
- **Multiple Phenotype Association Tests with and without considering the network modules**: test the association between K phenotypes with a SNP. The phenotypes can be either qualitative or binary, espectially the binary phenotypes with the extremely case-control ratio (the test statistics has been adjusted by the saddlepoint approximation)
  - ceCLC: Modified computational efficient clustering linear combination method
  - CLC: Modified clustering linear combination method
  - HCLC: Modified clustering linear combination method based on hierarchical clustering
  - MultiPhen: joint model of multiple phenotypes can increase discovery in GWAS
  - O'Brien: Modified O'Brien (SSU) method
  - Omnibus: Modified omnibus (Chi-square) method
  - TATES: Modified TATES method


## Installation

You can install the released version of APGD from Github with:

``` r
devtools::install_github("xueweic/GPN")
library(GPN)
```

## Reference
Xuewei Cao, Shuanglin Zhang, Qiuying Sha*. A novel method for multiple phenotype association studies based on genotype and phenotype network.

**Any questions**? xueweic_AT_mtu_DOT_edu


## Graphical Overview

<p align="center">
  <img src="Figure/Figure1 copy.png" width="1000">
</p>

- **a.** Construction of a signed bipartite Genotype and Phenotype Network (GPN). Each phenotype (yellow square) and each SNP form a directed edge which represents the strength of the association, where the red dashed line indicates that the minor allele of the SNP is a protective allele to the phenotype, and the blue dashed line indicates that the minor allele of the SNP is a risk allele to the phenotype.
- **b.** Construction of a signed network, named the Phenotype and Phenotype Network (PPN), which is the one-mode projection of GPN on phenotypes. 
- **c.** The powerful community detection method is used to partition phenotypes into disjoint network modules with different colors. 
- **d.** Multiple phenotype association tests are applied to test the association between phenotypes in each of the network modules and a SNP, then the Bonferroni correction is used to obtain the overall p-value.  
- **e.** GWAS signals are identified by a multiple phenotype association test with or without considering network modules.
- **f.** Functional enrichment analysis based on the detected GWAS signals and the publicly available functional database. 
- **g.** Colocalization of GWAS signals and eQTL analysis.





