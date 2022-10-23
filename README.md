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


