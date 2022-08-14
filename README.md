# Genotype and Phenotype Network (GPN)

A novel method for multiple phenotype association studies by constructing a bipartite signed network, linking phenotypes and genotypes into a Genotype and Phenotype Network (GPN), which is a new insight to investigate correlation among phenotypes. The GPN can be constructed by both quantitative and qualitative traits, especially phenotypes have extremely unbalanced case-control ratios.

This package "GPN" includes the following parts:

- **GPN construction**: Construction of the Genotype and Phenotype Network with the individual-level phenotype and genotype data;
- **Community detection method**: Novel community detection method to partition K phenotypes into disjoint network modules based on the similarity matrix from GPN;
- **Modified computational efficient clustering linear combination method**: Association test between phenotypes in each module and a SNP of interest;
- **Comparison methods**: CLC, HCLC, MultiPhen, O'Brien, and Omnibus.


## Installation

You can install the released version of APGD from Github with:

``` r
devtools::install_github("xueweic/GPN")
```

## Reference
Xuewei Cao, Shuanglin Zhang, Qiuying Sha*. A novel method for multiple phenotype association studies based on genotype and phenotype network.

**Any questions**? xueweic_AT_mtu_DOT_edu


