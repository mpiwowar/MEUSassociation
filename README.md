# MEUSassociation

This is an R library that implements Z coefficient, a measure of directional association for contingency tables.

## Prerequisites

To use this package one needs to have the R program installed. The package does not depend on any additional R packages.

## Installing

To install the package one can download the tar.gz file in tha latest version and install it from the local copy in R using install.package R command. For more information see R help.

## Running sample code

The followin presents sample code that can be used to see how the package works:
library(MEUSassociation)
> library(MEUSassociation)
> data("cancer_mutations")
> z_coefficient(cancer_mutations)
> data("cancer_mutations_gene_groups")
> z_coefficient(cancer_mutations, row_groups = cancer_mutations_gene_groups)
> head(z_coefficient_ranks(cancer_mutations))
> head(z_coefficient_ranks(cancer_mutations, col_groups = cancer_mutations_gene_groups)) 

## Authors

* **Monika Piwowar** - *Author*
* **Tomasz Kulaga** - *Maintainer, Author*

## License

This project is licensed under the GPL-3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The definition of the Z coefficient is acredited to Z. Meus.
