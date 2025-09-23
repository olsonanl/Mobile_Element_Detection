# Mobile Element Detection Service

## Overview

The Mobile Element Detection Service takes as input a set of reads and assembles them using the BV-BRC assembly service, or a set of assembled contigs.  It uses the geNomad tool (https://github.com/apcamargo/geNomad) to identify viral and plasmid contigs, as well as proviral regions within the assembled contigs.  The underlying geNomad tool employs a hybrid classification approach that combines an alignment-free neural network with a gene-based classifier that utilizes a database of marker protein profiles specific to chromosomes, plasmids, and viruses. Genomad also provides taxonomic assignment of viral contigs using ICTV taxonomy. Once identified, the viral and plasmid contigs are submitted to the BV-BRC annotation service.  If the contigs are viral, and they match one of the supported LowVan taxa, then the LowVan recpie is used for the annotation.  Likewise, if it is a phage or pro-virus (prophage) Phannotate is used. If the contigs are identifed as being plasmid, then recipe is default RAST.  

###  Contigs that will be annotated by the Mobile Element Detection Service 

| ICTV taxon          | NCBI Taxonomy ID | Annotation recipe |
|---------------------|------------------|-------------------|
| Bunyaviricetes      | 1980410          | LowVan            |
| Paramyxoviridae     | 11158            | LowVan            |
| Filoviridae         | 11266            | LowVan            |
| Pneumoviridae       | 11244            | LowVan            |
| Orthomyxoviridae    | 11308            | LowVan            |
|                     |                  |                   |
| Caudoviricetes      | 2731619          | Phannotate        |
| Loebvirae           | 2732090          | Phannotate        |
| Sangervirae         | 2732091          | Phannotate        |
| Vidaverviricetes    | 2732460          | Phannotate        |
| Leviviricetes       | 2842243          | Phannotate        |
| Picobirnaviridae    | 585893           | Phannotate        |
| Partitiviridae      | 11012            | Phannotate        |
| Matsushitaviridae   | 3152134          | Phannotate        |
| Ainoaviricetes      | 3044425          | Phannotate        |
| Vinavirales         | 2732560          | Phannotate        |
| Prepoliviricotina   | 3412695          | Phannotate        |
| Obscuriviridae      | 3424590          | Phannotate        |
| Plasmaviridae       | 10472            | Phannotate        |
|                     |                  |                   |
| Plasmids*           | 2  (Bacteria)    | RAST (default)    |

*Identified plasmids are assumed to be bacterial, but if this is not the case,
no bacterial-specific assumptions are made by the default RAST workflow.


The app currently does a simple run of geNomad, to identify plasmids and viruses.  In future versions we plan to add our own plasmid classification tools. 


## Inputs and outputs

Inputs:  Reads or Assemblies

### Mobile Element Detection Service Output Files:
Directory     | File Name | Description |
|----------------|---------------|-------------|
|Project/assembly| ./assembly    | All output files of the optional assembly run|
|Project/        | contigs.fasta | Your uploaded contigs|
|Project/        | final.contigs_plasmid.fna | FASTA file containing plasmid contigs |
|Project/        | final.contigs_plasmid_summary.tsv | Tab-separated summary file for plasmid contigs |
|Project/        | final.contigs_summary.json | JSON file containing overall summary of all contigs |
|Project/        | final.contigs_virus.fna | FASTA file containing viral and pro-viral contigs |
|Project/        | final.contigs_virus_summary.tsv | Tab-separated summary file for viral contigs |
|Project/        | geNomad_run.stderr | Standard error log file from the geNomad run |
|Project/prokka  |                    | Files with prokka data used by the geNomad classifers |
|Project/Viral_Annotation|   contig_name/ | One contig fasta and one annotation output directory |
|Project/Plasmid_Annotation| contig_name/ | One contig fasta and one annotation output directory |
|Project/        | Analysis_Summary.html | HTML file of the output |


The html output file provides linkouts to the genome landing pages for the viral and plasmid contigs that were anntotated.  It also provides basic stistics and other helpful information.  

It is worth noting that the geNomad output files are greatly reduced compared to all of the output the geNomad command line tool produces.  We could add back more of the results files if it is justifiable.  

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).



| Script name | Purpose |
| ----------- | ------- |
| [App-Genomad.pl](service-scripts/App-Genomad.pl) | App script for the Genomad Service|

`p3_assembly` and `p3_genome_annotation` must be checked out in the environment in order for App-Genomad.pl to find their respective app-specs.


## See Also

## References

Identification of mobile genetic elements with geNomad Camargo, A.P., Roux, S., Schultz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C., - Nature Biotechnology (2023), DOI: 10.1038./s41587-023-01953-y.

https://portal.nersc.gov/geNomad/

https://github.com/apcamargo/geNomad?tab=readme-ov-file

