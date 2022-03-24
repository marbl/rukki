# Rukki
Rukki (“spinning wheel” in Finnish) is a tool companion tool of Verkko assembler for extracting longer ‘scaffold’ paths from assembly graphs.

It's primary purpose is to utilize parental (trio) information attributed to the Verkko assembly graph nodes for extraction of longer haplotypes reconstructions in diploid organisms.

Rukki first assigns parental (maternal/paternal) classes to the nodes with prevalence of corresponding parental-specific markers, tries to identify homozygous nodes (belonging to both haplotypes), and then performs heuristic search of haplotype-paths starting from long nodes of the graph.

Plans are to turn it into a tool for comprehensive analysis of assembly graphs, in particular support extraction of ‘primary’ and ‘alt’ scaffolds is under development. 

## Some useful features

* Can exclude suspicious nodes (having high prevalence of both types of markers) from the traversals. 
* Prevents the re-use of long nodes (unless assigned as homozygous), which can happen if the graph has missing connections.
* Can scaffold across gaps in one haplotype if the other haplotype is intact.
* Can scaffold across ambiguous regions (e.g. tandem repeat expansions).
* Can deal with ambiguous bubble structures (either scaffold across of force-pick a path).

## Requirements

TODO

## Usage

TODO
