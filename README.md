# Rukki
Rukki ([“spinning wheel”](https://en.wikipedia.org/wiki/Spinning_wheel) in Finnish) is a companion tool of Verkko assembler for extracting longer ‘scaffold’ paths from assembly graphs.

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

Rust 2021 edition compiler.
Try building with `cargo build --release`.

## Usage

Basic usage example
```
./target/release/rukki trio -g <graph.gfa> -m <marker_cnts.tsv> -p <out_paths.tsv> [--final-assign <node_assign.tsv>] [--try-fill-bubbles]
```

* `graph.gfa` -- graph in GFA format. Sequences are ignored and optiona. 
Node coverage values will be used for various purposes if provided (as `RC:i:`, `FC:i:`, and/or `ll:f:` tags for `S` records).
* `marker_cnts.tsv` -- TSV file, where first three columns of every line are interpreted as
`node_name\tmaternal\tpaternal`, where 'maternal'/'paternal' are parental-specific marker counts.
All columns after the third in TSV are ignored.
* `out_paths.tsv` -- TSV output containing haplo-paths (one per line).
Lines have format `path_name\tpath\tassignment`.
By default paths are formatted as (`<node>[+-](,<node>[+-])*`).
Also supports GAF path format, i.e. `([<>]<node>)+`, via the `--gaf-format` option.
The path can also include gaps in the `[NXXXN]` format, where `XXX` is the integer giving an estimate gap size.
Estimators are currently work in progress and not available for all cases.
Default gap size (for cases where estimator is not yet available) is 5kb.
Minimal reported value is currently fixed at 1kb (if an estimated value is lower than 1kb, 1kb will be reported instead).
Gaps represent either an absense of the appropriate connections or a localized ambiguity within the graph.
Assignment categories are `MATERNAL`, `PATERNAL` or `NA` (for _unassigned_). 
`NA` can only be associated with paths consisting of a single node.
Every node of the graph is guaranteed to be covered by one or more output paths.
* `--try-fill-bubbles` -- enables more agressive filling of ambiguous regions with one of available alternatives (recommended).
* `node_assign.tsv` -- assignments of individual nodes, reflecting their usage by haplo-paths (`MATERNAL`, `PATERNAL` or `HOMOZYGOUS`). Nodes forming _unassigned_ paths are excluded.

To see all options use:
```
./target/release/rukki trio --help
```
