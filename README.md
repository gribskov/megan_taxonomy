# MEGAN 6 Taxonomy
create megan .tre, .map, and map.db files from NCBI taxonomy dump files, should 
correspond to current nr files

The .tre file is a newick formatted tree with the taxa represented by the taxid from the NCBI 
tasonomy database. The .map file lists the numeric taxid, scientific name, a -1, and the 
taxonomic level: 0-5 and 98-101 corresponding to Root (and unassigned), Kingdom, Phylum, Class, 
Order, Family, Genus (98), group, genus+species, and subspecies, respectively. This translation 
is determined by inspection, so it may not be completely reliable. Here are the observed 
taxonomic ranks, their counts (June 2023), and their assignment to megan lvl.

| rank             |   count | lvl | &nbsp;  | rank                 |   count | lvl |
|------------------|--------:|----:|---------|----------------------|--------:|----:|
| no rank          |  237092 |   0 | &nbsp;  | **genus**            |  106865 |  98 |
| **superkingdom** |       4 |   1 | &nbsp;  | subgenus             |    1770 |  98 |
| kingdom          |      13 |   1 | &nbsp;  | series               |       9 |  98 |
| subkingdom       |       1 |   1 | &nbsp;  | section              |     532 |  98 |
| **phylum**       |     307 |   2 | &nbsp;  | subsection           |      41 |  98 |
| subphylum        |      31 |   2 | &nbsp;  | **species group**    |     358 |  99 |
| **superclass**   |       6 |   3 | &nbsp;  | species              | 2051563 |  99 | 
| class            |     505 |   3 | &nbsp;  | **species subgroup** |     134 | 100 |
| subclass         |     168 |   3 | &nbsp;  | subspecies           |   28292 | 100 |
| infraclass       |      19 |   3 | &nbsp;  | **varietas**         |    9606 | 101 |
| cohort           |       5 |   3 | &nbsp;  | forma                |     670 | 101 |
| subcohort        |       3 |   3 | &nbsp;  | forma specialis      |     750 | 101 | 
| **superorder**   |      57 |   4 | &nbsp;  | pathogroup           |       5 | 101 | 
| order            |    1833 |   4 | &nbsp;  | morph                |      12 | 101 |
| suborder         |     372 |   4 | &nbsp;  | biotype              |      17 | 101 |
| infraorder       |     133 |   4 | &nbsp;  | genotype             |      21 | 101 |
| parvorder        |      26 |   4 | &nbsp;  | serogroup            |     147 | 101 |
| **superfamily**  |     896 |   5 | &nbsp;  | clade                |     955 | 101 |
| family           |   10126 |   5 | &nbsp;  | serotype             |    1237 | 101 |
| subfamily        |    3237 |   5 | &nbsp;  | isolate              |    1322 | 101 |
| tribe            |    2342 |   5 | &nbsp;  | strain               |   45790 | 101 |
| subtribe         |     583 |   5 | &nbsp;  |                      |         |     |

## proposed procedure
- create tree from taxid and parent class id in nodes.dmp
- add scientific names from names.dmp
- write out .tre file (newick format)
- write out mapping file (name, taxid, -1, lvl)
- write new mapping.db file (sqlite3 database)

### To Do
- [x] read in nodes.dmp and see what ranks are actually present so that they can be mapped on the 
  megan lvl
- [x] write function to convert text rank to numerical rank
- [x] implement non-recursive newick function, recursive function from Tree class hits recursion 
  limit. Seems like it shouldn't, but the tree is 2.5 M taxa.
- [x] write out tree
- [x] write out map file
- [ ] read megan_db
- [ ] megan_db entries with new map information]