# MEGAN 6 Taxonomy
create megan .tre and .map (and maybe map.db) files from NCBI taxonomy dump files, should 
correspond to current nr files

The .tre file is a newick formatted tree with the taxa represented by the taxid from the NCBI 
tasonomy database. The .map file lists the numeric taxid, scientific name, a -1, and the 
taxonomic level: 0-5 and 98-101 corresponding to Root (and unassigned), Kingdom, Phylum, Class, 
Order, Family, genus (98), group, genus+species, and subspecies, respectively. This translation 
is determined by inspection so it may not be completely reliable.

## proposed procedure
- create tree from taxid and parent class id in nodes.dmp
- add scientific names from names.dmp

### step 1
- read in nodes.dmp and see what ranks are actually present so that they can be mapped on the megan 
  lvls