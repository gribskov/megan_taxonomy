# MEGAN 6 Taxonomy
create from NCBI taxonomy dump files, should correspond to current nr files

## proposed procedure
- create tree from taxid and parent class id in nodes.dmp
- add scientific names from names.dmp

### step 1
- read in nodes.dmp and see what ranks are actually present so that they can be mapped on the megan 
  lvls