"""=================================================================================================
megan_taxonomy:make_megandb.py

Add up to to date accession -> taxonomy information to an existing megan_mab.db

db is an sqlite3 database with two tables, accessions and mappings
    accessions simply lists the names of the columns in mappings
    mappings
        Accession
        Taxonomy
        GTDB
        EGGNOG
        INTERPRO2GO
        SEED
        ED

 Mapping for prtein accession to taxonomy is in taxonomy db prot.accession2taxid.gz
 https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

 current megan-map.db file downloaded from https://software-ab.cs.uni-tuebingen.de/download/megan6

solution may be to use the polars library. plan
read nr and get list of accessions since the number in nr is far less than in prot.accession2taxid
get corresponding accessions and taxid from prot.accession2taxid (expect about 500G)
filter out accessions alread in megan-map
add new accessions and taxid to mapping table

13 November 2023     gribskov
================================================================================================="""
import sqlite3 as sq
import polars as pl

# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    datadir = '/scratch/negishi/mgribsko/megan/data/'
    megan = 'megan-map-Feb2022.db'
    prot2acc = datadir + 'prot.accession2taxid.FULL'

    # read nr data into dataframe df_nr
    nrid_file = datadir + 'nr.id.txt'
    df_nr = pl.read_csv(
        source=nrid_file,
        has_header=False,
        n_rows=10,
        truncate_ragged_lines=True
        )
    df_nr = df_nr.with_columns(pl.col("column_1").str.strip_chars_start('>').str.extract(r"(\S+) ",1))
    df_nr = df_nr.rename({"column_1": "accession"})
    print(df_nr)

    # connect to megan-map.db using sqlite3, and read taxonomy ids into dataframe df_mdb
    mdbfile = datadir + megan
    mdb_connection = sq.connect(mdbfile)

    df_mdb = pl.read_database(
        query=f'SELECT Accession, Taxonomy From mappings LIMIT 100',
        connection=mdb_connection
        )
    print( df_mdb)

    mdb = mdb_connection.cursor()
    # sql = 'SELECT Accession FROM mappings LIMIT 1000'
    # result = mdb.execute(sql)
    # accession_current = result.fetchall()

    # Read protein.accession to taxonomy data into dataframe df_pt
    df_pt = pl.read_csv(
        source=prot2acc,
        has_header=True,
        separator='\t',
        n_rows=100
        )
    print(df_pt)

    # user1 = {"id": 100, "name": "Rumpelstiltskin", "dob": "12/12/12"}
    # c.execute("INSERT INTO users VALUES (:id, :name, :dob)", user1)

    exit(0)
