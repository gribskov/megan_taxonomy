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
    megan = 'megan-map-Feb2022.db'
    cstr = f'sqlite:{megan}'
    megan_new = 'megan_new.db'
    prot2acc = 'prot.accession2taxid.FULL'

    mdb_connection = sq.connect(megan)
    # df =pl.read_sql(f'SELECT Accesion, Taxonomy From mappings, cstr')
    # mdbold_connection = sq.connect(megan)
    df = pl.read_database(
        query=f'SELECT Accession, Taxonomy From mappings',
        connection = mdb_connection
        )
    # mdb_connection = sq.connect(megan_new)
    # with mdb_connection:
    #     mdbold_connection.backup(mdb_connection)
    # mdbold_connection.close()

    mdb = mdb_connection.cursor()
    sql = 'SELECT Accession FROM mappings LIMIT 1000'
    result = mdb.execute(sql)
    accession_current = result.fetchall()

    p2a = open(prot2acc, 'r')
    p2a.readline() # strip header
    accession_new = {}
    n_new = 0
    n = 0
    t = 0
    for line in p2a:
        t += 1
        acc, tax = line.rstrip().split()
        if tax == '0':
            continue
        if not n%1000000:
            print(f'{n:12d}\t{t:12d}\t{acc}\t{tax}')
        n += 1

        if acc in accession_current:
            n_new += 1
            accession_new[acc] = tax

    # user1 = {"id": 100, "name": "Rumpelstiltskin", "dob": "12/12/12"}
    # c.execute("INSERT INTO users VALUES (:id, :name, :dob)", user1)

    exit(0)
