import pandas as pd
import sqlite3


def read_chunks(filename, connection, table, columns, limit, separator='\t', chunksize=1000000):
    '''---------------------------------------------------------------------------------------------
    :param filename:    data file name
    :param connection:  database connection
    :param table:       database table to create and load into
    :param columns:     column names for table
    :param limit:       maximum rows of data to load
    :param separator:   data column separator
    :param chunksize:   data chunksize (number of rows to load per batch)
    :return:
    ---------------------------------------------------------------------------------------------'''
    df = pd.read_csv(filename, sep=separator, header=None, names=columns, chunksize=chunksize)
    nchunk = 0
    for chunk in df:
        # Load the DataFrame into the database by chunks
        nchunk += 1
        print(f'processing chunk {nchunk}\t{nchunk * chunksize}')
        chunk.to_sql(table, conn, if_exists='append', index=False)

        # Commit the changes to the database
        conn.commit()
        if nchunk >= limit:
            break

    return nchunk


# ===================================================================================================
# main program
# ===================================================================================================
if __name__ == '__main__':
    # Create a connection to the SQLite3 database
    conn = sqlite3.connect('/scratch/negishi/mgribsko/megan/data/pandatest.db')
    pt = conn.cursor()

    # attach the megan map database and copy the mappings table
    pt.execute("ATTACH DATABASE '/scratch/negishi/mgribsko/megan/data/megan-map-work.db' as old")
    sql = '''
    CREATE TABLE IF NOT EXISTS mappings (
        Accession PRIMARY KEY , 
        Taxonomy INT, 
        GTDB INT, 
        EGGNOG INT, 
        INTERPRO2GO INT, 
        SEED INT, 
        EC INT
    ) WITHOUT ROWID
    '''
    pt.execute(sql)

    sql = '''
    INSERT INTO mappings
    SELECT * FROM old.mappings
    '''
    pt.execute(sql)

    # read in the new data
    pt.execute('DROP TABLE IF EXISTS mappings')
    pt.execute('CREATE TABLE IF NOT EXISTS mappings (Accession PRIMARY KEY, Taxonomy INT)')
    prot2acc = '/scratch/negishi/mgribsko/megan/data/prot.accession.notitle'
    read_chunks(prot2acc, conn, 'mappings', ['Accession', 'Taxonomy'], limit=1, )

    pt.execute('DROP TABLE IF EXISTS pro_tax')
    pt.execute('CREATE TABLE IF NOT EXISTS pro_tax (Accession PRIMARY KEY, Taxonomy INT)')

    pt.execute('DROP TABLE IF EXISTS nr_id')
    pt.execute('CREATE TABLE IF NOT EXISTS nr_id (Accession PRIMARY KEY)')

    # Read the CSV file into a Pandas DataFrame
    prot2acc = '/scratch/negishi/mgribsko/megan/data/prot.accession.notitle'
    read_chunks(prot2acc, conn, 'pro_tax', ['Accession', 'Taxonomy'], limit=5, )

    nr = '/scratch/negishi/mgribsko/megan/data/nr.id.ed2.txt'
    read_chunks(nr, conn, 'nr_id', ['Accession'], limit=5, )

    # find the set of protein id/taxonomy entries that are present in nr
    sql = '''
    CREATE TABLE nr_mappings AS
    SELECT p.Accession, p.Taxonomy FROM pro_tax as p, nr_id as n
    WHERE p.Accession=n.Accession;
    '''

    # Close the connection to the database
    conn.close()

# possibly useful sql
# CREATE TABLE nr_mappings AS
#     SELECT p.Accession, p.Taxonomy from pro_tax as m, nr_id as n
#     WHERE p.Accession=n.Accession;
#
# select count(*) from nrmappings n where not exists ( select * from mappings m where n.Accession=m.Accession);
#
# INSERT INTO mappings (Accession, Taxonomy)
# SELECT a.Accession, b.Taxonomy FROM mappings AS a INNER JOIN newmappings as b ON a.Accession=b.Accession
# ON CONFLICT(Accession) DO UPDATE SET Taxonomy=excluded.Taxonomy;
#
# insert into mappings (Accession, Taxonomy)
#    ...> select * from nrmappings n where not exists ( select * from mappings m where n.Accession=m.Accession)
#
#  insert into mappings (Accession,Taxonomy)
#    ...> select a.Accession,b.Taxonomy
#    ...> from mappings as a
#    ...> inner join nn as b
#    ...> on a.Accession=b.Accession
#    ...> on conflict(Accession) do update set Taxonomy=excluded.Taxonomy;
