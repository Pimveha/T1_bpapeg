import psycopg2
import datetime

conn_string = """
host='145.97.18.224' dbname='s1136289_db'
user='s1136289' password='s1136289'
"""
conn = psycopg2.connect(conn_string)
cursor = conn.cursor()

sql = """CREATE TABLE proteoom(
    proteoom_ID VARCHAR(50) PRIMARY KEY,
    lengte INT(50),
    aantal_genen INT(50),
    database_oorsprong VARCHAR(50)
    );
    CREATE TABLE eiwit(
    eiwit_ID VARCHAR(50) PRIMARY KEY,
    eiwit_naam VARCHAR(50),
    eiwit_lengte INT(20),
    EC_nummer VARCHAR(10),
    reactie_ID VARCHAR(50),
    reactie_naam VARCHAR(50),
    proteoom_ID VARCHAR(50) REFERENCES proteoom(proteoom_ID)
    substrate VARCHAR(50),
    product VARCHAR(50),
    location_in_org VARCHAR(50),
    sex_org VARCHAR(50)
    );
    CREATE TABLE organism(
    organism_ID VARCHAR(50) PRIMARY KEY,
    scientific_name VARCHAR(50),
    common_name_EN VARCHAR(50),
    common_name_NL VARCHAR(50),
    aantal_chromosomen INT(3)
    );
    CREATE TABLE gen(
    genID VARCHAR(50) PRIMARY KEY,
    gen_name VARCHAR(50),
    locatie VARCHAR(50),
    position_exon VARCHAR(50),
    position_intron VARCHAR(50),
    organism VARCHAR(50) REFERENCES organism(organism_ID)
    );
    CREATE TABLE mRNA_brokstukken(
    seqID VARCHAR(50) PRIMARY KEY,
    genID VARCHAR(50) REFERENCES gen(genID),
    lengte INT(50),
    seq TEXT,
    afgelezen_frame VARCHAR(20),
    gc_waarde VARCHAR(10)
    );
    CREATE TABLE genID_eiwitID(
    genID VARCHAR(20) PRIMARY KEY REFERENCES gen(genID),
    eiwitID VARCHAR(20) PRIMARY KEY REFERENCES eiwit(eiwitID)
    );"""

cursor.execute(sql)

conn.commit()
conn.close()

