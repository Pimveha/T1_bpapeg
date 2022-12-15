import psycopg2
import datetime

def database_maken():
    # DB_NAME = "bpapget1_db"  # name database
    # DB_USER = "bpapget1"  # name user
    # DB_PASS = "bpapget1"  # password
    # DB_HOST = "145.97.18.224"  # host
    # DB_PORT = "8081"  # default
    # connection = psycopg2.connect(database=DB_NAME, user=DB_USER,
    #                               password=DB_PASS, host=DB_HOST, port=DB_PORT)

    conn_string = """host='145.97.18.224' dbname='bpapget1_db' user='bpapget1' password='bpapget1'"""
    connection = psycopg2.connect(conn_string)

    print("Database maken en ermee verbinden is succesvol")
    return connection


dropper = """
    DROP TABLE IF EXISTS proteoom, eiwit, organism, gen, mrna_brokstukken, gen_id_eiwit_id CASCADE;
    """
    
sql = """CREATE TABLE proteoom(
    proteoom_ID VARCHAR(50) PRIMARY KEY,
    lengte INT,
    aantal_genen INT,
    database_oorsprong VARCHAR(50)
    );
    CREATE TABLE eiwit(
    eiwit_id VARCHAR(50) PRIMARY KEY,
    eiwit_naam VARCHAR(50),
    eiwit_lengte INT,
    ec_nummer VARCHAR(10),
    reactie_id VARCHAR(50),
    reactie_naam VARCHAR(50),
    proteoom_id VARCHAR(50) REFERENCES proteoom(proteoom_id),
    substrate VARCHAR(50),
    product VARCHAR(50),
    location_in_org VARCHAR(50),
    sex_org VARCHAR(50)
    );
    CREATE TABLE organism(
    organism_id VARCHAR(50) PRIMARY KEY,
    scientific_name VARCHAR(50),
    common_name_en VARCHAR(50),
    common_name_nl VARCHAR(50),
    aantal_chromosomen INT
    );
    CREATE TABLE gen(
    gen_id VARCHAR(50) PRIMARY KEY,
    gen_name VARCHAR(50),
    locatie VARCHAR(50),
    position_exon VARCHAR(50),
    position_intron VARCHAR(50),
    organism VARCHAR(50) REFERENCES organism(organism_id)
    );
    CREATE TABLE mrna_brokstukken(
    seq_id VARCHAR(50) PRIMARY KEY,
    gen_id VARCHAR(50) REFERENCES gen(gen_id),
    lengte INT,
    seq TEXT,
    afgelezen_frame VARCHAR(20),
    gc_waarde VARCHAR(10)
    );
    CREATE TABLE gen_id_eiwit_id(
    gen_id VARCHAR(20) REFERENCES gen(gen_id),
    eiwit_id VARCHAR(20) REFERENCES eiwit(eiwit_id)
    );"""

connection = database_maken()
cursor = connection.cursor()
cursor.execute(dropper)
cursor.execute(sql)

connection.commit()
connection.close()