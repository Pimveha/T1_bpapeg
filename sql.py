import psycopg2
import datetime
import requests

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


def ec_getter(connection, row):
    cur = connection.cursor()
    content = requests.get(f"http://rest.genome.jp/get/oaa:{row[0]}")
    for line in content.text.split("\n"):
        if "[EC:" in line:
            ecs = line.split("[")[1].split()
            name = " ".join(line.split("[EC:")[0].split()[2:])
            for ec in ecs:
                # ec_lijst.append([ec.strip("EC:").strip("]"), row.split(';')[3]])
                cur.execute("""
                    INSERT INTO reactie ("ec_nummer", reactie_naam) 
                    VALUES (%s, %s) ON CONFLICT (ec_nummer) DO NOTHING;
                    """,
                    (ec.strip("EC:").strip("]"), name))
                
                cur.execute("""
                    INSERT INTO ec_nummer_eiwit_id ("ec_nummer", "eiwit_id") 
                    VALUES (%s, %s) ;
                    """,
                    (ec.strip("EC:").strip("]"), row[3]))
        if "NCBI-GeneID" in line:
            print(line.split(":")[1].strip())
        if "NCBI-ProteinID" in line:
            print(line.split(":")[1].strip())

def appender(connection):
    cur = connection.cursor()
    with open("output.txt", "r") as f:
        for count ,row in enumerate(f.readlines(), 1):
            row = row.split(";")

            cur.execute("""
            INSERT INTO gen ("gen_id", "gen_name") 
            VALUES (%s, %s) ON CONFLICT (gen_id) DO NOTHING;
            """,
            (row[0], row[1]))

            cur.execute("""
            INSERT INTO eiwit ("eiwit_id", eiwit_naam, gen_id, chromosome, eiwit_lengte) 
            VALUES (%s, %s, %s, %s, %s) ON CONFLICT (eiwit_id) DO NOTHING;
            """,
            (row[3], row[4], row[0], row[2], row[5]))   

            cur.execute("""
            INSERT INTO mrna_brokstukken ("seq_id", "eiwit_id", sequentie_len ,sequentie) 
            VALUES (%s, %s, %s, %s) ON CONFLICT (seq_id) DO NOTHING;
            """,
            (count, row[3], int(row[6]), row[7].strip()))   

            ec_getter(connection, row)

    connection.commit()

# def appender2(connection):
#     cur = connection.cursor()
#     with open("output.txt", "r") as f:
#         for count ,row in enumerate(f.readlines()):
                      
#     connection.commit()

dropper = """
    DROP TABLE IF EXISTS proteoom, eiwit, organism, gen, mrna_brokstukken, reactie ,ec_nummer_eiwit_id CASCADE;
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
    gen_id VARCHAR(50),
    ec_nummer VARCHAR(10),
    chromosome INT,
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
    eiwit_id VARCHAR(50) REFERENCES eiwit(eiwit_id),
    sequentie TEXT,
    sequentie_len INT,
    afgelezen_frame VARCHAR(20),
    gc_waarde VARCHAR(10)
    );
    CREATE TABLE reactie(
    ec_nummer VARCHAR(50) PRIMARY KEY,
    reactie_id VARCHAR(50),
    reactie_naam VARCHAR(50),
    substraat VARCHAR(50),
    product VARCHAR(50)
    );
    CREATE TABLE ec_nummer_eiwit_id(
    ec_nummer VARCHAR(50) REFERENCES reactie(ec_nummer),
    eiwit_id VARCHAR(20) REFERENCES eiwit(eiwit_id)
    );"""

def main():
    
    connection = database_maken()
    cursor = connection.cursor()
    cursor.execute(dropper)
    cursor.execute(sql)
    appender(connection)
    connection.commit()
    connection.close()

main()