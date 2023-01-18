import psycopg2
import datetime
import requests
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

Entrez.email = "qridderpl@gmail.com"


def database_maken():
    conn_string = """host='145.97.18.224' dbname='s1136289_db' user='s1136289' password='s1136289'"""
    connection = psycopg2.connect(conn_string)
    print("Database maken en ermee verbinden is succesvol")
    return connection


def extractor(genbank):
    for info in genbank.splitlines():
        if "LOCUS" in info:
            eiwit_len = " ".join(info.split()[2:3])

        if "/chromosome" in info:
            chrom_num = info.strip().split("=")[1].strip("\"")

        if '/db_xref=\"GeneID:' in info:
            NCBI_gene_id = info.strip().split(":")[1].strip("\"")
            NCBI_gene_id_found = True

        if "/sex" in info:
            sex_ = info.strip().split("=")[1].strip("\"")

        if "/gene" in info:
            gene_name = info.strip().split("=")[1].strip("\"")

        if "/tissue_type" in info:
            tissue_type = info.strip().split("=")[1].strip("\"")
            tissue_type = tissue_type.replace(";", ",")

        if "ACCESSION" in info:
            protein_id = info.split()[1].strip("\"")
        if "/product" in info:
            protein = info.split("=")[1].strip().strip("\"")

        if "/calculated_mol_wt" in info:
            mol_wt = info.strip().split("=")[1].strip("\"")

        if "/site_type" in info:
            site_type = info.strip().split("=")[1].strip("\"")
    with open("output.txt", "a") as output:
        output.write(f"{NCBI_gene_id};{gene_name};{chrom_num};{protein_id};{protein};{eiwit_len};{sex_};{tissue_type}\n")


def record_getter():
    records_list = []
    count = 0
    with open("proteoom_ids.txt") as f:
        f = f.readlines()
        protein_ids = []
        for x in f:
            if "XP_" in x:
                protein_ids.append(x.strip())
                if len(protein_ids) > 9000:
                    count += 1
                    handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="gp")
                    records_list.append(handle.read().split("//"))
                    protein_ids = []
                    print(f"{count}/{len(f)/9000}")
                    break
        handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="gp")
        records_list.append(handle.read().split("//"))
    return records_list


def output_writer(records_list):
    for records in records_list:
        records.pop()
        for bank in records:
            try:
                extractor(bank)
            except:
                continue


def ec_getter(connection, row, gen_done):
    cur = connection.cursor()
    flag = False
    pathway = ""
    if row[0] not in gen_done:
        content = requests.get(f"http://rest.genome.jp/get/oaa:{row[0]}")
        for line in content.text.split("\n"):
            if "MODULE" in line or "BRITE" in line:
                flag = False
                pathway = ""
            if flag is True:
                pathway = " ".join(line.split()[1:])
            if "PATHWAY" in line:
                flag = True
                pathway = " ".join(line.split()[2:])
            if len(pathway) > 2:
                cur.execute("""
                    INSERT INTO pathway_gen_id ("gen_id", pathway_name) 
                    VALUES (%s, %s);
                    """,
                    (row[0], pathway))
            if "[EC:" in line:
                ecs = line.split("[EC:")[1].split()
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
    gen_done.add(row[0])                 
    return gen_done


def appender(connection):
    cur = connection.cursor()
    gen_done = set()
    with open("output.txt", "r") as f:
        f = f.readlines()
        for count ,row in enumerate(f, 1):
            row = row.split(";")
            cur.execute("""
            INSERT INTO gen ("gen_id", "gen_name") 
            VALUES (%s, %s) ON CONFLICT (gen_id) DO NOTHING;
            """,
            (row[0], row[1]))

            cur.execute("""
            INSERT INTO eiwit (eiwit_id, eiwit_naam, gen_id, chromosome, eiwit_lengte, sex_org, location_in_org) 
            VALUES (%s, %s, %s, %s, %s, %s, %s) ON CONFLICT (eiwit_id) DO NOTHING;
            """,
            (row[3], row[4], row[0], row[2], row[5] ,row[6], row[7]))
            gen_done = ec_getter(connection, row, gen_done)
            print(f"{count}/{len(f)}")
    connection.commit()


def seq_appender(connection):
    cur = connection.cursor()
    sequentie_lijst = []
    count = 0
    with open("T1.fa") as sequenties, open("ensemble_ids.txt") as eiwit_id:
        eiwit = eiwit_id.readlines()
        for seq in sequenties:
            if len(seq) > 60:
                count += 1
                cur.execute("""
            INSERT INTO mrna_brokstukken ("seq_id", "eiwit_id", sequentie_len ,sequentie) 
            VALUES (%s, %s, %s, %s) ON CONFLICT (seq_id) DO NOTHING;
            """,
            (count, eiwit[count-1].split(".")[0], len(seq), seq))


def blast_seq():
    with open("T1.fa", "r") as f:
        f = f.readlines()
        seq = f"{f[1]}"

    result_handle = NCBIWWW.qblast(program="blastn", database="nr" ,sequence=seq)
    blast_records = NCBIXML.read(result_handle)
    for count,x in enumerate(blast_records.descriptions):
        if count == 1:
            break
        gi = x.title.split("|")[3]

    handle = Entrez.efetch(db="nucleotide", id=f"{gi}", rettype="gp")
    for x in handle.readlines():
        if "ORGANISM" in x:
            organisme = " ".join(x.split()[1:])
            break
    return organisme


def assembly_getter(naam):
    getal = 0
    lijst = {}
    handle = Entrez.esearch(db="assembly", term=f"{naam}[Organism]")
    record = Entrez.read(handle)
    organisme = record["IdList"]
    for id in organisme:
        esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
        esummary_record = Entrez.read(esummary_handle)
        if len(esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']) > 0:
            getal += 1
            print(getal, esummary_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'], esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Organism'], esummary_record['DocumentSummarySet']['DocumentSummary'][0]['SubmissionDate'], esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq'])
            lijst[str(getal)] = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']

    if len(lijst) == 0:
        return False
    keuze = input(f"Welke assembly wil je gebruiken? (Typ een getal tussen 1-{getal}) ")

    if lijst.get(keuze)!=None:
        output = f"{lijst[keuze]}/{lijst[keuze].split('/')[-1]}_protein.faa.gz"
        with open("temp.txt", "w") as out:
                out.write(output)
        return True
    else:
        return False



dropper = """
    DROP TABLE IF EXISTS proteoom, eiwit, organism, gen, mrna_brokstukken, reactie ,ec_nummer_eiwit_id, pathway_gen_id CASCADE;
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
    chromosome VARCHAR(10),
    reactie_id VARCHAR(200),
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
    reactie_naam VARCHAR(200),
    substraat VARCHAR(50),
    product VARCHAR(50)
    );
    CREATE TABLE ec_nummer_eiwit_id(
    ec_nummer VARCHAR(50) REFERENCES reactie(ec_nummer),
    eiwit_id VARCHAR(20) REFERENCES eiwit(eiwit_id)
    );
    CREATE TABLE pathway_gen_id(
    pathway_name VARCHAR(200),
    gen_id VARCHAR(20) REFERENCES gen(gen_id)
    );"""


def main():
    print(sys.argv[1])
    if sys.argv[1] == "1":
        flag = False
        while flag == False:
            naam = input("Van welk organisme wil je de assembly? Druk enter als je niet weet van welk organisme de brokstukken afkomstig van zijn. ")
            if naam == "":
                print("BLASTen van 1 van de sequenties om te kijken om welk organisme het gaat. ")
                flag = assembly_getter(blast_seq())
            else:
                flag = assembly_getter(naam)
    if sys.argv[1] == "2":
        records_list = record_getter()
        output_writer(records_list)
        connection = database_maken()
        cursor = connection.cursor()
        cursor.execute(dropper)
        cursor.execute(sql)
        appender(connection)
        seq_appender(connection)
        connection.commit()
        connection.close()


main()
