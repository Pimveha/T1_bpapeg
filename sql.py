import os
import psycopg2
import requests
import sys
import time
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

Entrez.email = "qridderpl@gmail.com"


def database_maken():
    """
    Connectie maken met de database
    """
    conn_string = """host='145.97.18.224' dbname='c' 
    user='s1135502' password='s1135502'"""
    connection = psycopg2.connect(conn_string)
    print("Database maken en ermee verbinden is succesvol")
    return connection


def extractor(genbank):
    """
    Functie om de informatie van de protein pagina te verzamelen.
    """
    for info in genbank.splitlines():
        if "LOCUS" in info:
            eiwit_len = " ".join(info.split()[2:3])

        if "/chromosome" in info:
            chrom_num = info.strip().split("=")[1].strip("\"")

        if '/db_xref=\"GeneID:' in info:
            ncbi_gene_id = info.strip().split(":")[1].strip("\"")

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
    with open("output.txt", "a") as output:
        output.write(
            f"{ncbi_gene_id};{gene_name};{chrom_num};{protein_id};{protein};"
            f"{eiwit_len};{sex_};{tissue_type}\n")


def record_getter():
    """
    Deze functie haalt alle webpaginas van het gehele proteoom op en returned
    alle pagina's in een grote lijst.
    """
    records_list = []
    count = 0
    with open("proteoom_ids.txt") as f:
        f = f.readlines()
        protein_ids = []
        for x in f:
            if "XP_" in x:
                protein_ids.append(x.strip())
                # De eiwitten worden per 9000 aangevraagd.
                # Niet alle eiwitten kunnen tegelijk opgehaald worden.
                if len(protein_ids) > 9000:
                    count += 1
                    handle = Entrez.efetch(
                        db="protein", id=",".join(protein_ids), rettype="gp")
                    records_list.append(handle.read().split("//"))
                    protein_ids = []
                    print(f"{count}/{round(len(f)/9000)}")
        handle = Entrez.efetch(
            db="protein", id=",".join(protein_ids), rettype="gp")
        records_list.append(handle.read().split("//"))
    return records_list


def output_writer(records_list):
    """
    Loopt door de grote lijst met informatie en 
    geeft ze door aan de extractor functie.
    """
    for records in records_list:
        records.pop()
        for bank in records:
            try:
                extractor(bank)
            except:
                continue


def ec_getter(connection, row, gen_done):
    """
    Haalt de informatie op over de EC nummers en de pathways.
    Insert ze daarna in de SQL database
    """
    cur = connection.cursor()
    flag = False
    pathway = ""
    if row[0] not in gen_done:
        # De KEGG API wordt gebruikt om informatie per gen op te halen.
        content = requests.get(f"http://rest.genome.jp/get/oaa:{row[0]}", timeout=5)
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
    """
    Leest de informatie uit de output.txt file. 
    Daarna INSERT hij het in de SQL database.
    """
    start_time = time.time()
    cur = connection.cursor()
    gen_done = set()
    with open("output.txt", "r") as f:
        f = f.readlines()
        for count, row in enumerate(f, 1):
            # Om de gemiddelde tijd te berekenen.
            current_time = time.time()
            elapsed_time = current_time - start_time
            remaining_iterations = len(f) - count
            estimated_remaining_time = remaining_iterations \
                * (elapsed_time / count)
            remaining_time_h = int(estimated_remaining_time // 3600)
            remaining_time_m = int((estimated_remaining_time % 3600) // 60)
            remaining_time_s = int((estimated_remaining_time % 3600) % 60)

            row = row.split(";")

            cur.execute("""
            INSERT INTO gen ("gen_id", "gen_name") 
            VALUES (%s, %s) ON CONFLICT (gen_id) DO NOTHING;
            """,
                        (row[0], row[1]))

            cur.execute("""
            INSERT INTO eiwit (eiwit_id, eiwit_naam, gen_id, chromosome, 
            eiwit_lengte, sex_org, location_in_org) 
            VALUES (%s, %s, %s, %s, %s, %s, %s) 
            ON CONFLICT (eiwit_id) DO NOTHING;
            """,
                    (row[3], row[4], row[0], row[2], row[5], row[6], row[7]))
            gen_done = ec_getter(connection, row, gen_done)
            os.system("clear")
            print(f"{count}/{len(f)}.\nEstimated time remaining: "
                f"{str(remaining_time_h).zfill(2)}:"
                f"{str(remaining_time_m).zfill(2)}:"
                f"{str(remaining_time_s).zfill(2)}")

    connection.commit()


def seq_appender(connection):
    """
    Leest de sequenties van de brokstukken samen met de eiwit_id
    en zet deze in de SQL database
    """
    cur = connection.cursor()
    count = 0
    with open("T1.fa") as sequenties, open("totaal_sort.txt") as eiwit_id:
        eiwit = eiwit_id.readlines()
        for seq in sequenties:
            if len(seq) > 60:
                count += 1
                cur.execute("""
            INSERT INTO mrna_brokstukken ("seq_id", "eiwit_id", 
                sequentie_len, sequentie) 
            VALUES (%s, %s, %s, %s) ON CONFLICT (seq_id) DO NOTHING;
            """,
                        (count, eiwit[count-1].split(".")[0], len(seq), seq))


def blast_seq():
    """
    Als het organisme onbekend is wordt hier een van de sequenties geblast om het organisme te vinden.
    """
    with open("T1.fa", "r") as f:
        f = f.readlines()
        seq = f"{f[1]}"

    result_handle = NCBIWWW.qblast(
        program="blastn", database="nr", sequence=seq)
    blast_records = NCBIXML.read(result_handle)
    for count, x in enumerate(blast_records.descriptions):
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
    """
    Zoekt in de NCBI database naar het proteoom voor het opgegeven organisme.
    """
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
            print(getal, esummary_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'], esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Organism'],
                  esummary_record['DocumentSummarySet']['DocumentSummary'][0]['SubmissionDate'], esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq'])
            lijst[str(
                getal)] = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']

    if len(lijst) == 0:
        return False
    keuze = input(
        f"Welke assembly wil je gebruiken? (Typ een getal tussen 1-{getal}) ")

    if lijst.get(keuze) != None:
        output = f"{lijst[keuze]}/{lijst[keuze].split('/')[-1]}_protein.faa.gz"
        with open("temp.txt", "w") as out:
            out.write(output)
        return True
    else:
        return False


# Verwijdert de databases als deze al bestaan.
DROPPER = """
    DROP TABLE IF EXISTS proteoom, eiwit, organism, gen, mrna_brokstukken, 
        reactie ,ec_nummer_eiwit_id, pathway_gen_id CASCADE;
    """


# Hier worden alle tables aan gemaakt met de juiste connecties.
SQL = """CREATE TABLE proteoom(
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
    """
    Wordt aangeroepen via het download.sh script. 
    Eerst wordt het eerste gedeelte opgroepen door sys.argv 1 mee te geven. 
    Hier wordt het proteoom opgevraagd via assembly_getter.
    Hierna wordt sys.argv 2 meegegeven. 
    Hier wordt alle informatie opgehaald en in de tables gevuld.
    """
    if sys.argv[1] == "1":
        flag = False
        while flag is False:
            naam = input("Van welk organisme wil je de assembly? "
                "Druk enter als je niet weet van welk organisme de "
                "brokstukken afkomstig van zijn.")
            if naam == "":
                print("BLASTen van 1 van de sequenties om te "
                "kijken om welk organisme het gaat. ")
                flag = assembly_getter(blast_seq())
            else:
                flag = assembly_getter(naam)
    if sys.argv[1] == "2":
        print("Starten van het vullen van de database. ")
        records_list = record_getter()
        output_writer(records_list)
        connection = database_maken()
        cursor = connection.cursor()
        cursor.execute(DROPPER)
        cursor.execute(SQL)
        appender(connection)
        seq_appender(connection)
        connection.commit()
        connection.close()


main()
