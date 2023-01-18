from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

Entrez.email = "qridderpl@gmail.com"


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


def main():
    flag = False
    while flag == False:
        naam = input("Van welk organisme wil je de assembly? Druk enter als je niet weet van welk organisme de brokstukken afkomstig van zijn. ")
        if naam == "":
            print("BLASTen van 1 van de sequenties om te kijken om welk organisme het gaat. ")
            flag = assembly_getter(blast_seq())
        else:
            flag = assembly_getter(naam)


main()
