from Bio import Entrez
from Bio import SeqIO

Entrez.email = "qridderpl@gmail.com"


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
    with open("proteoom_ids.txt") as f:
        protein_ids = []
        for x in f:
            if "XP_" in x:
                protein_ids.append(x.strip())
                if len(protein_ids) > 9000:
                    handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="gp")
                    records_list.append(handle.read().split("//"))
                    protein_ids = []
                    print("Hey")
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

def main():
    records_list = record_getter()
    output_writer(records_list)

main()