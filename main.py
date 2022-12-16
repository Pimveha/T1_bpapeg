from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import requests
import os
import subprocess
from pathlib import Path
from collections import defaultdict
import json

#change to input
Entrez.email = input("mail: ")


def protein_id_to_info():
        # protein_id_info_dict = defaultdict(dict)
        protein_id_info_dict = defaultdict(lambda: defaultdict(set))
        # todo: getting gene ID's needs to be incorperated into main
        sequentie_lijst = []
        with open("T1.fa") as sequenties:
            for x in sequenties:
                if len(x) > 60:
                    sequentie_lijst.append(x)
        with open("totaal_sort.txt", "r") as f:
            protein_id_list = f.read().split("\n")
        with open("output.txt", "w") as output:
            for protein_id in protein_id_list:
                # might raise errors?:
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="gp")
                geneinfo_list = handle.read().split("\n")
                NCBI_gene_id_found = False
                for info in geneinfo_list:
                    if "LOCUS" in info:
                        eiwit_len = " ".join(info.split()[2:3])

                    if "chromosome" in info:
                        chrom_num = info.strip().split("=")[1].strip("\"")  
                        protein_id_info_dict[protein_id]["Chromosome"].add(chrom_num)

                    if "/db_xref=\"GeneID:" in info:
                        NCBI_gene_id = info.strip().split(":")[1].strip("\"")
                        protein_id_info_dict[protein_id]["NCBI gene id"].add(NCBI_gene_id)
                        NCBI_gene_id_found = True

                    if "sex" in info:
                        sex_ = info.strip().split("=")[1].strip("\"")  
                        protein_id_info_dict[protein_id]["Sex"].add(sex_)

                    if "/gene" in info:
                        gene_name = info.strip().split("=")[1].strip("\"") 
                        protein_id_info_dict[gene_name]["gene_name"].add(tissue_type)

                    if "tissue_type" in info:
                        tissue_type = info.strip().split("=")[1].strip("\"")
                        tissue_type = tissue_type.replace(";", ",")
                        protein_id_info_dict[protein_id]["Tissue type"].add(tissue_type)

                    if "/product" in info:
                        protein = info.split("=")[1].strip().strip("\"")  
                        protein_id_info_dict[protein_id]["Protein"].add(protein)

                    if "/calculated_mol_wt" in info:
                        mol_wt = info.strip().split("=")[1].strip("\"")
                        protein_id_info_dict[protein_id]["mol weight"].add(mol_wt)

                    if "/site_type" in info:
                        site_type = info.strip().split("=")[1].strip("\"")
                        protein_id_info_dict[protein_id]["site type"].add(site_type)
                output.write(f"{NCBI_gene_id};{gene_name};{chrom_num};{protein_id};{protein};{eiwit_len};{sex_};{tissue_type};{len(sequentie_lijst[0])};{sequentie_lijst.pop(0)}")




                print(protein_id_info_dict[protein_id])

                # maybe create a new method for this?

                while NCBI_gene_id_found:
                    try:
                        handle = Entrez.efetch(db="gene", id=NCBI_gene_id, rettype="gene_table", format="xml")
                        record = Entrez.read(handle)
                        break
                    except Exception as e:
                        print(e)
                        continue

                for x in record:
                    try:
                        for i in x["Entrezgene_comments"][1]["Gene-commentary_comment"][0]["Gene-commentary_products"][0]["Gene-commentary_products"]:
                            isoform = i["Gene-commentary_products"][0]["Gene-commentary_source"][0]["Other-source_post-text"]
                            print(isoform)
                            protein_id_info_dict[protein_id]["Protein_isoforms"].add(isoform)
                    except:
                        for i in x["Entrezgene_comments"][2]["Gene-commentary_comment"][0]["Gene-commentary_products"][0]["Gene-commentary_products"]:
                            isoform = i["Gene-commentary_products"][0]["Gene-commentary_source"][0]["Other-source_post-text"]
                            print(isoform)
                            protein_id_info_dict[protein_id]["Protein_isoforms"].add(isoform)


        # print(protein_id_info_dict)
        # print(protein_id_info_dict.keys())
        # print(protein_id_info_dict["XP_007655542"])
        # print(json.dumps(protein_id_info_dict, indent=2, default=tuple))
        return protein_id_info_dict


def main():
    # fasta_dict = read_fasta_file("T1.fa")
    # # fasta_dict = read_fasta_file("T1.fa")
    # inf1 = mRNA_py(fasta_dict)
    # print(inf1.get_organism())
    # print(inf1.get_organism_proteoom())
    # # print(inf1.protein_dict())
    # print(inf1.match_mrna_protein())
    protein_id_to_info()

    # print(inf1.protein_id_to_info())
    

if __name__ == '__main__':
    main()

