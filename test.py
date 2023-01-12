from Bio import Entrez
from Bio import SeqIO

Entrez.email = "qridderpl@gmail.com"

def extractor(genbank):
    for info in genbank.splitlines():
        if "LOCUS" in info:
            eiwit_len = " ".join(info.split()[2:3])
        if "/chromosome" in info:
            chrom_num = info.strip().split("=")[1].strip("\"")  
            # protein_id_info_dict[protein_id]["Chromosome"].add(chrom_num)

        if '/db_xref=\"GeneID:' in info:
            NCBI_gene_id = info.strip().split(":")[1].strip("\"")
            # protein_id_info_dict[protein_id]["NCBI gene id"].add(NCBI_gene_id)
            NCBI_gene_id_found = True

        if "/sex" in info:
            sex_ = info.strip().split("=")[1].strip("\"")  
            # protein_id_info_dict[protein_id]["Sex"].add(sex_)

        if "/gene" in info:
            gene_name = info.strip().split("=")[1].strip("\"") 
            # protein_id_info_dict[gene_name]["gene_name"].add(tissue_type)

        if "/tissue_type" in info:
            tissue_type = info.strip().split("=")[1].strip("\"")
            tissue_type = tissue_type.replace(";", ",")
            # protein_id_info_dict[protein_id]["Tissue type"].add(tissue_type)
        if "ACCESSION" in info:
            protein_id = info.split()[1].strip("\"") 
        if "/product" in info:
            protein = info.split("=")[1].strip().strip("\"")  
            # protein_id_info_dict[protein_id]["Protein"].add(protein)

        if "/calculated_mol_wt" in info:
            mol_wt = info.strip().split("=")[1].strip("\"")
            # protein_id_info_dict[protein_id]["mol weight"].add(mol_wt)

        if "/site_type" in info:
            site_type = info.strip().split("=")[1].strip("\"")
            # protein_id_info_dict[protein_id]["site type"].add(site_type)

    output.write(f"{NCBI_gene_id};{gene_name};{chrom_num};{protein_id};{protein};{eiwit_len};{sex_};{tissue_type}\n")

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
                print(protein_ids)
        # if len(protein_ids) > 1000:
        #     break
    print(protein_ids)







with open("output.txt", "a") as output:
    for records in records_list:
        records.pop()
        for bank in records:
            try:
                extractor(bank)
            except:
                continue
        
# # Parse the records to extract the information you are interested in
# protein_id_info_dict = {}
# for record in records:
#     protein_id = record['Id']
#     protein_id_info_dict[protein_id] = {}
#     for feature in record['Features']:
#         if feature['Type'] == 'gene':
#             protein_id_info_dict[protein_id]['gene_name'] = feature['Qualifiers']['gene'][0]
#         elif feature['Type'] == 'product':
#             protein_id_info_dict[protein_id]['product'] = feature['Qualifiers']['product'][0]
#         # ...

# # Print the information for each protein
# for protein_id, info in protein_id_info_dict.items():
#     print(f"Protein ID: {protein_id}")
#     print(f"Gene name: {info['gene_name']}")
#     print(f"Product: {info['product']}")