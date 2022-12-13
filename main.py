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


class mRNA_py:
    def __init__(self, fasta_dict: dict, organism: str = None) -> None:
        self.fasta_dict = fasta_dict
        self.organism = organism or self.set_organism()

    def make_blastn_xml(self) -> None:
        input_str = (list(self.fasta_dict.values())[0])
        result_handle = NCBIWWW.qblast("blastn", "nt", input_str, hitlist_size=1)
        with open("NCBIWWW_blastn_out.xml", "w") as out_handle: 
            out_handle.write(result_handle.read())
            out_handle.close()

    def set_organism(self) -> None:
        my_file = Path("./NCBIWWW_blastn_out.xml")
        if not my_file.exists():
            self.make_blastn_xml()
        with open("NCBIWWW_blastn_out.xml", "r") as f:
            xml_read_list = f.readlines()
        for line in xml_read_list:
            if "PREDICTED:" in line:
                bad_format_organism_str = line
        organism = " ".join(bad_format_organism_str.split()[1:3])
        return organism

    def get_organism(self) -> str:
        return self.organism

    def get_organism_proteoom(self):
        # 
        url = f'https://ftp.ensembl.org/pub/release-108/fasta/{(self.organism).replace(" ", "_").lower()}/pep/'
        # print(url)
        r = requests.get(url)
        # print(r.status_code())
        for line in r.text.split():
            if ".all.fa.gz" in line:
                sub_url = line.split('"')[1]

                print(sub_url)

        # bashCommand = f"{url}{sub_url}"
        # bashCommand = f"mkdir proteoom"
        # print(bashCommand)
        # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        # output, error = process.communicate()
        # print(output)
        return f"{url}{sub_url}"

        # return (r.text)

    def cg_fract_dict(self) -> dict:
        cg_fract_dict = {}
        for k, v in self.fasta_dict.items():
            cg_perc_dict[k] = (v.count("C") + v.count("G"))/len(v)
        # self.cg_fract_dict = cg_fract_dict
        return cg_fract_dict

    def seq_len_dict(self) -> dict:
        seq_len_dict = {}
        for k, v in self.fasta_dict.items():
            seq_len_dict[k] = len(v)
        # self.seq_len_dict = seq_len_dict
        return seq_len_dict

    def protein_dict(self) -> dict:
        protein_dict = {}
        for k, v in self.fasta_dict.items():
            # protein_dict[k] = v.translate()
            # print(v.translate(to_stop=True))
            # print(k, v.translate())
            # not sure if it makes sense to check different reading frames
            frame_list = [v.translate(), v[1:].translate(), v[2:].translate()]
            frame_list.sort(key=lambda x:x.count("*"))
            best_frame = frame_list[0]
            split_frame = best_frame.split("*")
            split_frame.sort(key=lambda x:len(x), reverse=True)
            # print(k, split_frame[0])
            protein_dict[k] = split_frame[0]
            # print(f"{k}:\n {v.translate()}\n {v[1:].translate()}\n {v[2:].translate()}\n\n")
        return protein_dict

    def match_mrna_protein(self) -> dict:
        protein_dict = self.protein_dict()
        subprocess.call(f'bash alter_fasta.sh {"_".join(self.organism.split(" "))}', shell=True)
        print("files were created")
        # with open(f'{"_".join(self.organism.split(" "))}.ensembl.fa', 'r') as f:
        with open(f'{"_".join(self.organism.split(" "))}.uniprot.fa', 'r') as f:
            file_list = f.read().split("\n")
        # print(file_list)
        match_dict = {}
        for k, v in protein_dict.items():
            for line in file_list:
                if str(v) in line:
                    match_dict[k] = line
        return match_dict
        # prot = rna.translate(to_stop=True)

    def gene_id_to_proteins(self) -> dict:
        # gene_id_info_dict = defaultdict(dict)
        gene_id_info_dict = defaultdict(lambda: defaultdict(set))
        # todo: getting gene ID's needs to be incorperated into main
        with open("gene_ids", "r") as f:
            gene_id_list = f.read().split("\n")
        for gene_id in gene_id_list:
            # might raise errors?:
            handle = Entrez.efetch(db="protein", id=gene_id, rettype="gp")
            geneinfo_list = handle.read().split("\n")
            NCBI_gene_id_found = False
            for info in geneinfo_list:

                if "chromosome" in info:
                    chrom_num = info.strip().split("=")[1].strip("\"")  
                    gene_id_info_dict[gene_id]["Chromosome"].add(chrom_num)

                if "/db_xref=\"GeneID:" in info:
                    NCBI_gene_id = info.strip().split(":")[1].strip("\"")
                    gene_id_info_dict[gene_id]["NCBI gene id"].add(NCBI_gene_id)
                    NCBI_gene_id_found = True

                if "sex" in info:
                    sex_ = info.strip().split("=")[1].strip("\"")  
                    gene_id_info_dict[gene_id]["Sex"].add(sex_)

                if "tissue_type" in info:
                    tissue_type = info.strip().split("=")[1].strip("\"")
                    gene_id_info_dict[gene_id]["Tissue type"].add(tissue_type)

                if "/product" in info:
                    protein = info.split("=")[1].strip().strip("\"")  
                    gene_id_info_dict[gene_id]["Protein"].add(protein)

                if "/calculated_mol_wt" in info:
                    mol_wt = info.strip().split("=")[1].strip("\"")
                    gene_id_info_dict[gene_id]["mol weight"].add(mol_wt)

                if "/site_type" in info:
                    site_type = info.strip().split("=")[1].strip("\"")
                    gene_id_info_dict[gene_id]["site type"].add(site_type)




            print(gene_id_info_dict[gene_id])

            # maybe create a new method for this?
            while NCBI_gene_id_found:
                try:
                    handle = Entrez.efetch(db="gene", id=NCBI_gene_id, rettype="gene_table", format="xml")
                    record = Entrez.read(handle)
                    NCBI_gene_id_found = 0
                except Exception as e:
                    print(e)
                    continue

            for x in record:
                try:
                    for i in x["Entrezgene_comments"][1]["Gene-commentary_comment"][0]["Gene-commentary_products"][0]["Gene-commentary_products"]:
                        isoform = i["Gene-commentary_products"][0]["Gene-commentary_source"][0]["Other-source_post-text"]
                        print(isoform)
                        gene_id_info_dict[gene_id]["Protein_isoforms"].add(isoform)
                except:
                    for i in x["Entrezgene_comments"][2]["Gene-commentary_comment"][0]["Gene-commentary_products"][0]["Gene-commentary_products"]:
                        isoform = i["Gene-commentary_products"][0]["Gene-commentary_source"][0]["Other-source_post-text"]
                        print(isoform)
                        gene_id_info_dict[gene_id]["Protein_isoforms"].add(isoform)


        # print(gene_id_info_dict)
        # print(gene_id_info_dict.keys())
        # print(gene_id_info_dict["XP_007655542"])
        print(json.dumps(gene_id_info_dict, indent=2, default=tuple))
        return gene_id_info_dict




def read_fasta_file(file_name: str) -> dict:
    fasta_dict = {}
    with open(file_name) as f_read:
        for record in SeqIO.parse(f_read, "fasta"):
            fasta_dict[record.id] = record.seq
    return fasta_dict

def main():
    fasta_dict = read_fasta_file("T1.fa")
    # fasta_dict = read_fasta_file("T1.fa")
    inf1 = mRNA_py(fasta_dict)
    print(inf1.get_organism())
    print(inf1.get_organism_proteoom())
    # print(inf1.protein_dict())
    print(inf1.match_mrna_protein())
    inf1.gene_id_to_proteins()

    # print(inf1.gene_id_to_proteins())
    

if __name__ == '__main__':
    main()

