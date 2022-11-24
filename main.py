from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import requests
import os
import subprocess
from pathlib import Path

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

    def set_cg_fract_dict(self) -> None:
        cg_fract_dict = {}
        for k, v in self.fasta_dict.items():
            cg_perc_dict[k] = (v.count("C") + v.count("G"))/len(v)
        self.cg_fract_dict = cg_fract_dict

    def get_cg_fract_dict(self) -> dict:
        return self.cg_fract_dict

    def set_seq_len_dict(self) -> None:
        seq_len_dict = {}
        for k, v in self.fasta_dict.items():
            seq_len_dict[k] = len(v)
        self.seq_len_dict = seq_len_dict

    def get_seq_len_dict(self) -> dict:
        return self.seq_len_dict


def read_fasta_file(file_name: str) -> dict:
    fasta_dict = {}
    with open(file_name) as f_read:
        for record in SeqIO.parse(f_read, "fasta"):
            fasta_dict[record.id] = record.seq
    return fasta_dict

def main():
    fasta_dict = read_fasta_file("T1.fa")
    inf1 = mRNA_py(fasta_dict)
    print(inf1.get_organism())
    

if __name__ == '__main__':
    main()
