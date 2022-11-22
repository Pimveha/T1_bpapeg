from Bio import SeqIO
import requests

class mRNA_py:
    def __init__(self, fasta_dict: dict, organism: str = None) -> None:
        self.fasta_dict = fasta_dict
        self.organism = organism or set_organism()

    def set_organism(self) -> None:
        # check blast 100% alignment
        self.organism = "froggo" # temp

    def get_organism(self) -> str:
        return self.organism

    def set_cg_fract_dict(self) -> None:
        cg_fract_dict = {}
        for k, v in self.fasta_dict.items():
            cg_perc_dict[k] = (v.count("C") + v.count("G"))/len(v)
        self.cg_fract_dict = cg_fract_dict

    def get_cg_fract_dict(self) -> dict:
        return self.cg_fract_dict


def read_fasta_file(file_name: str) -> dict:
    fasta_dict = {}
    with open(file_name) as f_read:
        for record in SeqIO.parse(f_read, "fasta"):
            fasta_dict[record.id] = record.seq
    return fasta_dict

def main():
    fasta_dict = read_fasta_file("T1_c.fa")
    


if __name__ == '__main__':
    main()

