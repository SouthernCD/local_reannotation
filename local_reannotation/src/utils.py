# import pysam
from yxseq import read_fasta, read_gff_file, Gene
from yxutil import mkdir, pickle_dump, have_file


def read_sequence_from_file(file_path):
    seq_dict, seq_id_list = read_fasta(file_path)

    if len(seq_id_list) == 0:
        return ""
    else:
        return seq_dict[seq_id_list[0]].seq


def gff_preparser(gff_file, genome_file, output_dir):
    mkdir(output_dir)
    gene_dict = read_gff_file(gff_file)['gene']
    for gene_id in gene_dict:
        mRNA_pkl = "%s/%s.pkl" % (output_dir, gene_id)
        if have_file(mRNA_pkl):
            continue
        print(gene_id)
        gene = gene_dict[gene_id]
        gene = Gene(from_gf=gene)
        gene.build_gene_seq(genome_file)
        mRNA = gene.model_mRNA
        pickle_dump(mRNA, mRNA_pkl)


if __name__ == "__main__":
    gff_file = "Sbicolor_730_v5.1.gene_exons.gff3"
    genome_file = "Sbicolor_730_v5.0.fa"
    output_dir = "Sbicolor_730_v5.1.gene_exons.gff3.pkls"
    gff_preparser(gff_file, genome_file, output_dir)
