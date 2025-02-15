from yxutil import cmd_run, have_file, rmdir
from yxseq import read_gff_file, Gene
import numpy as np


def get_all_cds_depth(bam_file, gff_file, depth_file):
    cds_bed_file = gff_file + ".cds.bed"
    if not have_file(cds_bed_file):
        gene_dict = read_gff_file(gff_file)['gene']
        cds_list = []
        for gene_id in gene_dict:
            gene = gene_dict[gene_id]
            gene = Gene(from_gf=gene)
            # gene.build_gene_seq(genome_file)
            cds_list += [cds for cds in gene.sub_features[0].sub_features if cds.type == 'CDS']

        sorted_cds_list = sorted(
            cds_list, key=lambda x: (x.chr_id, x.start, x.end))

        with open(cds_bed_file, 'w') as f:
            for cds in sorted_cds_list:
                f.write("%s\t%s\t%s\n" % (cds.chr_id, cds.start, cds.end))

    if have_file(depth_file):
        rmdir(depth_file)

    cmd_string = "samtools depth -Q 30 -aa -b %s %s > %s" % (
        cds_bed_file, bam_file, depth_file)
    cmd_run(cmd_string)

    return depth_file


def get_mRNA_depth(bam_file, mRNA, depth_file):
    cds_list = [cds for cds in mRNA.sub_features if cds.type == 'CDS']

    if have_file(depth_file):
        rmdir(depth_file)

    for cds in cds_list:
        cmd_string = "samtools depth -Q 30 -aa -r %s:%s-%s %s >> %s" % (
            mRNA.chr_id, cds.start, cds.end, bam_file, depth_file)
        cmd_run(cmd_string)

    return depth_file


def parse_depth_file(depth_file):
    with open(depth_file, 'r') as f:
        depth_list = [int(line.strip().split('\t')[-1]) for line in f]
    depth_list = np.array(depth_list)
    coverage = np.sum(depth_list > 0) / len(depth_list)
    depth = np.mean(depth_list)
    depth_sd = np.std(depth_list)
    return {
        'coverage': coverage,
        'depth': depth,
        'depth_sd': depth_sd
    }


if __name__ == '__main__':
    genome_file = 'Sbicolor_730_v5.0.fa'
    gff_file = 'Sbicolor_730_v5.1.gene_exons.gff3'
    bam_file = "IPDE.sorted.markdup.bam"
    gene_id = 'Sobic.005G213600.v5.1'
    work_dir = 'test'

    from yxseq import read_gff_file, Gene
    from yxutil import cmd_run, have_file, rmdir, pickle_dump

    gene_dict = read_gff_file(gff_file)['gene']
    gene = gene_dict[gene_id]
    gene = Gene(from_gf=gene)
    gene.build_gene_seq(genome_file)
    mRNA = gene.model_mRNA

    depth_file = "%s/depth.txt" % work_dir
    get_mRNA_depth(bam_file, mRNA, depth_file)
    parse_depth_file(depth_file)

    depth_file = bam_file + ".cds.depth"
    get_all_cds_depth(bam_file, gff_file, depth_file)
    parse_depth_file(depth_file)
