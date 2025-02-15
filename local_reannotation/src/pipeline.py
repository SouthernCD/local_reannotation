from local_reannotation.src.utils import read_sequence_from_file
from local_reannotation.src.miniprot import miniprot_job
from local_reannotation.src.bamstat import get_mRNA_depth, parse_depth_file
from yxutil import mkdir, rmdir, pickle_load
from yxseq import read_fasta
import json
import shutil
import warnings
import os
warnings.filterwarnings('ignore')


def gen_annotation(gene_id, mRNA_pkl, work_dir, bam_file, haplo_seq_fa, ref_sequence_fa):
    work_dir = os.path.abspath(work_dir)
    mkdir(work_dir)
    
    seq_dict, seq_id_list = read_fasta(haplo_seq_fa)
    haplo_seq_list = [seq_dict[i].seq for i in seq_dict]
    ref_sequence = read_sequence_from_file(ref_sequence_fa)

    mRNA = pickle_load(mRNA_pkl)

    ref_aa_seq, ref_cds_seq = mRNA.aa_seq, mRNA.cds_seq

    # state reseq depth and coverage in the CDS region
    depth_file = "%s/cds_depth.txt" % work_dir
    get_mRNA_depth(bam_file, mRNA, depth_file)
    depth_report = parse_depth_file(depth_file)

    hap_gene_dict = {}
    for i, seq in enumerate(haplo_seq_list):
        hap_gene_dict[i] = "%s/haplotype%d_genome.fasta" % (
            work_dir, i+1)
        with open(hap_gene_dict[i], 'w') as f:
            f.write(f">haplotype{i+1}\n{seq}\n")

    ref_seq_file = "%s/reference_genome.fasta" % work_dir
    with open(ref_seq_file, 'w') as f:
        f.write(f">ref_genome\n{ref_sequence}\n")

    # clustalw nucleotide alignment
    genome_aln_file = "%s/all_genome.fasta" % work_dir
    with open(genome_aln_file, 'w') as f:
        f.write(f">ref_genome\n{ref_sequence}\n")
        for i, seq in enumerate(haplo_seq_list):
            f.write(f">haplotype{i+1}\n{seq}\n")

    # cmd_string = "clustalw2 -INFILE=%s -ALIGN -OUTPUT=FASTA -OUTFILE=%s.aln -type=DNA" % (
    #     genome_aln_file, genome_aln_file)
    # cmd_run(cmd_string, cwd=gene_dir)

    # re-annotation
    anno_work_dir = work_dir + "/annotation"
    mkdir(anno_work_dir)

    # write reference query protein sequence
    ref_prot_file = "%s/reference_protein.fasta" % anno_work_dir
    with open(ref_prot_file, 'w') as f:
        f.write(">%s\n%s\n" % (gene_id, ref_aa_seq))

    ref_cds_file = "%s/reference_CDS.fasta" % anno_work_dir
    with open(ref_cds_file, 'w') as f:
        f.write(">%s\n%s\n" % (gene_id, ref_cds_seq))

    # run miniprot
    reanno_ref_pseudo_dict, reanno_ref_pt_seq, reanno_ref_cds_seq = miniprot_job(
        ref_prot_file, ref_seq_file, anno_work_dir + "/reference_miniprot")

    reanno_ref_pt_file = "%s/reanno_ref_protein.fasta" % anno_work_dir
    with open(reanno_ref_pt_file, 'w') as f:
        f.write(">reanno_ref\n%s\n" % reanno_ref_pt_seq)

    reanno_ref_cds_file = "%s/reanno_ref_CDS.fasta" % anno_work_dir
    with open(reanno_ref_cds_file, 'w') as f:
        f.write(">reanno_ref\n%s\n" % reanno_ref_cds_seq)

    hap_miniport_results_dict = {}
    for i, hap_gene_file in hap_gene_dict.items():
        pseudo_dict, pt_seq, cds_seq = miniprot_job(
            ref_prot_file, hap_gene_file, anno_work_dir + f"/haplotype{i+1}_miniprot")
        hap_miniport_results_dict[i] = (pseudo_dict, pt_seq, cds_seq)

        hap_pt_file = anno_work_dir + f"/haplotype{i+1}_protein.fasta"
        with open(hap_pt_file, 'w') as f:
            f.write(f">haplotype{i+1}\n{pt_seq}\n")

        hap_cds_file = anno_work_dir + f"/haplotype{i+1}_CDS.fasta"
        with open(hap_cds_file, 'w') as f:
            f.write(f">haplotype{i+1}\n{cds_seq}\n")

    # run clustalw protein alignment
    all_protein_file = anno_work_dir + "/all_protein.fasta"
    with open(all_protein_file, 'w') as f:
        f.write(">%s\n%s\n" % (gene_id, ref_aa_seq))
        f.write(f">reanno_ref\n{reanno_ref_pt_seq}\n")
        for i, (_, pt_seq, _) in hap_miniport_results_dict.items():
            f.write(f">haplotype{i+1}\n{pt_seq}\n")

    # cmd_string = "clustalw2 -INFILE=%s -ALIGN -OUTPUT=FASTA -OUTFILE=%s.aln -type=PROTEIN" % (
    #     all_protein_file, all_protein_file)
    # cmd_string = "mafft --preservecase --auto %s > %s.aln" % (
    #     all_protein_file, all_protein_file)
    # cmd_run(cmd_string, cwd=anno_work_dir)

    # run clustalw cds alignment
    all_cds_file = anno_work_dir + "/all_CDS.fasta"
    with open(all_cds_file, 'w') as f:
        f.write(">%s\n%s\n" % (gene_id, ref_cds_seq))
        f.write(
            f">reanno_ref\n{read_sequence_from_file(reanno_ref_cds_file)}\n")
        for i, (_, _, cds_seq) in hap_miniport_results_dict.items():
            f.write(f">haplotype{i+1}\n{cds_seq}\n")

    # cmd_string = "clustalw2 -INFILE=%s -ALIGN -OUTPUT=FASTA -OUTFILE=%s.aln -type=DNA" % (
    #     all_cds_file, all_cds_file)
    # cmd_string = "mafft --preservecase --auto %s > %s.aln" % (
    #     all_cds_file, all_cds_file)
    # cmd_run(cmd_string, cwd=anno_work_dir)

    # merge all results
    results_dict = {'gene_id': gene_id, 'bam_file': bam_file}
    results_dict['reseq_stat'] = depth_report
    results_dict['ref_redo'] = reanno_ref_pseudo_dict

    for i in hap_miniport_results_dict:
        results_dict[f"haplotype{i+1}"] = hap_miniport_results_dict[i][0]

    low_coverage = True if results_dict['reseq_stat']['coverage'] < 0.5 else False
    low_depth = True if results_dict['reseq_stat']['depth'] < 5 else False

    fatal_mut = True
    for i in results_dict:
        if i.startswith('haplotype'):
            if results_dict[i]['frameshift'] is False and results_dict[i]['stopcodon'] is False and results_dict[i]['headmissed'] is False and results_dict[i]['coverage'] > 0.9 and results_dict[i]['identity'] > 0.9:
                fatal_mut = False
                break

    ref_bad_mut = results_dict['ref_redo']['frameshift'] or results_dict['ref_redo'][
        'stopcodon'] or results_dict['ref_redo']['headmissed'] or results_dict['ref_redo']['coverage'] < 0.9
    fatal_mut = False if ref_bad_mut else fatal_mut
    fatal_mut = False if low_depth else fatal_mut

    if low_coverage or fatal_mut:
        results_dict['LoF'] = "Yes"
    else:
        if ref_bad_mut or low_depth:
            results_dict['LoF'] = "Not Sure"
        else:
            results_dict['LoF'] = "No"

    result_json = "%s/result.json" % work_dir
    with open(result_json, 'w') as f:
        f.write(json.dumps(results_dict, indent=4))

    # zip work_dir
    shutil.make_archive(work_dir, 'zip', work_dir)
    rmdir(work_dir)
