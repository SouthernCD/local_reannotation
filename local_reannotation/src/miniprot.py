from yxutil import cmd_run, mkdir
from yxseq import read_gff_file, Gene, ChrLoci
from local_reannotation.src.utils import read_sequence_from_file


def miniprot_job(query_file, subject_file, work_dir):
    mkdir(work_dir, False)
    miniprot_out = "%s/miniprot.out" % work_dir
    run_miniprot(query_file, subject_file, miniprot_out, work_dir)
    miniprot_out_gff, miniprot_out_aln = split_miniprot_out(miniprot_out)
    hit_gf = parse_miniprot_gff(miniprot_out_gff, subject_file)

    if hit_gf:
        pseudo_flag = is_pseudo(hit_gf, query_file, work_dir)
        return pseudo_flag, hit_gf.model_aa_seq, hit_gf.model_cds_seq
    else:
        pseudo_flag = {
            'frameshift': False,
            'stopcodon': False,
            'headmissed': False,
            'identity': 0.0,
            'coverage': 0.0
        }
        target_pt_seq = ""
        target_cds_seq = ""
        return pseudo_flag, target_pt_seq, target_cds_seq


def run_miniprot(query_file, subject_file, miniprot_out, work_dir=None):
    cmd_string = "miniprot --gff --aln %s %s > %s" % (
        subject_file, query_file, miniprot_out)
    # print(cmd_string)
    f,o,e = cmd_run(cmd_string, cwd=work_dir)
    # print(f,o,e)
    return miniprot_out


def split_miniprot_out(miniprot_out):
    miniprot_out_gff = miniprot_out + '.gff'
    miniprot_out_aln = miniprot_out + '.aln'

    gff_fh = open(miniprot_out_gff, 'w')
    aln_fh = open(miniprot_out_aln, 'w')

    with open(miniprot_out, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('##gff-version') or not line.startswith('#'):
                gff_fh.write(line + '\n')
            else:
                aln_fh.write(line + '\n')

    gff_fh.close()
    aln_fh.close()

    return miniprot_out_gff, miniprot_out_aln


def parse_miniprot_gff(miniprot_out_gff, miniprot_subject_file):
    """
    Warning: only work for one query miniprot output
    """

    gff_dict = read_gff_file(miniprot_out_gff)
    if len(gff_dict) == 0:
        return None

    gff_dict = gff_dict['mRNA']

    best_hit_id = None
    best_score = 0
    for hit_id in gff_dict:
        gf = gff_dict[hit_id]
        score = float(gf.qualifiers['score'][0])
        if score > best_score:
            best_score = score
            best_hit_id = hit_id

    best_gf = gff_dict[best_hit_id]

    gene_loci = ChrLoci(best_gf.chr_id, best_gf.strand,
                        best_gf.start, best_gf.end)
    gene = Gene(id=best_gf.id+".gene", chr_loci=gene_loci,
                sub_features=[best_gf])

    gene.build_gene_seq(miniprot_subject_file)
    gene.qualifiers = best_gf.qualifiers

    return gene


def is_pseudo(gene, query_file, work_dir):
    frameshift_flag = False
    stopcodon_flag = False
    headmissed_flag = False

    if 'Frameshift' in gene.qualifiers:
        frameshift_flag = True
    if 'StopCodon' in gene.qualifiers:
        stopcodon_flag = True

    start = int(gene.qualifiers['Target'][0].split()[1])
    end = int(gene.qualifiers['Target'][0].split()[2])
    if start != 1:
        if gene.model_aa_seq is not None:
            if gene.model_aa_seq.startswith('M'):
                headmissed_flag = False
            else:
                headmissed_flag = True
        else:
            headmissed_flag = True

    # identity and coverage
    anno_pt_fasta = work_dir + '/anno_pt.fasta'
    with open(anno_pt_fasta, 'w') as f:
        f.write(">anno_pt\n%s\n" % gene.model_aa_seq)

    cmd_run("blastp -query %s -subject %s -outfmt 6 -out %s.blast" %
            (anno_pt_fasta, query_file, work_dir + '/anno_pt'), cwd=work_dir)

    identity = 0.0
    coverage = 0.0
    query_len = len(read_sequence_from_file(query_file))
    with open(work_dir + '/anno_pt.blast', 'r') as f:
        for line in f:
            line = line.strip()
            identity = (
                float(line.split('\t')[2]) * float(line.split('\t')[3]))/query_len/100
            coverage = (
                abs(int(line.split('\t')[9]) - int(line.split('\t')[8])) + 1) / query_len
            break

    return {
        'frameshift': frameshift_flag,
        'stopcodon': stopcodon_flag,
        'headmissed': headmissed_flag,
        'identity': identity,
        'coverage': coverage
    }


def parse_miniprot_aln(miniprot_out_aln):
    with open(miniprot_out_aln, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('##ATA'):
                target_pt_str = line.split('\t')[1]
            elif line.startswith('##AAS'):
                similarity_str = line.split('\t')[1]
            elif line.startswith('##AQA'):
                query_pt_str = line.split('\t')[1]
            elif line.startswith('##ATN'):
                target_nt_str = line.split('\t')[1]
            else:
                pass

    target_pt_seq = target_pt_str.replace('.', '').replace(
        ' ', '').replace('!!', 'X').replace('!', 'X')
    target_pt_seq = target_pt_str.replace('.', '').replace(
        ' ', '').replace('$$', 'X').replace('$', 'X')
    endwith_star = False
    if target_pt_seq.endswith('*'):
        endwith_star = True
        target_pt_seq = target_pt_seq[:-1]
    target_pt_seq = target_pt_seq.replace('*', 'X')
    if endwith_star:
        target_pt_seq += '*'

    # pseudo_flag = False
    # if 'X' in target_pt_seq or '*' in target_pt_seq:
    #     pseudo_flag = True

    target_cds_seq = ''.join(
        [char for char in target_nt_str if not char.islower()])

    return target_pt_seq, target_cds_seq


if __name__ == '__main__':
    query_file = 'ref_prot.fasta'
    subject_file = 'hap1.pseudo.fasta'
    miniprot_out = 'hap1.pseudo.miniprot.out'
    run_miniprot(query_file, subject_file, miniprot_out)
    miniprot_out_gff, miniprot_out_aln = split_miniprot_out(miniprot_out)
    hit_gf = parse_miniprot_gff(miniprot_out_gff)
    target_pt_seq, target_cds_seq = parse_miniprot_aln(miniprot_out_aln)
    is_pseudo(hit_gf, target_pt_seq)
