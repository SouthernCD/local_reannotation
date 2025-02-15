import argparse

class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse
        parser = argparse.ArgumentParser(
            prog='reanno',
            description="Locally reannotate gene structures from haplotype sequences"
        )

        parser.add_argument('gene_id', type=str, help='Gene ID')
        parser.add_argument('mRNA_pkl', type=str, help='mRNA pickle file')
        parser.add_argument('work_dir', type=str, help='Working directory')
        parser.add_argument('bam_file', type=str, help='Bam file')
        parser.add_argument('haplo_seq_fa', type=str, help='Haplotype sequence fasta file')
        parser.add_argument('ref_sequence_fa', type=str, help='Reference sequence fasta file')

        self.arg_parser = parser
        self.args = parser.parse_args()

    def run(self):
        self.run_arg_parser()
        from local_reannotation.src.pipeline import gen_annotation
        gen_annotation(self.args.gene_id, self.args.mRNA_pkl, self.args.work_dir, self.args.bam_file, self.args.haplo_seq_fa, self.args.ref_sequence_fa)

def main():
    job = Job()
    job.run()

if __name__ == '__main__':
    main()