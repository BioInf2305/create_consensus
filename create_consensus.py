import sys
import pysam
import argparse


class ReferenceAlleleMismatchError(Exception):
    """Raised when the reference allele does not match the VCF allele."""

    pass

class GenotypeFormatError(Exception):
    """Raised when the genotypes in the input vcf files are in Unexpected format"""

    pass


def extract_ref_seq(ref_fasta, chrom, start, end):
    """extract the reference sequence given the coordinates and chromosome"""

    fa = pysam.FastaFile(ref_fasta)
    seq = fa.fetch(chrom, start - 1, end)
    return seq


def vcf_to_consensus(seq, vcf, chrom, start, end, sample_list, plink_vcf=False):
    """create dictionary of the consensus sequences using vcf and ref seq:{"sample1":{10:"A",20:"T},"sample2":{10:"A",20:"G"}}"""

    seq_dict = {}
    sample_seq_dict = {}
    for i, v in enumerate(list(range(start - 1, end))):
        seq_dict[v + 1] = seq[i].capitalize()
    vcf_in = pysam.VariantFile(vcf)
    if len(sample_list) == 0:
        sample_list = list(vcf_in.header.samples)
    # Dictionary mapping heterozygous genotypes to IUPAC codes (without "/")
    heterozygous_to_iupac = {
        "AG": "R",  # purine
        "GA": "R",
        "CT": "Y",  # pyrimidine
        "TC": "Y",
        "GC": "S",
        "CG": "S",
        "AT": "W",
        "TA": "W",
        "GT": "K",
        "TG": "K",
        "AC": "M",
        "CA": "M",
    }

    for sample in sample_list:
        sample_seq_dict[sample] = seq_dict.copy()
    for rec in vcf_in.fetch(chrom, start, end):
        if plink_vcf:
            if (
                seq_dict[rec.pos].capitalize() != rec.ref.capitalize()
                and seq_dict[rec.pos].capitalize() != rec.alts[0].capitalize()
            ):
                raise ReferenceAlleleMismatchError(
                    f"reference base does not match for {rec.pos}, ref: {seq_dict[rec.pos]}, vcf:{rec.ref}"
                )
        else:
            if seq_dict[rec.pos].capitalize() != rec.ref.capitalize():
                raise ReferenceAlleleMismatchError(
                    f"reference base does not match for {rec.pos}, ref: {seq_dict[rec.pos]}, vcf:{rec.ref}"
                )
        for sample in sample_list:
            gt = rec.samples[sample]["GT"]
            if gt == (None, None):
                sample_seq_dict[sample][rec.pos] = "N"
            elif gt == (0, 0):
                sample_seq_dict[sample][rec.pos] = rec.ref.capitalize()
            elif gt == (1, 1):
                sample_seq_dict[sample][rec.pos] = rec.alts[0].capitalize()
            elif gt == (1, 0) or gt == (0, 1):
                sample_seq_dict[sample][rec.pos] = heterozygous_to_iupac[
                    rec.ref + rec.alts[0]
                ].capitalize()
            else:
                raise GenotypeFormatError(f"genotypes in unexpected format, {gt} for {sample} at {rec.pos}")
    return sample_seq_dict


def seq_dict_to_fasta(sample_seq_dict, output):
    """write dictionary to fasta file"""

    with open(f"{output}", "w") as dest:
        for sample in sample_seq_dict:
            dest.write(f">{sample}\n")
            seq = "".join(list(sample_seq_dict[sample].values()))
            wrapped = "\n".join([seq[i : i + 60] for i in range(0, len(seq), 60)])
            dest.write(f"{wrapped}\n")


def main():
    parser = argparse.ArgumentParser(
        description="create consensus fasta files given vcf and reference fasta file"
    )
    parser.add_argument("-v", "--vcf", required=True, help="input vcf file")
    parser.add_argument("-f", "--ref_fasta", required=True, help="reference fasta file")
    parser.add_argument("-c", "--chrom", required=True, help="chromosome name")
    parser.add_argument("-s", "--start", required=True, help="start coordinates")
    parser.add_argument("-e", "--end", required=True, help="end coordinates")
    parser.add_argument("-S", "--sample_list", required=False, default=None, help="chromosome name")
    parser.add_argument(
        "-o", "--output", required=False, default="output.fa", help="output_file"
    )
    parser.add_argument(
        "--plink_vcf", action="store_true", help="vcf file is generated using plink."
    )

    args = parser.parse_args()

    sample_list = []
    if args.sample_list:
        with open(args.sample_list) as source:
            for line in source:
                line = line.rstrip().split()
                sample_list.append(line[0])


    seq = extract_ref_seq(args.ref_fasta, args.chrom, int(args.start), int(args.end))
    sample_seq_dict = vcf_to_consensus(
        seq,
        args.vcf,
        args.chrom,
        int(args.start),
        int(args.end),
        sample_list,
        plink_vcf=args.plink_vcf,
    )
    if args.output == "output.fa":
        args.output = f"{args.chrom}_{args.start}_{args.end}.fa"

    seq_dict_to_fasta(sample_seq_dict, output=args.output)


if __name__ == "__main__":
    main()
