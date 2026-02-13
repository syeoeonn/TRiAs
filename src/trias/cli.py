# src/trias/cli.py
from __future__ import annotations
import argparse
from .bed.prepare import bed_prepare_cmd
from .allele.extract import vcf_alleles_cmd
from .motif.count import motif_count_cmd
from .summary.motif_summary import motif_summary_cmd
from .dntre.dntre import run_dntre as dntre_cmd, DEFAULT_FMT, default_jobs
from .annotation.overlap_gene import annotate_cmd
from .annotation.moi_annot import annotate_moi_cmd

# from .mc_ratio.mc_ratio import mc_ratio_cmd
# from .allele.report_dntre import report_dntre_cmd
    
def main():
    parser = argparse.ArgumentParser(
        prog="trias",
        description="TRiAs: Tandem Repeat analysis toolkit"
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # bed-prepare
    p = sub.add_parser(
        "bed-prepare",
        help="Produce per-chromosome TR BED files from input BED(+gz) and a reference FASTA"
    )
    p.add_argument("--bed", required=True, help="Input BED file (.bed or .bed.gz)")
    p.add_argument("--fasta", required=True, help="Reference FASTA file")
    p.add_argument("--outdir", required=True, help="Output directory for per-chromosome files")
    p.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    p.add_argument("--chunk", type=int, default=1000, help="Chunk size per thread (default: 1000)")
    p.add_argument("--mapping", default=None, help="Chromosome name mapping file (BED_name<TAB>FASTA_name)")
    p.add_argument("--motif-col", type=int, default=10, help="1-based column index of motif in the input BED (default: 10)")
    p.add_argument("--has-header", action="store_true", help="Set if the first line of the BED is a header to be skipped")
    p.set_defaults(func=bed_prepare_cmd)

    # vcf-alleles
    q = sub.add_parser(
        "vcf-alleles",
        help="Extract population TR alleles per interval using a VCF and group definitions"
    )
    q.add_argument("--vcf", required=True, help="Indexed VCF(.vcf.gz) with .tbi/.csi")
    q.add_argument("--input", required=True, help="Input TSV with columns: CHROM, START, END, MOTIF, TR_allele")
    q.add_argument("--output", required=True, help="Output TSV path")
    q.add_argument(
        "--groups", nargs="+", required=True,
        help=(
            "Group TSV files (one or more). The FIRST file is treated as the REF group. "
            "Format per line: base_sample_name<TAB>ploidy"
        )
    )
    q.add_argument("--processes", type=int, default=32, help="Number of worker processes (default: 32)")
    q.add_argument("--chunk-size", type=int, default=10000, help="Number of intervals per chunk (default: 10000)")
    q.add_argument("--tmpdir", default="./tmp", help="Base temporary directory (default: ./tmp)")
    q.set_defaults(func=vcf_alleles_cmd)

    # motif-count
    r = sub.add_parser(
        "motif-count",
        help="Count per-motif repeats from TR alleles using WIS+gap compression"
    )
    r.add_argument("--input", required=True, help="Input TSV containing TR_allele with motif and grouping columns")
    r.add_argument("--output", required=True, help="Output TSV path")
    r.add_argument("--max-workers", type=int, default=100, help="Max worker processes (default: 100)")
    r.add_argument("--group-chunk-size", type=int, default=8000, help="Number of groups per chunk (default: 8000)")
    r.set_defaults(func=motif_count_cmd)

    # motif-summary
    s = sub.add_parser(
        "motif-summary",
        help="Summarize motif counts per TR across Population and custom groups"
    )
    s.add_argument("--input", required=True, help="Input TSV with TR rows and per-group count columns")
    s.add_argument("--output", required=True, help="Output TSV path")
    s.add_argument(
        "--groups", action="append", default=None,
        help="Path to a text file listing column names for a group; "
             "group label is derived from the file basename. "
             "Use multiple --groups to add multiple groups."
    )
    s.add_argument("--processes", type=int, default=None, help="Number of processes (default: CPU count)")
    s.add_argument("--chunk-size", type=int, default=None, help="Groups per chunk (default: auto)")
    s.set_defaults(func=motif_summary_cmd)

    # report-dntre: Extract dnTRE
    p = sub.add_parser(
        "dntre",
        help="Parallel dnTRE caller (per-chromosome workers)."
    )
    # Glob-enabled inputs
    p.add_argument("--loci", required=True, help="TR loci TSV path or glob (e.g. 'MERGE_hg38/chr*_merge.tsv')")
    p.add_argument("--alleles", required=True, help="Per-allele TSV path or glob (e.g. 'TR_population_chr*_mc.tsv')")
    # Multi-group (comma-separated)
    p.add_argument("--case", required=True, help="Case column name(s), comma-separated (allele files)")
    p.add_argument("--control", required=True, help="Control column name(s), comma-separated (allele files)")

    # Loci prefixes (default = first of each list)
    p.add_argument("--loci-case-prefix", default=None, help="Case prefix used in LOCI columns (defaults to first of --case)")
    p.add_argument("--loci-control-prefix", default=None, help="Control prefix used in LOCI columns (defaults to first of --control)")

    # Thresholds / formats
    p.add_argument("--threshold", type=float, default=1.0, help="MC ratio threshold (default: 1.0)")
    p.add_argument("--min-pop-mean", type=float, default=1.0, help="Minimum population mean for ratio calc (default: 1.0)")
    p.add_argument("--fmt", default=DEFAULT_FMT, help=f"Numeric format for strings (default: {DEFAULT_FMT})")
    # Trio and hap labels
    p.add_argument("--trio", default=None, help="Optional trio TSV: child<TAB>father<TAB>mother (header ok)")
    p.add_argument("--hap1-label", default="haplotype1", help="Label for .1 haplotype (default: haplotype1)")
    p.add_argument("--hap2-label", default="haplotype2", help="Label for .2 haplotype (default: haplotype2)")
    # Outputs
    p.add_argument("--out-dntre", default="dnTRE.tsv", help="Output TSV for dnTRE alleles (default: dnTRE.tsv)")
    p.add_argument("--out-loci", default="dnTRE_loci.tsv", help="Output TSV for loci that have dnTRE (default: dnTRE_loci.tsv)")
    p.add_argument("--out-alleles", default="all_alleles_in_dnTRE_loci.tsv", help="Output TSV for all alleles in dnTRE loci (default: all_alleles_in_dnTRE_loci.tsv)")
    # Parallel / temp
    p.add_argument("--jobs", type=int, default=default_jobs(), help="Number of parallel workers (default: CPU count)")
    p.add_argument("--tmpdir", default="./_dntre_tmp", help="Temporary directory to write per-chr files (default: ./_dntre_tmp)")
    p.add_argument("--keep-temp", action="store_true", help="Keep temporary per-chr files (default: remove after merge)")
    p.set_defaults(func=dntre_cmd)
    
    # Annotation with Gene & Region
    p = sub.add_parser(
        "annotate",
        help="Annotate regions by adding overlapping gene information (keeps original columns)."
    )
    p.add_argument("--annotation", required=True, help="Annotation file with columns: chrom, start, end, region_type, strand, gene. Supports mixed space/tab delimiters.")
    p.add_argument("--regions", required=True, help="Regions TSV. If a header exists, must include Chrom/Start/End. " "Otherwise, the first three columns are treated as Chrom/Start/End.")
    p.add_argument("--output", required=True, help="Output TSV file path.")
    p.set_defaults(func=annotate_cmd)
    
    # Annotation with MOI (Mode Of Inheritance) & disease name
    p = sub.add_parser(
        "annotate-moi",
        help="Append MOI and disease annotations by mapping genes (from 'Ann_Genes' or a specified column) to a MOI/disease table."
    )
    p.add_argument("--regions", required=True, help="Input TSV (e.g., the output of 'trias annotate' with 'Ann_Genes' column).")
    p.add_argument("--mapping", required=True, help="Gene-to-MOI/disease mapping TSV. Header must include a gene column "
            "(gene|gene_symbol|symbol). MOI candidates: (moi|moi_title|mode_of_inheritance). "
            "Disease candidates: (disease|diseases|phenotype|disorder|condition|omim_disease).")
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--gene-col-candidates", default="Ann_Genes,Gene,Genes", help="Comma-separated candidate names to locate the gene column in --regions. ""Default: 'Ann_Genes,Gene,Genes'.")
    p.add_argument("--out-moi-col", default="Ann_MOI", help="Name of the output MOI column. Default: Ann_MOI")
    p.add_argument("--out-dis-col", default="Ann_Diseases", help="Name of the output disease column. Default: Ann_Diseases")
    p.set_defaults(func=annotate_moi_cmd)
    
        # # mc-ratio
    # sp = sub.add_parser(
    # "mc-ratio",
    # help="Compute MC ratios with optional dual filtering for case/control groups."
    # )
    # sp.add_argument("--input", "-i", required=True, help="Input TSV file")
    # sp.add_argument("--output", "-o", default="-", help="Output TSV file (default: stdout)")
    # sp.add_argument("--case", required=True, help="Case group name (e.g., Our)")
    # sp.add_argument("--control", required=True, help="Control group name (e.g., Unaffected)")
    # sp.add_argument("--threshold", type=float, default=2.0, help="MC ratio threshold (default: 2.0)")
    # sp.add_argument("--min-pop-mean", type=float, default=1.0,
    #                 help="Minimum population mean; values below are treated as NaN (default: 1.0)")
    # sp.add_argument("--diff", action="store_true", help="Add Case-Control max MC difference columns")
    # sp.add_argument("--dual-filter", action="store_true",
    #                 help="Enable dual filtering: Case MC > threshold AND Control MC ≤ threshold")
    # sp.add_argument("--max", action="store_true",
    #     help="Enable additional check: Case max > Control max at the same motif position")
    # sp.add_argument("--fmt", default=".2f", help="Numeric format (default: .2f)")
    # sp.set_defaults(func=mc_ratio_cmd)
    
    # # report-dntre
    # t = sub.add_parser(
    #     "report-dntre",
    #     help="Report dnTRE alleles and ratios (trio mode or auto-sample mode)"
    # )
    # t.add_argument("--mc", required=True, help="mc-ratio TSV (flags, group ratios, population means)")
    # t.add_argument("--population", required=True,
    #                help="Glob pattern for TR population TSVs, e.g. 'TR_population_chr*_mc.tsv'")
    # t.add_argument("--output", "-o", default="-", help="Output TSV (default: stdout)")
    # t.add_argument("--trio", help="Optional TSV: child<TAB>father<TAB>mother (no/header)")
    # t.add_argument("--threshold", type=float, default=2.0, help="MC ratio threshold (>=)")
    # t.add_argument("--min-pop-mean", type=float, default=1.0, help="Min population mean for ratios")
    # t.add_argument("--fmt", default=".2f", help="Numeric format for ratio strings (default: .2f)")
    # t.add_argument("--hap1-label", default="haplotype1", help="Label for .1 haplotype (e.g., 'maternal')")
    # t.add_argument("--hap2-label", default="haplotype2", help="Label for .2 haplotype (e.g., 'paternal')")
    # t.set_defaults(func=report_dntre_cmd)
    
    args = parser.parse_args()
    return args.func(args)
