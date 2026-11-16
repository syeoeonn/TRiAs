"""Merging TR allele tables produced by separate `vcf-alleles` runs."""

from .merge_alleles import merge_alleles_cmd

__all__ = ["merge_alleles_cmd"]
