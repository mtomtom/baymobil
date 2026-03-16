import gzip
import pandas as pd

def _open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def _parse_ad(ad_val):
    """
    Parse FORMAT/AD and return (n_ref, n_alt).

    AD formats:
      ref,alt
      ref
      .
    """
    if not ad_val or ad_val == ".":
        return None, None

    parts = ad_val.split(",")

    try:
        if len(parts) >= 2:
            return int(parts[0]), int(parts[1])
        elif len(parts) == 1:
            return int(parts[0]), 0
    except ValueError:
        return None, None

    return None, None


def parse_vcf_ad(vcf_path, foreign=False, snp_only=True):
    """
    Parse a single-sample VCF and return a DataFrame with:
      SNP, N, n, n_ref, n_alt, GeneID, Effect

    Definitions:
      n_ref = AD[0]
      n_alt = AD[1]
      N     = n_ref + n_alt

    foreign:
      False → n = n_alt
      True  → n = n_ref
    snp_only:
      True → keep only single-base REF/ALT (exclude indels and multi-allelics)
             but allow ALT="." for reference-only calls at those sites
    """
    records = []

    with _open_maybe_gzip(vcf_path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom, pos = fields[0], fields[1]
            ref, alt = fields[3], fields[4]
            fmt = fields[8]
            sample = fields[9]

            if snp_only:
                if len(ref) != 1:
                    continue
                if alt != "." and len(alt) != 1:
                    continue
                if "," in alt:
                    continue

            snp = f"{chrom}_{pos}"

            # FORMAT parsing
            fmt_keys = fmt.split(":")
            sample_vals = sample.split(":")
            fmt_map = dict(zip(fmt_keys, sample_vals))

            # AD parsing (authoritative)
            n_ref, n_alt = _parse_ad(fmt_map.get("AD", "."))

            if n_ref is None:
                N = None
                n = None
            else:
                N = n_ref + n_alt
                n = n_ref if foreign else n_alt


            records.append({
                "SNP": snp,
                "CHROM": chrom,
                "POS": pos,
                "N": N,
                "n": n,
                "n_ref": n_ref,
                "n_alt": n_alt
            })

    return pd.DataFrame(records)


def extract_all_snp_gene_table(vcf_path):
    rows = []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt, qual, flt, info = fields[:8]

            gene = None
            for entry in info.split(";"):
                if entry.startswith("GENE="):
                    gene = entry.split("=", 1)[1]
                    break

            snp_id = f"{chrom}_{pos}"
            rows.append({
                "SNP": snp_id,
                "CHROM": chrom,
                "POS": int(pos),
                "REF": ref,
                "ALT": alt,
                "GENE": gene,
            })
    return pd.DataFrame(rows)
