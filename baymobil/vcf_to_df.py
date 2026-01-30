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


def parse_vcf_ad(vcf_path, foreign=False):
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
            fmt = fields[8]
            sample = fields[9]

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