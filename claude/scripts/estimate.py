#!/usr/bin/env python3
"""Resource model for kmindex build / merge / query.

Stand-alone (no kmhelpers dependency). Formulas and their sources are listed
in claude/performances_review.md. All sizes are in bytes unless stated.

Examples:
  # one chunk of 5000 samples, biggest sample 5e8 k-mers, 64 threads
  ./estimate.py --kmers 5e8 --samples 5000 --total-samples 100000 \
      --threads 64 --ulimit 1048576 --ram 237e9 --fpr 0.25

  # sweep the number of chunks for a fixed total
  ./estimate.py --kmers 5e8 --total-samples 100000 --threads 64 \
      --ulimit 1048576 --ram 237e9 --sweep 1,2,4,8,16,32,64 --csv
"""

import argparse
import csv
import math
import sys

# constants read in the code (see the review report for file:line)
BYTE_PER_KMER = 8            # kmtricks count: uint64 hash per k-mer
COUNT_EXTRA = 8192           # get_required_memory_hash
HASH_READER_BUF = 5 * 32768  # HashReader<32768>: IFile buf + 4 arrays
HASH_RECORD = 12             # uint64 hash + uint32 count (default MAX_C) in hash files
SUPERK_BUF = 32768           # SuperKStorageWriter buffer per partition
HEADER = 49                  # .cmbf header
WINDOW_ROUND = 64            # bloom window rounded to 64 bits per partition
SMER_ENTRY = 24              # (smer{i,p,h}, qid) in batch_query
MIN_PARTITIONS = 4


def f_value(fpr):
    return -math.log(fpr) / (math.log(2) ** 2)


def bf_size_from_fpr(kmers, fpr):
    # kmhelpers SpanManager.get_bf_size, rounded to a byte
    return ((int(f_value(fpr) * kmers) + 7) // 8) * 8


def effective_bf(bf_size, partitions):
    window = math.ceil(bf_size / partitions)
    window = ((window + WINDOW_ROUND - 1) // WINDOW_ROUND) * WINDOW_ROUND
    return window * partitions, window


def storage(bf_size, partitions, samples, bw=1):
    bf_eff, window = effective_bf(bf_size, partitions)
    row = (samples * bw + 7) // 8
    return partitions * (HEADER + window * row)


def nb_partitions_ram(kmers, ram, threads, margin=1.05):
    # kmparams.nb_partitions: RAM floor of the count stage
    per_thread = ram / threads
    kmers_per_part = per_thread / (BYTE_PER_KMER * margin)
    return max(MIN_PARTITIONS, math.ceil(kmers / kmers_per_part))


def open_files(threads, partitions, samples, focus=0.5):
    nw = max(1, math.floor(threads * focus))
    return {
        "superk_exact": nw * partitions + 2 * (threads - nw),
        "superk_kmhelpers": threads * partitions + nw,
        "count": 2 * threads,
        "merge_exact": min(threads, partitions) * (samples + 1),
        "merge_kmhelpers": threads * (samples + 1),
    }


def ram_build(kmers, partitions, threads, samples, focus=0.5, gamma=1.2):
    nw = max(1, math.floor(threads * focus))
    per_count = kmers * gamma / partitions * BYTE_PER_KMER + COUNT_EXTRA
    return {
        "count_stage": threads * per_count,
        "count_stage_kmhelpers": kmers * BYTE_PER_KMER * 1.05 * threads / partitions,
        "superk_stage": nw * partitions * SUPERK_BUF,
        "merge_stage": min(threads, partitions) * samples * HASH_READER_BUF,
    }


def temp_disk(mean_kmers, samples, partitions, cpr_ratio=1.0, superk_frac=1.0,
              k=25, m=10):
    # peak at the end of the count stage: all hash files, plus the fraction of
    # super-k-mer files not yet consumed (superk_frac = 1 is the safe bound)
    ell = (k - m + 2) / 2                      # k-mers per super-k-mer
    superk_per_kmer = 0.25 + ((k - 1) / 4 + 1) / ell
    total = mean_kmers * samples
    hash_bytes = HASH_RECORD * total * cpr_ratio
    superk_bytes = superk_per_kmer * total * superk_frac
    return {
        "hash_files_bytes": hash_bytes,
        "superk_files_bytes": superk_bytes,
        "peak_bytes": hash_bytes + superk_bytes,
        "temp_files": 2 * samples * partitions,
    }


def ram_query(samples, seq_len, smer_size, n_seq, threads, batch, bw=1):
    row = (samples * bw + 7) // 8
    smers = max(0, seq_len - smer_size + 1)
    per_seq = smers * (row + SMER_ENTRY)
    live = n_seq if batch == 0 else min(n_seq, threads * batch)
    return {"bytes_per_smer_row": row, "per_sequence": per_seq, "peak": live * per_seq}


def batch_size(ram_budget, samples, max_seq_len, smer_size, threads, bw=1):
    # largest -b such that T * b sequences of the longest length fit in ram_budget
    per_seq = ram_query(samples, max_seq_len, smer_size, 1, 1, 1, bw)["per_sequence"]
    return int(ram_budget // (threads * per_seq))


def time_build(bases, kmers_total, index_bytes, hash_bytes, threads, partitions,
               parse_bps, count_kps, rd_bps, wr_bps, bitop_ps):
    # crude stage model, coefficients are inputs (fit them with E3)
    superk = bases / (parse_bps * threads)
    count = kmers_total / (count_kps * threads) + hash_bytes / wr_bps
    par = min(threads, partitions)
    merge_io = hash_bytes / rd_bps + index_bytes / wr_bps
    merge_cpu = index_bytes * 8 / (bitop_ps * par)
    return {"superk": superk, "count": count, "kmtricks_merge": merge_io + merge_cpu,
            "total": superk + count + merge_io + merge_cpu}


def time_kmindex_merge(index_bytes, first_chunk_bytes, threads, partitions,
                       rd_bps, wr_bps, bitop_ps, levels=1):
    par = min(threads, partitions)
    io = levels * (index_bytes / rd_bps + index_bytes / wr_bps)
    cpu = levels * (index_bytes - first_chunk_bytes) * 8 / (bitop_ps * par)
    return {"io": io, "cpu": cpu, "total": io + cpu}


def fmt_bytes(b):
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(b) < 1000:
            return f"{b:.3g} {unit}"
        b /= 1000
    return f"{b:.3g} EB"


def fmt_time(s):
    if s < 120:
        return f"{s:.0f} s"
    if s < 7200:
        return f"{s / 60:.1f} min"
    if s < 172800:
        return f"{s / 3600:.1f} h"
    return f"{s / 86400:.1f} d"


def one_config(a, samples, chunks, out):
    bf = a.bf_size or bf_size_from_fpr(a.kmers, a.fpr)
    P = a.partitions or nb_partitions_ram(a.kmers, a.ram, a.threads)
    bf_eff, window = effective_bf(bf, P)
    files = open_files(a.threads, P, samples, a.focus)
    ram = ram_build(a.kmers, P, a.threads, samples, a.focus, a.gamma)
    tmp = temp_disk(a.mean_kmers, samples, P, a.cpr_ratio, a.superk_frac, a.k, a.m)
    size_chunk = storage(bf, P, samples)
    size_total = storage(bf, P, a.total_samples)
    tb = time_build(a.bases_per_sample * samples, a.mean_kmers * samples, size_chunk,
                    tmp["hash_files_bytes"], a.threads, P, a.parse_bps, a.count_kps,
                    a.read_bps, a.write_bps, a.bitop_ps)
    tm = time_kmindex_merge(size_total, size_chunk, a.merge_threads or a.threads, P,
                            a.read_bps, a.write_bps, a.bitop_ps, a.merge_levels)
    T_eff = min(a.merge_threads or a.threads, P)
    q = ram_query(a.total_samples, a.query_len, a.smer_size, a.query_nseq,
                  a.query_threads, a.query_batch)
    nodes = a.nodes or chunks
    waves = math.ceil(chunks / nodes)
    row = {
        "chunks": chunks, "samples_per_chunk": samples, "partitions": P,
        "bf_bits": bf, "bf_bits_effective": bf_eff, "bits_per_kmer": bf / a.kmers,
        "storage_chunk": size_chunk, "storage_total": size_total,
        "files_superk": files["superk_exact"], "files_merge": files["merge_exact"],
        "files_merge_bound": files["merge_kmhelpers"],
        "files_kmindex_merge": T_eff * (chunks + 1),
        "ram_count": ram["count_stage"], "ram_kmtricks_merge": ram["merge_stage"],
        "ram_superk": ram["superk_stage"],
        "temp_hash_bytes": tmp["hash_files_bytes"], "temp_peak_bytes": tmp["peak_bytes"],
        "temp_files": tmp["temp_files"],
        "t_build_chunk": tb["total"], "t_build_waves": waves * tb["total"],
        "t_kmindex_merge": tm["total"], "t_wall": waves * tb["total"] + tm["total"],
        "node_hours": (chunks * tb["total"] + tm["total"]) / 3600,
        "query_peak_ram": q["peak"], "query_row_bytes": q["bytes_per_smer_row"],
        "ok_files": files["merge_exact"] <= a.ulimit,
        "ok_ram": max(ram.values()) <= a.ram,
        "ok_kmindex_merge_files": T_eff * (chunks + 1) <= a.ulimit,
        "ok_walltime": tb["total"] <= a.walltime,
        "ok_tmp": tmp["peak_bytes"] <= a.scratch,
    }
    if out is not None:
        out.append(row)
    return row, files, ram, tmp, tb, tm, q


def print_report(a, row, files, ram, tmp, tb, tm, q):
    p = print
    p("# parameters")
    p(f"K (max k-mers/sample) = {a.kmers:.3g}, mean k-mers/sample = {a.mean_kmers:.3g}")
    p(f"S (samples/chunk) = {row['samples_per_chunk']}, S_total = {a.total_samples}, "
      f"chunks = {row['chunks']}")
    p(f"P = {row['partitions']}, T = {a.threads}, focus = {a.focus}, "
      f"ulimit = {a.ulimit}, RAM = {fmt_bytes(a.ram)}")
    p(f"bf_size = {row['bf_bits']:.4g} bits ({row['bits_per_kmer']:.3f} bits/k-mer), "
      f"effective after 64-bit window rounding = {row['bf_bits_effective']:.4g}")
    p()
    p("# storage (exact: P*49 + window*ceil(S/8) per partition)")
    p(f"chunk index = {fmt_bytes(row['storage_chunk'])}, merged index = "
      f"{fmt_bytes(row['storage_total'])}")
    p(f"temp hash files per chunk (12 B/k-mer, cpr ratio {a.cpr_ratio}) = "
      f"{fmt_bytes(tmp['hash_files_bytes'])}, super-k-mer files (fraction {a.superk_frac}) = "
      f"{fmt_bytes(tmp['superk_files_bytes'])}, peak = {fmt_bytes(tmp['peak_bytes'])}, "
      f"temp files = {tmp['temp_files']}")
    p()
    p("# open files (per process)")
    for k, v in files.items():
        p(f"{k:18s} {v:>12d} {'OK' if v <= a.ulimit else 'EXCEEDS ulimit'}")
    p(f"{'kmindex_merge':18s} {row['files_kmindex_merge']:>12d} "
      f"(min(T,P) * (chunks+1))")
    p()
    p("# peak RAM by stage")
    for k, v in ram.items():
        p(f"{k:22s} {fmt_bytes(v):>10s} {'OK' if v <= a.ram else 'EXCEEDS RAM'}")
    p(f"{'kmindex_merge (anon)':22s} {fmt_bytes(row['partitions'] * math.ceil(a.total_samples / 8)):>10s}")
    p(f"{'query per sequence':22s} {fmt_bytes(q['per_sequence']):>10s} "
      f"({q['bytes_per_smer_row']} B per s-mer row, L={a.query_len})")
    p(f"{'query peak':22s} {fmt_bytes(q['peak']):>10s} "
      f"(nseq={a.query_nseq}, T={a.query_threads}, b={a.query_batch})")
    b_max = batch_size(a.query_ram or a.ram, a.total_samples, a.query_len, a.smer_size,
                       a.query_threads)
    p(f"{'query max -b':22s} {b_max:>10d} "
      f"(budget {fmt_bytes(a.query_ram or a.ram)}, T={a.query_threads}, L={a.query_len}; "
      f"0 means reduce T)")
    p()
    p("# time (coefficients are inputs, fit them with experiment E3)")
    for k, v in tb.items():
        p(f"build {k:16s} {fmt_time(v):>10s}")
    p(f"build per chunk = {fmt_time(tb['total'])}, walltime limit "
      f"{'OK' if row['ok_walltime'] else 'EXCEEDED'}")
    p(f"kmindex merge: io {fmt_time(tm['io'])}, cpu {fmt_time(tm['cpu'])}, "
      f"total {fmt_time(tm['total'])} (levels={a.merge_levels})")
    p(f"wall (chunks on {a.nodes or row['chunks']} nodes) = {fmt_time(row['t_wall'])}, "
      f"node-hours = {row['node_hours']:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kmers", type=float, required=True, help="max k-mers in one sample")
    ap.add_argument("--mean-kmers", type=float, help="mean k-mers per sample (default: kmers)")
    ap.add_argument("--bases-per-sample", type=float, help="input bases per sample (default: mean k-mers)")
    ap.add_argument("--samples", type=int, help="samples per chunk (default: total/chunks)")
    ap.add_argument("--total-samples", type=int, required=True)
    ap.add_argument("--chunks", type=int, default=0, help="number of chunks (default: from samples)")
    ap.add_argument("--partitions", type=int, default=0, help="0 = RAM-driven kmparams value")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--merge-threads", type=int, default=0)
    ap.add_argument("--merge-levels", type=int, default=1)
    ap.add_argument("--focus", type=float, default=0.5)
    ap.add_argument("--gamma", type=float, default=1.2, help="partition imbalance factor")
    ap.add_argument("--bf-size", type=float, default=0, help="bits, overrides --fpr")
    ap.add_argument("--fpr", type=float, default=0.25)
    ap.add_argument("--cpr-ratio", type=float, default=1.0, help="temp size ratio with --cpr")
    ap.add_argument("--superk-frac", type=float, default=1.0,
                    help="fraction of super-k-mer files present at the temp-disk peak (1 = bound)")
    ap.add_argument("-k", type=int, default=25)
    ap.add_argument("-m", type=int, default=10, help="minimizer size")
    ap.add_argument("--ulimit", type=int, default=1048576)
    ap.add_argument("--ram", type=float, default=237e9)
    ap.add_argument("--scratch", type=float, default=286e9, help="temp disk available")
    ap.add_argument("--walltime", type=float, default=48 * 3600)
    ap.add_argument("--nodes", type=int, default=0, help="nodes building chunks in parallel (default: chunks)")
    ap.add_argument("--read-bps", type=float, default=1e9)
    ap.add_argument("--write-bps", type=float, default=1e9)
    ap.add_argument("--parse-bps", type=float, default=2e7, help="bases/s/thread in superk")
    ap.add_argument("--count-kps", type=float, default=1e7, help="k-mers/s/thread in count")
    ap.add_argument("--bitop-ps", type=float, default=5e8, help="bit ops/s/thread in merges")
    ap.add_argument("--query-len", type=int, default=1000)
    ap.add_argument("--query-nseq", type=int, default=100)
    ap.add_argument("--query-threads", type=int, default=8)
    ap.add_argument("--query-batch", type=int, default=0)
    ap.add_argument("--query-ram", type=float, default=0,
                    help="RAM budget for the -b recommendation (default: --ram)")
    ap.add_argument("--smer-size", type=int, default=25)
    ap.add_argument("--sweep", help="comma-separated chunk counts")
    ap.add_argument("--csv", action="store_true", help="CSV output for --sweep")
    a = ap.parse_args()
    a.kmers = int(a.kmers)
    a.mean_kmers = int(a.mean_kmers or a.kmers)
    a.bases_per_sample = a.bases_per_sample or a.mean_kmers
    a.bf_size = int(a.bf_size)

    if a.sweep:
        rows = []
        for c in (int(x) for x in a.sweep.split(",")):
            one_config(a, math.ceil(a.total_samples / c), c, rows)
        if a.csv:
            w = csv.DictWriter(sys.stdout, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        else:
            keys = ["chunks", "samples_per_chunk", "files_merge", "ram_kmtricks_merge",
                    "t_build_chunk", "t_kmindex_merge", "t_wall", "node_hours",
                    "ok_files", "ok_ram", "ok_walltime", "ok_tmp"]
            print("\t".join(keys))
            for r in rows:
                vals = []
                for k in keys:
                    v = r[k]
                    if k.startswith("ram") or k.startswith("storage"):
                        v = fmt_bytes(v)
                    elif k.startswith("t_"):
                        v = fmt_time(v)
                    elif isinstance(v, float):
                        v = f"{v:.1f}"
                    vals.append(str(v))
                print("\t".join(vals))
        return

    chunks = a.chunks or (math.ceil(a.total_samples / a.samples) if a.samples else 1)
    samples = a.samples or math.ceil(a.total_samples / chunks)
    print_report(a, *one_config(a, samples, chunks, None))


if __name__ == "__main__":
    main()
