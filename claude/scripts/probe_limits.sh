#!/usr/bin/env bash
# Print the system limits that bound a kmindex/kmtricks build on this host.
# Run on a login node and on a compute node, e.g.:
#   sbatch -p gg -N 1 -t 00:05:00 --wrap "bash probe_limits.sh" -o probe_%j.out
# Output is key=value, one per line, easy to diff between nodes.

set -u

kv() { printf '%s=%s\n' "$1" "$2"; }

kv host "$(hostname)"
kv arch "$(uname -m)"
kv kernel "$(uname -r)"
kv date "$(date -Is)"
kv slurm_job "${SLURM_JOB_ID:-none}"
kv slurm_partition "${SLURM_JOB_PARTITION:-none}"

# open files: soft/hard per-process limit and kernel ceilings
kv nofile_soft "$(ulimit -S -n)"
kv nofile_hard "$(ulimit -H -n)"
kv fs_nr_open "$(cat /proc/sys/fs/nr_open 2>/dev/null || echo NA)"
kv fs_file_max "$(cat /proc/sys/fs/file-max 2>/dev/null || echo NA)"
kv fs_file_nr "$(cut -f1 /proc/sys/fs/file-nr 2>/dev/null || echo NA)"

# mmap count: kmindex merge and query map one file per (thread, sub-index)
kv vm_max_map_count "$(cat /proc/sys/vm/max_map_count 2>/dev/null || echo NA)"

# processes / threads
kv nproc_soft "$(ulimit -S -u)"
kv nproc_hard "$(ulimit -H -u)"
kv threads_max "$(cat /proc/sys/kernel/threads-max 2>/dev/null || echo NA)"

# memory
kv cpus "$(nproc)"
kv mem_total_kb "$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
kv mem_avail_kb "$(awk '/MemAvailable/ {print $2}' /proc/meminfo)"
kv rss_soft "$(ulimit -S -m)"
kv as_soft "$(ulimit -S -v)"
kv cgroup_mem_max "$(cat /sys/fs/cgroup/memory.max 2>/dev/null || cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo NA)"

# disks: space and inodes on the candidate run directories
for d in /tmp "${SCRATCH:-}" "${WORK:-}" "${HOME:-}" "$PWD"; do
  [ -n "$d" ] && [ -d "$d" ] || continue
  tag=$(echo "$d" | tr '/' '_')
  kv "df_bytes_avail$tag" "$(df -B1 --output=avail "$d" 2>/dev/null | tail -1 | tr -d ' ')"
  kv "df_inodes_avail$tag" "$(df --output=iavail "$d" 2>/dev/null | tail -1 | tr -d ' ')"
  kv "fstype$tag" "$(df --output=fstype "$d" 2>/dev/null | tail -1 | tr -d ' ')"
done

# quotas, when the tools exist
command -v lfs >/dev/null 2>&1 && [ -n "${WORK:-}" ] && lfs quota -u "$USER" "$WORK" 2>/dev/null | sed 's/^/lfs_quota_work: /'
command -v quota >/dev/null 2>&1 && quota -s 2>/dev/null | sed 's/^/quota: /'

# binaries
for b in kmindex kmtricks kmhelpers; do
  p=$(command -v "$b" 2>/dev/null || echo NA)
  kv "bin_$b" "$p"
  [ "$p" != NA ] && kv "ver_$b" "$("$b" --version 2>&1 | head -1)"
done

# quick sequential disk throughput on /tmp and $SCRATCH (1 GiB, optional)
if [ "${PROBE_DISK:-0}" = 1 ]; then
  for d in /tmp "${SCRATCH:-}"; do
    [ -n "$d" ] && [ -d "$d" ] || continue
    f="$d/probe_$$.bin"
    w=$( { dd if=/dev/zero of="$f" bs=1M count=1024 oflag=direct 2>&1 || true; } | tail -1)
    r=$( { dd if="$f" of=/dev/null bs=1M iflag=direct 2>&1 || true; } | tail -1)
    rm -f "$f"
    kv "dd_write$(echo "$d" | tr '/' '_')" "$w"
    kv "dd_read$(echo "$d" | tr '/' '_')" "$r"
  done
fi
