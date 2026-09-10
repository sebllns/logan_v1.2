#!/usr/bin/env bash
# Build C chunks with kmindex, merge them, query the result, and record
# wall time, peak RSS, peak open fds, output size, temp file count and temp
# disk peak of every step as one CSV line. A timeline "epoch fds rss_kb" of
# the process tree is written to logs/STEP_NAME.fds (stage profiles).
# Usable directly under SLURM:
#   sbatch -p gg -N 1 -t 08:00:00 --wrap "bash run_experiment.sh -f fofs/ -o run1 -c 4 -t 32 -p 64 -b 3000000000 -q query.fa"
#
# Inputs
#   -f DIR   directory with chunk FOFs named chunk_0.fof, chunk_1.fof, ...
#            (kmtricks format "id : path", see make_fofs.py); ids unique across chunks
#   -o DIR   output directory (index registry + logs + results.csv)
#   -c INT   number of chunks to use (default: all FOFs found)
#   -t INT   threads for build
#   -m INT   threads for merge (default: -t)
#   -p INT   partitions (0 = kmtricks auto)
#   -b INT   bloom size in bits
#   -k INT   k-mer size (default 25)
#   -C       compress intermediate files (kmindex build --cpr)
#   -q FILE  query fasta (optional)
#   -z INT   z value for the query (default 6)
#   -B INT   query batch size (default 0)
#   -T INT   query threads (default 8)
#   -F       query with --fast
#   -n NAME  experiment name written in the CSV (default: basename of -o)
#   -M       skip the merge step (query the chunks together instead)
#   -D       delete chunk indexes after merge (kmindex merge -d)
#   -L       two-level merge: merge pairs of chunks first, then merge the results
#
# CSV columns: exp,step,name,samples,partitions,threads,bf_size,chunks,batch,z,
#   wall_s,max_rss_kb,max_fds,out_bytes,tmp_files_max,tmp_bytes_max,status
set -euo pipefail

fof_dir=""; out=""; chunks=0; threads=32; parts=0; bloom=0; k=25; cpr=0
query=""; z=6; qbatch=0; qthreads=8; mthreads=0; name=""; fast=0; skip_merge=0
delete_old=0; two_level=0
while getopts "f:o:c:t:m:p:b:k:Cq:z:B:T:Fn:MDL" opt; do
  case $opt in
    f) fof_dir=$OPTARG;; o) out=$OPTARG;; c) chunks=$OPTARG;; t) threads=$OPTARG;;
    m) mthreads=$OPTARG;; p) parts=$OPTARG;; b) bloom=$OPTARG;; k) k=$OPTARG;; C) cpr=1;;
    q) query=$OPTARG;; z) z=$OPTARG;; B) qbatch=$OPTARG;; T) qthreads=$OPTARG;; F) fast=1;;
    n) name=$OPTARG;; M) skip_merge=1;; D) delete_old=1;; L) two_level=1;;
    *) exit 2;;
  esac
done
[ -n "$fof_dir" ] && [ -n "$out" ] && [ "$bloom" -gt 0 ] || { echo "need -f -o -b" >&2; exit 2; }
[ "$mthreads" -gt 0 ] || mthreads=$threads
[ -n "$name" ] || name=$(basename "$out")
command -v kmindex >/dev/null || { echo "kmindex not in PATH" >&2; exit 2; }
command -v kmtricks >/dev/null || { echo "kmtricks not in PATH" >&2; exit 2; }

mkdir -p "$out/logs" "$out/tmp"
reg="$out/index"
csv="$out/results.csv"
[ -f "$csv" ] || echo "exp,step,name,samples,partitions,threads,bf_size,chunks,batch,z,wall_s,max_rss_kb,max_fds,out_bytes,tmp_files_max,tmp_bytes_max,status" > "$csv"

mapfile -t fofs < <(ls "$fof_dir"/chunk_*.fof | sort -V)
[ "$chunks" -gt 0 ] && fofs=("${fofs[@]:0:$chunks}")
chunks=${#fofs[@]}
echo "ulimit -n soft=$(ulimit -S -n) hard=$(ulimit -H -n), chunks=$chunks, threads=$threads"

# run "$@" under /usr/bin/time; sample fds and RSS of the process tree every
# second, temp files and temp bytes of WATCHDIR every 10 seconds.
# usage: run STEP NAME SAMPLES THREADS WATCHDIR OUTDIR -- cmd...
run() {
  local step=$1 sname=$2 samples=$3 tcol=$4 watch=$5 outdir=$6; shift 6; shift  # drop --
  local log="$out/logs/${step}_${sname}.log" tlog="$out/logs/${step}_${sname}.time"
  local fdlog="$out/logs/${step}_${sname}.fds"
  echo "[$(date +%T)] $step $sname: $*" | tee -a "$out/logs/commands.log"
  echo "epoch fds rss_kb" > "$fdlog"
  local t0; t0=$(date +%s.%N)
  /usr/bin/time -v -o "$tlog" "$@" > "$log" 2>&1 &
  local tpid=$!
  local maxfd=0 maxtmp=0 maxbytes=0 status=0 tick=0
  sleep 0.5
  local cpid; cpid=$(pgrep -P "$tpid" | head -1 || true)
  while kill -0 "$tpid" 2>/dev/null; do
    if [ -n "$cpid" ] && [ -d "/proc/$cpid/fd" ]; then
      local n=0 rss=0 pid
      # kmindex build forks kmtricks: count the whole tree
      for pid in $cpid $(pgrep -P "$cpid" 2>/dev/null || true); do
        n=$(( n + $(ls "/proc/$pid/fd" 2>/dev/null | wc -l) ))
        rss=$(( rss + $(awk '/VmRSS/ {print $2}' "/proc/$pid/status" 2>/dev/null || echo 0) ))
      done
      [ "$n" -gt "$maxfd" ] && maxfd=$n
      echo "$(date +%s) $n $rss" >> "$fdlog"
    fi
    if [ -n "$watch" ] && [ -d "$watch" ] && [ $((tick % 10)) = 0 ]; then
      local m b
      m=$(find "$watch" -type f 2>/dev/null | wc -l)
      b=$(du -sb "$watch" 2>/dev/null | cut -f1)
      [ "$m" -gt "$maxtmp" ] && maxtmp=$m
      [ "${b:-0}" -gt "$maxbytes" ] && maxbytes=$b
    fi
    tick=$((tick + 1))
    sleep 1
  done
  wait "$tpid" || status=$?
  local wall; wall=$(echo "$(date +%s.%N) - $t0" | bc)
  local rssmax; rssmax=$(awk -F': ' '/Maximum resident/ {print $2}' "$tlog")
  local bytes=0
  [ -d "$outdir" ] && bytes=$(du -sb "$outdir" 2>/dev/null | cut -f1)
  local p="$parts"
  [ "$step" = build ] && [ -f "$outdir/matrices/matrix_0.cmbf" ] && p=$(ls "$outdir/matrices" | wc -l)
  echo "$name,$step,$sname,$samples,$p,$tcol,$bloom,$chunks,$qbatch,$z,$wall,$rssmax,$maxfd,$bytes,$maxtmp,$maxbytes,$status" >> "$csv"
  echo "[$(date +%T)] $step $sname done: wall=${wall}s rss=${rssmax}kB fds=$maxfd tmp_files=$maxtmp tmp_bytes=$maxbytes status=$status"
  return $status
}

count_samples() { grep -c ':' "$1"; }
join_names() { local IFS=,; echo "$*"; }

# 1. build chunks; the first defines the repartition, the others reuse it (--from)
cpropt=(); [ "$cpr" = 1 ] && cpropt=(--cpr)
first=""
for i in "${!fofs[@]}"; do
  fof=${fofs[$i]}; cname="chunk_$i"; dir="$out/tmp/$cname"
  from=(); [ -n "$first" ] && from=(--from "$first")
  run build "$cname" "$(count_samples "$fof")" "$threads" "$dir" "$dir" -- \
    kmindex build -i "$reg" -f "$fof" -d "$dir" -r "$cname" -k "$k" --hard-min 1 \
      --bloom-size "$bloom" --nb-partitions "$parts" -t "$threads" "${cpropt[@]}" "${from[@]}"
  [ -n "$first" ] || first=$cname
done
total=0; for fof in "${fofs[@]}"; do total=$(( total + $(count_samples "$fof") )); done
names=(); for i in "${!fofs[@]}"; do names+=("chunk_$i"); done

# 2. merge
merged="merged"
if [ "$skip_merge" = 0 ] && [ "$chunks" -gt 1 ]; then
  del=(); [ "$delete_old" = 1 ] && del=(-d)
  if [ "$two_level" = 1 ] && [ "$chunks" -gt 2 ]; then
    lvl=(); j=0
    for ((i = 0; i < chunks; i += 2)); do
      if [ $((i + 1)) -lt "$chunks" ]; then
        run merge "pair_$j" "$total" "$mthreads" "$out/tmp" "$out/tmp/pair_$j" -- \
          kmindex merge -i "$reg" -n "pair_$j" -p "$out/tmp/pair_$j" -m "chunk_$i,chunk_$((i + 1))" -t "$mthreads" "${del[@]}"
        lvl+=("pair_$j")
      else
        lvl+=("chunk_$i")
      fi
      j=$((j + 1))
    done
    run merge "$merged" "$total" "$mthreads" "$out/tmp" "$out/tmp/$merged" -- \
      kmindex merge -i "$reg" -n "$merged" -p "$out/tmp/$merged" -m "$(join_names "${lvl[@]}")" -t "$mthreads" "${del[@]}"
  else
    run merge "$merged" "$total" "$mthreads" "$out/tmp" "$out/tmp/$merged" -- \
      kmindex merge -i "$reg" -n "$merged" -p "$out/tmp/$merged" -m "$(join_names "${names[@]}")" -t "$mthreads" "${del[@]}"
  fi
else
  merged="chunk_0"
fi

# 3. query the merged index, and the unmerged chunks together for comparison
if [ -n "$query" ]; then
  fastopt=(); [ "$fast" = 1 ] && fastopt=(--fast)
  run query "$merged" "$total" "$qthreads" "" "$out/q_$merged" -- \
    kmindex query -i "$reg" -n "$merged" -q "$query" -o "$out/q_$merged" -z "$z" -b "$qbatch" -t "$qthreads" "${fastopt[@]}"
  if [ "$skip_merge" = 1 ] && [ "$chunks" -gt 1 ]; then
    run query "all_chunks" "$total" "$qthreads" "" "$out/q_all_chunks" -- \
      kmindex query -i "$reg" -n "$(join_names "${names[@]}")" -q "$query" -o "$out/q_all_chunks" -z "$z" -b "$qbatch" -t "$qthreads" "${fastopt[@]}"
  fi
fi
echo "results in $csv"
