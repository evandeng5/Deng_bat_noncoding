#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys, re, os, json
from datetime import datetime

BASES = set("ACGTacgt")

def parse_species_from_header(header: str) -> str:
    """将 FASTA 标题解析为“纯物种名”"""
    h = header.lstrip('>').strip()
    if not h:
        return ""
    h = h.split()[0]      # >Name xxx -> Name
    h = h.split(':')[0]   # Name:pos  -> Name
    sp = h.split('.')[0]  # Name.chr1 -> Name
    return sp

def read_fasta_by_species(path):
    """返回 {species: seq}, 逐物种聚合"""
    seqs = {}
    name = None; buf = []
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\r\n')        # 重要：去掉 CRLF
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    sp = parse_species_from_header(name)
                    if sp:
                        seqs[sp] = ''.join(buf)
                name = line
                buf = []
            else:
                buf.append(line)
        if name is not None:
            sp = parse_species_from_header(name)
            if sp:
                seqs[sp] = ''.join(buf)
    return seqs

def cov_strict(ref: str, sp: str):
    """只在 ref 为 A/C/G/T 的列上计数（cov ≤ 1）"""
    L = min(len(ref), len(sp))
    den = 0; num = 0
    for i in range(L):
        r = ref[i]; s = sp[i]
        if r in BASES:
            den += 1
            if s in BASES:
                num += 1
    cov = (num / den) if den else 0.0
    return cov, num, den

def cov_species_insert_ok(ref: str, sp: str):
    """分母：ref 非缺口列；分子：sp 所有 A/C/G/T（允许插入，cov 可>1）"""
    L = min(len(ref), len(sp))
    den = sum(ref[i] in BASES for i in range(L))
    num = sum(sp[i]  in BASES for i in range(L))
    cov = (num / den) if den else 0.0
    return cov, num, den

def cov_naive(sp: str, ref_len_from_bed: int):
    """老口径：sp 的 A/C/G/T 总数 / BED 长度"""
    num = sum(c in BASES for c in sp)
    den = ref_len_from_bed
    cov = (num / den) if den else 0.0
    return cov, num, den

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--ref-len', required=True, type=int)
    ap.add_argument('--cov-low',  type=float, default=0.8)
    ap.add_argument('--cov-high', type=float, default=1.2)
    ap.add_argument('--min-outgroup', type=int, default=5)
    ap.add_argument('--min-yin',       type=int, default=2)
    ap.add_argument('--min-yang',      type=int, default=2)
    ap.add_argument('--min-bats-total',type=int, default=10)
    ap.add_argument('--outgroup', nargs='*', default=[])
    ap.add_argument('--yin',      nargs='*', default=[])
    ap.add_argument('--yang',     nargs='*', default=[])
    ap.add_argument('--log-tsv',  required=True)
    ap.add_argument('--ref-species', default='Human')
    ap.add_argument('--coverage-mode', choices=['strict','species_insert_ok','naive'],
                    default='strict')
    ap.add_argument('--per-species-dir', default=None)
    ap.add_argument('--debug', action='store_true', help='写出详细调试信息')
    ap.add_argument('--dump-all-present', action='store_true', default=True,
                    help='per_species 输出所有 present 物种（group=present_but_ungrouped）')
    args = ap.parse_args()

    # 读取 FASTA
    seqs = read_fasta_by_species(args.fasta)
    present = sorted(seqs.keys())

    # 解析坐标（用于命名 debug 文件）
    chrom, start, end = "NA","NA","NA"
    m = re.search(r'Human\.(.+):(\d+)-(\d+)\.fasta$', args.fasta)
    if m:
        chrom, start, end = m.group(1), m.group(2), m.group(3)
    region_tag = f"{chrom}:{start}-{end}" if chrom!="NA" else os.path.basename(args.fasta)

    # 参考序列
    if args.ref_species not in seqs:
        sys.stderr.write(f"[ERROR] ref species {args.ref_species} not in {args.fasta}\n")
        return 98
    ref = seqs[args.ref_species]
    aln_len = len(ref)
    ref_non_gap = sum(c in BASES for c in ref)

    # 构造有效外群集合（Human 永远计入）
    outgroup_eff = sorted(set(args.outgroup) | {args.ref_species})

    expected = sorted(set(args.outgroup + args.yin + args.yang))
    missing = [x for x in expected if x not in present]
    extra   = [x for x in present  if x not in expected]

    # 警告：ref 非缺口列 vs BED 长度
    warn_len = abs(ref_non_gap - args.ref_len) > max(5, int(0.05*args.ref_len))
    if warn_len:
        sys.stderr.write(
            f"[WARN] {os.path.basename(args.fasta)}: bed_len={args.ref_len}, "
            f"ref_non_gap={ref_non_gap}, aln_len={aln_len}\n"
        )

    # 选择口径
    def coverage_for_species(sp_name: str):
        sp = seqs.get(sp_name, '')
        if not sp:
            return 0.0, 0, 0
        if sp_name == args.ref_species:
            den = ref_non_gap
            return (1.0 if den>0 else 0.0), den, den
        if args.coverage_mode == 'strict':
            return cov_strict(ref, sp)
        elif args.coverage_mode == 'species_insert_ok':
            return cov_species_insert_ok(ref, sp)
        else:
            return cov_naive(sp, args.ref_len)

    # 计算各组达标情况
    out_ok = yin_ok = yang_ok = 0
    per_rows = []  # [(group, species, cov, ok)]

    for sp in outgroup_eff:
        cov, _, _ = coverage_for_species(sp)
        ok = int(args.cov_low <= cov <= args.cov_high)
        out_ok += ok
        per_rows.append(("outgroup", sp, cov, ok))

    for sp in args.yin:
        cov, _, _ = coverage_for_species(sp)
        ok = int(args.cov_low <= cov <= args.cov_high)
        yin_ok += ok
        per_rows.append(("yin", sp, cov, ok))

    for sp in args.yang:
        cov, _, _ = coverage_for_species(sp)
        ok = int(args.cov_low <= cov <= args.cov_high)
        yang_ok += ok
        per_rows.append(("yang", sp, cov, ok))

    # 额外：把 FASTA 里存在但不在任何组的物种也列出来，便于排查数组是否没传入
    grouped = set(outgroup_eff) | set(args.yin) | set(args.yang)
    for sp in present:
        if sp not in grouped and args.dump_all_present:
            cov, _, _ = coverage_for_species(sp)
            per_rows.append(("present_but_ungrouped", sp, cov, int(args.cov_low <= cov <= args.cov_high)))

    bats_ok = yin_ok + yang_ok

    ok_outgroup   = (out_ok >= args.min_outgroup)
    ok_yin        = (yin_ok >= args.min_yin)
    ok_yang       = (yang_ok >= args.min_yang)
    ok_bats_total = (bats_ok >= args.min_bats_total)

    verdict, reason, exit_code = "qualified", "all thresholds satisfied", 0
    if not ok_outgroup:
        verdict, reason, exit_code = "fail_outgroup", f"outgroup_ok={out_ok} < {args.min_outgroup}", 10
    elif not ok_yin:
        verdict, reason, exit_code = "fail_yin", f"yin_ok={yin_ok} < {args.min_yin}", 11
    elif not ok_yang:
        verdict, reason, exit_code = "fail_yang", f"yang_ok={yang_ok} < {args.min_yang}", 12
    elif not ok_bats_total:
        verdict, reason, exit_code = "fail_bats_total", f"bats_ok={bats_ok} < {args.min_bats_total}", 13

    # 追加总表
    os.makedirs(os.path.dirname(args.log_tsv) or ".", exist_ok=True)
    with open(args.log_tsv, 'a') as w:
        w.write("\t".join(map(str, [
            chrom, start, end, args.ref_len,
            int(ok_outgroup), int(ok_yin), int(ok_yang), int(ok_bats_total),
            out_ok, yin_ok, yang_ok, bats_ok,
            len(seqs), verdict, reason, args.fasta,
            args.coverage_mode, ref_non_gap
        ])) + "\n")

    # 每物种覆盖率明细
    per_dir = args.per_species_dir or os.path.join(os.path.dirname(args.log_tsv) or ".", "per_species")
    os.makedirs(per_dir, exist_ok=True)
    with open(os.path.join(per_dir, f"{region_tag}.tsv"), "w") as w:
        w.write("group\tspecies\tcoverage\tok\n")
        for g, sp, cov, ok in per_rows:
            w.write(f"{g}\t{sp}\t{cov:.6f}\t{ok}\n")

    # 调试信息
    if args.debug:
        dbg_dir = os.path.join(os.path.dirname(args.log_tsv) or ".", "debug")
        os.makedirs(dbg_dir, exist_ok=True)
        dbg = {
            "timestamp": datetime.now().isoformat(timespec='seconds'),
            "region": region_tag,
            "args": {
                "cov_low": args.cov_low, "cov_high": args.cov_high,
                "min_outgroup": args.min_outgroup, "min_yin": args.min_yin,
                "min_yang": args.min_yang, "min_bats_total": args.min_bats_total,
                "coverage_mode": args.coverage_mode,
                "ref_len": args.ref_len,
                "ref_species": args.ref_species,
                "outgroup_cli": args.outgroup,
                "yin_cli": args.yin,
                "yang_cli": args.yang,
                "outgroup_effective": outgroup_eff,
            },
            "present_species": present,
            "expected_species": expected,
            "missing_expected_in_fasta": missing,
            "extra_present_not_in_expected": extra,
            "ref_non_gap": ref_non_gap,
            "aln_len": aln_len,
            "warn_len": warn_len,
            "counts": {
                "out_ok": out_ok, "yin_ok": yin_ok, "yang_ok": yang_ok, "bats_ok": bats_ok,
                "ok_outgroup": bool(ok_outgroup), "ok_yin": bool(ok_yin),
                "ok_yang": bool(ok_yang), "ok_bats_total": bool(ok_bats_total)
            },
            "verdict": verdict, "reason": reason
        }
        with open(os.path.join(dbg_dir, f"{region_tag}.debug.txt"), "w") as w:
            w.write(json.dumps(dbg, indent=2, ensure_ascii=False))
        # 也同步一份简短 stderr
        sys.stderr.write(f"[DEBUG] {region_tag} out_ok={out_ok} yin_ok={yin_ok} "
                         f"yang_ok={yang_ok} bats_ok={bats_ok} verdict={verdict}\n")

    return exit_code

if __name__ == '__main__':
    sys.exit(main())
