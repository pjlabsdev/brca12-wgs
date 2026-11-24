//! Self-checks of the aligner: every alignment is compared with an exhaustive full-matrix alignment
//! on random sequences, every reported path is scored again, and mapping must not depend on the
//! number of threads.
use super::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n.max(1) as u64) as usize
    }
    /// random sequence of length lo + below(spread)
    fn rseq(&mut self, lo: usize, spread: usize) -> Vec<u8> {
        let n = lo + self.below(spread);
        self.seq(n)
    }
    fn seq(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[self.below(4)]).collect()
    }
    /// substitutions, small and large insertions and deletions at the given rate per base
    fn mutate(&mut self, s: &[u8], rate_per_mille: usize, max_indel: usize) -> Vec<u8> {
        let mut o = Vec::with_capacity(s.len() + 16);
        let mut i = 0;
        while i < s.len() {
            if self.below(1000) < rate_per_mille {
                match self.below(4) {
                    0 | 1 => {
                        o.push(b"ACGT"[(code(s[i]) as usize + 1 + self.below(3)) % 4]);
                        i += 1;
                    }
                    2 => {
                        let l = 1 + self.below(max_indel);
                        o.extend(self.seq(l));
                    }
                    _ => i += 1 + self.below(max_indel),
                }
            } else {
                o.push(s[i]);
                i += 1;
            }
        }
        o
    }
}

/// number of random cases; HLA_PROPS=n runs n times as many (a longer sweep before a release)
fn reps(n: usize) -> usize {
    n * std::env::var("HLA_PROPS").ok().and_then(|v| v.parse().ok()).unwrap_or(1)
}

fn subst(sc: &Score, a: u8, b: u8) -> i32 {
    if a == b'N' || b == b'N' {
        -1
    } else if a == b {
        sc.a
    } else {
        -sc.b
    }
}

fn gap(sc: &Score, l: i32) -> i32 {
    (sc.q + l * sc.e).min(sc.q2 + l * sc.e2)
}

/// exhaustive two-piece affine alignment from (0, 0): the full matrix H
fn full_matrix(q: &[u8], r: &[u8], sc: &Score) -> Vec<Vec<i32>> {
    let (n, m) = (q.len(), r.len());
    let neg = i32::MIN / 4;
    let mut h = vec![vec![neg; m + 1]; n + 1];
    let mut e1 = vec![vec![neg; m + 1]; n + 1];
    let mut e2 = e1.clone();
    let mut f1 = e1.clone();
    let mut f2 = e1.clone();
    h[0][0] = 0;
    for i in 0..=n {
        for j in 0..=m {
            if i == 0 && j == 0 {
                continue;
            }
            if j > 0 {
                e1[i][j] = (h[i][j - 1] - sc.q - sc.e).max(e1[i][j - 1] - sc.e);
                e2[i][j] = (h[i][j - 1] - sc.q2 - sc.e2).max(e2[i][j - 1] - sc.e2);
            }
            if i > 0 {
                f1[i][j] = (h[i - 1][j] - sc.q - sc.e).max(f1[i - 1][j] - sc.e);
                f2[i][j] = (h[i - 1][j] - sc.q2 - sc.e2).max(f2[i - 1][j] - sc.e2);
            }
            let d = if i > 0 && j > 0 { h[i - 1][j - 1] + subst(sc, q[i - 1], r[j - 1]) } else { neg };
            h[i][j] = d.max(e1[i][j]).max(e2[i][j]).max(f1[i][j]).max(f2[i][j]);
        }
    }
    h
}

/// score of a path, and the bases of q and r it uses; panics on an op that does not fit the bases
fn path_score(q: &[u8], r: &[u8], ops: &[(u8, u32)], sc: &Score) -> (i32, usize, usize) {
    let (mut i, mut j, mut s) = (0usize, 0usize, 0i32);
    for &(op, l) in ops {
        let l = l as usize;
        match op {
            b'=' | b'X' => {
                for _ in 0..l {
                    assert_eq!(op == b'=', q[i] == r[j], "op {} at q {i} r {j}", op as char);
                    s += subst(sc, q[i], r[j]);
                    i += 1;
                    j += 1;
                }
            }
            b'I' => {
                s -= gap(sc, l as i32);
                i += l;
            }
            b'D' => {
                s -= gap(sc, l as i32);
                j += l;
            }
            _ => panic!("op {op}"),
        }
    }
    // two gaps of one kind next to each other would have been merged by push()
    for w in ops.windows(2) {
        assert_ne!(w[0].0, w[1].0, "unmerged ops");
    }
    (s, i, j)
}

fn presets() -> [(&'static str, Score); 4] {
    [("sr", Opts::sr().sc), ("asm10", Opts::asm10().sc), ("asm5", Opts::asm5().sc), ("hifi", Opts::hifi().sc)]
}

#[test]
fn global_is_optimal_on_small_pieces() {
    let mut g = Rng(1);
    for (name, sc) in presets() {
        for t in 0..reps(400) {
            let q = g.rseq(0, 120);
            let r = match t % 3 {
                0 => g.rseq(0, 120),
                1 => g.mutate(&q, 30, 6),
                _ => g.mutate(&q, 100, 25),
            };
            let (s, qn, rn, ops) = dp_raw(&q, &r, &sc, false);
            let best = full_matrix(&q, &r, &sc)[q.len()][r.len()];
            assert_eq!(s, best, "{name} case {t}: score {s}, optimum {best}");
            assert_eq!((qn, rn), (q.len(), r.len()));
            assert_eq!(path_score(&q, &r, &ops, &sc), (s, q.len(), r.len()), "{name} case {t}: path");
        }
    }
}

#[test]
fn banded_global_is_optimal_when_the_path_fits_the_band() {
    // pieces over 100k cells use a band of the length difference plus 64 around the diagonal
    let mut g = Rng(2);
    for (name, sc) in presets() {
        for t in 0..reps(20) {
            let q = g.rseq(400, 500);
            let r = g.mutate(&q, 8, 8);
            let (s, _, _, ops) = dp_raw(&q, &r, &sc, false);
            let best = full_matrix(&q, &r, &sc)[q.len()][r.len()];
            assert_eq!(s, best, "{name} case {t} ({}x{})", q.len(), r.len());
            assert_eq!(path_score(&q, &r, &ops, &sc).0, s);
        }
    }
}

#[test]
fn extension_finds_the_best_end() {
    // without z-drop and with a band wider than r, the extension must end at the best cell of the
    // full matrix (with the bonus for the end of q) and report the best score anywhere
    let mut g = Rng(3);
    for (name, mut sc) in presets() {
        sc.zdrop = 1 << 20;
        sc.bw = 1000;
        if sc.end_bonus == 0 {
            sc.end_bonus = 10;
        }
        for t in 0..reps(300) {
            let q = g.rseq(1, 150);
            let mut r = match t % 3 {
                0 => g.mutate(&q, 20, 5),
                1 => g.mutate(&q, 80, 10),
                _ => g.rseq(0, 160),
            };
            r.extend(g.rseq(0, 30));
            let h = full_matrix(&q, &r, &sc);
            let any = h.iter().flatten().copied().max().unwrap();
            let at_end = h[q.len()].iter().copied().max().unwrap() + sc.end_bonus;
            let (s, qn, rn, ops) = dp_raw(&q, &r, &sc, true);
            assert_eq!(s, any, "{name} case {t}: score");
            let (ps, pi, pj) = path_score(&q, &r, &ops, &sc);
            assert_eq!((pi, pj), (qn, rn), "{name} case {t}: path end");
            assert_eq!(ps, h[qn][rn], "{name} case {t}: path not optimal for its end");
            let got = ps + if qn == q.len() { sc.end_bonus } else { 0 };
            assert_eq!(got, any.max(at_end), "{name} case {t}: end chosen {qn},{rn} q {} r {}", String::from_utf8_lossy(&q), String::from_utf8_lossy(&r));
        }
    }
}

#[test]
fn split_fill_matches_the_full_alignment() {
    // large gaps between anchors are split at unique 13-mers; on allele-like pairs this must cost nothing
    let mut g = Rng(4);
    let sc = Opts::asm10().sc;
    let (mut same, mut total) = (0, 0);
    for t in 0..reps(12) {
        let q = g.rseq(600, 600);
        let r = g.mutate(&q, if t % 2 == 0 { 5 } else { 15 }, 20);
        let (s, ops) = fill(&q, &r, &sc, 0);
        assert_eq!(path_score(&q, &r, &ops, &sc), (s, q.len(), r.len()), "case {t}: path");
        let best = full_matrix(&q, &r, &sc)[q.len()][r.len()];
        assert!(s <= best);
        total += 1;
        if s == best {
            same += 1;
        } else {
            eprintln!("fill case {t}: {s} vs optimum {best}");
        }
    }
    assert!(same * 10 >= total * 9, "split fill optimal in {same} of {total}");
}

#[test]
fn left_align_keeps_the_score_and_moves_gaps_left() {
    let mut g = Rng(5);
    let sc = Opts::sr().sc;
    for t in 0..reps(300) {
        // repeats make gap positions ambiguous
        let unit = g.rseq(1, 4);
        let mut q = g.seq(20);
        for _ in 0..3 + g.below(6) {
            q.extend_from_slice(&unit);
        }
        q.extend(g.seq(20));
        let r = g.mutate(&q, 40, 3);
        let (s, _, _, ops) = dp_raw(&q, &r, &sc, false);
        let la = left_align(&ops, &q, &r, 0, 0);
        assert_eq!(path_score(&q, &r, &la, &sc), (s, q.len(), r.len()), "case {t}");
        assert_eq!(left_align(&la, &q, &r, 0, 0), la, "case {t}: not stable");
        // no single gap can move one base further left
        let (mut i, mut j) = (0usize, 0usize);
        for (k, &(op, l)) in la.iter().enumerate() {
            let l = l as usize;
            if (op == b'I' || op == b'D') && k > 0 && la[k - 1].0 == b'=' {
                let last = if op == b'D' { r[j + l - 1] } else { q[i + l - 1] };
                let before = if op == b'D' { r[j - 1] } else { q[i - 1] };
                assert_ne!(last, before, "case {t}: gap at q {i} r {j} could move left");
            }
            match op {
                b'=' | b'X' => {
                    i += l;
                    j += l;
                }
                b'I' => i += l,
                _ => j += l,
            }
        }
    }
}

/// a gene family: one ancestor, alleles with few differences each
fn family(g: &mut Rng, len: usize, alleles: usize) -> (Vec<String>, Vec<Vec<u8>>) {
    let base = g.seq(len);
    let mut names = vec![];
    let mut seqs = vec![];
    for a in 0..alleles {
        names.push(format!("G*{:02}", a + 1));
        seqs.push(if a == 0 { base.clone() } else { g.mutate(&base, 6, 3) });
    }
    (names, seqs)
}

#[test]
fn reads_map_to_their_allele_and_paths_hold() {
    let mut g = Rng(6);
    let (names, seqs) = family(&mut g, 3000, 12);
    for win in [0usize, 11] {
        let idx = Index::new(names.clone(), seqs.clone(), Opts { win, full: 0.8, per_group: 25, ..Opts::sr() });
        let (mut ok, mut n) = (0, 0);
        for t in 0..reps(300) {
            let a = g.below(seqs.len());
            let s = &seqs[a];
            let p = g.below(s.len() - 150);
            let mut read = g.mutate(&s[p..p + 150], 3, 2);
            let fwd = t % 2 == 0;
            if !fwd {
                read = revcomp(&read);
            }
            let hs = idx.map(&read);
            let Some(top) = hs.first() else { panic!("read {t} without hit") };
            for h in &hs {
                let qo = if h.fwd { read.clone() } else { revcomp(&read) };
                let (qa, qb) = if h.fwd { (h.q_st, h.q_en) } else { (read.len() - h.q_en, read.len() - h.q_st) };
                let (ps, pi, pj) = path_score(&qo[qa..qb], &seqs[h.seq][h.r_st..h.r_en], &h.ops, &Opts::sr().sc);
                assert_eq!((pi, pj), (qb - qa, h.r_en - h.r_st), "read {t}: path length");
                assert!(ps <= h.score, "read {t}: path {ps} over score {}", h.score);
                assert_eq!(h.nm, h.ops.iter().filter(|o| o.0 != b'=').map(|o| o.1 as usize).sum::<usize>());
            }
            n += 1;
            let mine = hs.iter().filter(|h| h.seq == a).map(|h| h.score).max();
            assert!(top.fwd == fwd, "read {t}: strand");
            if mine == Some(top.score) {
                ok += 1;
            }
        }
        // a read of an allele scores best on its own allele (ties with alleles equal there count)
        assert!(ok * 100 >= n * 99, "win {win}: own allele best for {ok} of {n} reads");
    }
}

#[test]
fn mapping_is_the_same_with_any_number_of_threads() {
    let mut g = Rng(7);
    let (names, seqs) = family(&mut g, 2500, 10);
    let idx = Index::new(names, seqs.clone(), Opts { win: 11, full: 0.8, per_group: 25, ..Opts::sr() });
    let reads: Vec<Vec<u8>> = (0..400)
        .map(|_| {
            let s = &seqs[g.below(seqs.len())];
            let p = g.below(s.len() - 150);
            g.mutate(&s[p..p + 150], 10, 4)
        })
        .collect();
    let key = |hs: Vec<Hit>| hs.iter().map(|h| (h.seq, h.fwd, h.q_st, h.q_en, h.r_st, h.r_en, h.score, h.ops.clone())).collect::<Vec<_>>();
    let one: Vec<_> = reads.iter().map(|r| key(idx.map(r))).collect();
    for th in [2usize, 4, 8] {
        let out: Vec<_> = std::thread::scope(|s| {
            let hs: Vec<_> = (0..th)
                .map(|k| {
                    let (idx, reads) = (&idx, &reads);
                    s.spawn(move || reads.iter().enumerate().filter(|(i, _)| i % th == k).map(|(i, r)| (i, key(idx.map(r)))).collect::<Vec<_>>())
                })
                .collect();
            let mut v: Vec<_> = hs.into_iter().flat_map(|h| h.join().unwrap()).collect();
            v.sort_by_key(|x| x.0);
            v.into_iter().map(|x| x.1).collect()
        });
        assert!(out == one, "{th} threads differ");
    }
}

#[test]
fn memo_gives_the_same_result() {
    let mut g = Rng(8);
    let sc = Opts::asm10().sc;
    for _ in 0..4 {
        let q = g.seq(700);
        let r = g.mutate(&q, 10, 10);
        let a = dp_raw(&q, &r, &sc, false);
        assert_eq!(dp(&q, &r, &sc, false), a);
        assert_eq!(dp(&q, &r, &sc, false), a);
        let a = dp_raw(&q, &r, &sc, true);
        assert_eq!(dp(&q, &r, &sc, true), a);
    }
}

#[test]
fn quality_lowers_only_the_cost_of_uncertain_mismatches() {
    let mut g = Rng(9);
    let r = g.seq(400);
    let idx = Index::single(&r, Opts::sr_single());
    let mut q = r[100..250].to_vec();
    for p in [20usize, 140] {
        q[p] = if q[p] == b'A' { b'C' } else { b'A' };
    }
    let b = Opts::sr().sc.b;
    let mut qual = vec![30u8; 150];
    let h = idx.map(&q)[0].clone();
    assert_eq!(qual_score(&h, &qual, b, 20), h.score, "high qualities change nothing");
    qual[140] = 2;
    assert_eq!(qual_score(&h, &qual, b, 20), h.score + b - b * 2 / 20);
    // on the reverse strand the same base keeps its quality
    let rq = revcomp(&q);
    let rqual: Vec<u8> = qual.iter().rev().copied().collect();
    let hr = idx.map(&rq)[0].clone();
    assert!(!hr.fwd);
    assert_eq!(qual_score(&hr, &rqual, b, 20), hr.score + b - b * 2 / 20);
    assert_eq!(qual_score(&h, &[], b, 20), h.score);
}

/// speed of the extension alignment: cargo test --release bench_ext -- --ignored --nocapture
#[test]
#[ignore]
fn bench_ext() {
    let mut g = Rng(10);
    let sc = Opts::sr().sc;
    let cases: Vec<(Vec<u8>, Vec<u8>)> = (0..reps(20000))
        .map(|_| {
            let q = g.rseq(60, 80);
            let mut r = g.mutate(&q, 40, 3);
            r.extend(g.rseq(32, 1));
            (q, r)
        })
        .collect();
    let t = std::time::Instant::now();
    let (mut cells, mut sum) = (0usize, 0i64);
    for (q, r) in &cases {
        let o = dp_raw(q, r, &sc, true);
        cells += q.len() * (2 * sc.bw + 1).min(r.len());
        sum += o.0 as i64;
    }
    let e = t.elapsed().as_secs_f64();
    eprintln!("ext: {:.1} ms, {:.2} ns per band cell, checksum {sum}", e * 1e3, e * 1e9 / cells as f64);
}

#[test]
fn banded_local_matches_local() {
    let mut g = Rng(11);
    let sc = Opts::sr().sc;
    for t in 0..reps(300) {
        let r = g.rseq(150, 200);
        let p = g.below(r.len() - 100);
        let q = g.mutate(&r[p..(p + 150).min(r.len())], 30, 6);
        let full = local(&q, &r, &sc);
        // a band over every cell is the full problem
        assert_eq!(local_band(&q, &r, -(q.len() as i64), r.len() as i64, &sc), full, "case {t}: wide band");
        // a band around the true diagonal with room for the indels finds the same score
        let d = p as i64;
        let b = local_band(&q, &r, d - 40, d + 40, &sc);
        assert!(b.0 <= full.0);
        assert_eq!(b.0, full.0, "case {t}: band around the diagonal");
    }
}

/// speed of the mate rescue: cargo test --release bench_rescue -- --ignored --nocapture
#[test]
#[ignore]
fn bench_rescue() {
    let mut g = Rng(12);
    let sc = Opts::sr().sc;
    let cases: Vec<(Vec<u8>, Vec<u8>)> = (0..reps(3000))
        .map(|_| {
            let r = g.rseq(1750, 1);
            let p = g.below(1500);
            let q = g.mutate(&r[p..p + 150], 20, 3);
            (q, r)
        })
        .collect();
    let t = std::time::Instant::now();
    let mut sum = 0i64;
    for (q, r) in &cases {
        sum += local_seeded(q, r, &sc).0 as i64;
    }
    let e = t.elapsed().as_secs_f64();
    eprintln!("rescue: {:.1} us per call, checksum {sum}", e * 1e6 / cases.len() as f64);
    let t = std::time::Instant::now();
    for (q, r) in &cases {
        let p: i64 = 0;
        sum += local_band(q, &r[..230.min(r.len())], p - 40, p + 40, &sc).0 as i64;
    }
    let e = t.elapsed().as_secs_f64();
    eprintln!("band only: {:.1} us per call", e * 1e6 / cases.len() as f64);
}

/// banded local alignment cell by cell: (best, i, j) with the smallest (i, j) among the best cells
fn band_ref(q: &[u8], r: &[u8], dlo: i64, dhi: i64, sc: &Score) -> (i32, usize, usize) {
    let (n, m) = (q.len(), r.len());
    let neg = i32::MIN / 4;
    let mut h = vec![vec![0i32; m + 1]; n + 1];
    let (mut e1, mut e2, mut f1, mut f2) = (vec![vec![neg; m + 1]; n + 1], vec![vec![neg; m + 1]; n + 1], vec![vec![neg; m + 1]; n + 1], vec![vec![neg; m + 1]; n + 1]);
    let mut best = (0, 0, 0);
    for i in 1..=n {
        for j in 1..=m {
            let d = j as i64 - i as i64;
            if d < dlo || d > dhi {
                continue;
            }
            e1[i][j] = (h[i][j - 1] - sc.q - sc.e).max(e1[i][j - 1] - sc.e);
            e2[i][j] = (h[i][j - 1] - sc.q2 - sc.e2).max(e2[i][j - 1] - sc.e2);
            f1[i][j] = (h[i - 1][j] - sc.q - sc.e).max(f1[i - 1][j] - sc.e);
            f2[i][j] = (h[i - 1][j] - sc.q2 - sc.e2).max(f2[i - 1][j] - sc.e2);
            let v = (h[i - 1][j - 1] + subst(sc, q[i - 1], r[j - 1])).max(e1[i][j]).max(e2[i][j]).max(f1[i][j]).max(f2[i][j]).max(0);
            h[i][j] = v;
            if v > best.0 {
                best = (v, i, j);
            }
        }
    }
    best
}

#[test]
fn antidiagonal_band_matches_cell_by_cell() {
    let mut g = Rng(13);
    let sc = Opts::sr().sc;
    for t in 0..reps(400) {
        let r = g.rseq(20, 250);
        let q = if t % 3 == 0 { g.rseq(5, 160) } else { let p = g.below(r.len()); g.mutate(&r[p..(p + 150).min(r.len())], 40, 6) };
        if q.is_empty() {
            continue;
        }
        let c = g.below(r.len()) as i64 - g.below(q.len()) as i64;
        let w = 2 + g.below(60) as i64;
        let (dlo, dhi) = (c - w, c + w);
        assert_eq!(local_band_ad(&q, &r, dlo, dhi, &sc), band_ref(&q, &r, dlo, dhi, &sc), "case {t} n {} m {} band {dlo}..{dhi}", q.len(), r.len());
    }
}

#[test]
fn antidiagonal_extension_equals_row_wise() {
    let mut g = Rng(14);
    for (name, base) in presets() {
        for t in 0..reps(500) {
            let mut sc = base;
            sc.bw = [8usize, 16, 32, 100, 500][t % 5];
            sc.zdrop = [20, 50, 100, 400][(t / 5) % 4];
            if t % 7 == 0 {
                sc.end_bonus = 0;
            }
            let q = g.rseq(1, 180);
            let mut r = match t % 4 {
                0 => g.mutate(&q, 20, 4),
                1 => g.mutate(&q, 120, 12),
                2 => g.rseq(0, 200),
                _ => {
                    let mut x = g.mutate(&q[..q.len() / 2], 10, 3);
                    x.extend(g.rseq(10, 100));
                    x
                }
            };
            r.extend(g.rseq(0, 40));
            if t % 11 == 0 && !r.is_empty() {
                let p = g.below(r.len());
                r[p] = b'N';
            }
            let a = dp_raw(&q, &r, &sc, true);
            FORCE_ROWS.with(|c| c.set(true));
            let b = dp_raw(&q, &r, &sc, true);
            FORCE_ROWS.with(|c| c.set(false));
            assert_eq!(a, b, "{name} case {t} n {} m {} bw {} zdrop {}", q.len(), r.len(), sc.bw, sc.zdrop);
        }
    }
}

#[test]
fn antidiagonal_global_equals_row_wise() {
    let mut g = Rng(15);
    for (name, sc) in presets() {
        for t in 0..reps(300) {
            let q = g.rseq(1, if t % 3 == 0 { 900 } else { 200 });
            let mut r = match t % 3 {
                0 => g.mutate(&q, 10, 10),
                1 => g.mutate(&q, 60, 20),
                _ => g.rseq(1, 250),
            };
            if t % 13 == 0 && !r.is_empty() {
                let p = g.below(r.len());
                r[p] = b'N';
            }
            let a = dp_raw(&q, &r, &sc, false);
            FORCE_ROWS.with(|c| c.set(true));
            let b = dp_raw(&q, &r, &sc, false);
            FORCE_ROWS.with(|c| c.set(false));
            assert_eq!(a, b, "{name} case {t} n {} m {}", q.len(), r.len());
        }
    }
}
