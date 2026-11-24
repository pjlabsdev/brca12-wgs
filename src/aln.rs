//! Own sequence aligner in Rust (taken over from hla-wgs src/aln.rs, main 7b22b0f): exact k-mer anchors, chaining of the anchors, and an alignment with
//! two-piece affine gap costs between the anchors (global) and beyond them (extension with z-drop).
//! It replaces the C library minimap2; scores follow the minimap2 presets, so thresholds keep their
//! meaning.
use ahash::AHashMap;

#[derive(Clone, Copy, Debug)]
pub struct Score {
    pub a: i32,
    pub b: i32,
    /// gap of length l costs min(q + l*e, q2 + l*e2)
    pub q: i32,
    pub e: i32,
    pub q2: i32,
    pub e2: i32,
    pub zdrop: i32,
    /// bonus when an extension reaches the end of the query (not part of the reported score)
    pub end_bonus: i32,
    /// band beyond the query length for end extensions
    pub bw: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct Opts {
    pub k: usize,
    pub sc: Score,
    /// k-mers seen more often than this in the index are not used as anchors
    pub max_occ: usize,
    /// largest distance between two chained anchors
    pub max_gap: i64,
    /// largest difference of the diagonals of two chained anchors
    pub chain_bw: i64,
    pub min_chain: f32,
    /// at most this many chains are aligned per query strand pair
    pub max_align: usize,
    /// chains below this fraction of the best chain are not aligned
    pub chain_ratio: f32,
    pub min_score: i32,
    /// at most this many chains per group of sequences (gene: name before '*'); 0 = no limit
    pub per_group: usize,
    /// k-mers that occur more often than this within one sequence give no anchors on it (repeats)
    pub max_occ_seq: usize,
    /// how many earlier anchors chaining looks at
    pub lookback: usize,
    /// only every stride-th query k-mer is looked up (long reads)
    pub stride: usize,
    /// minimizer window: of every win consecutive query k-mers only the one with the smallest hash
    /// is looked up (0 or 1: all); the index keeps every k-mer, so each looked-up k-mer finds all its matches
    pub win: usize,
    /// with minimizers: reseed with every k-mer when the best chain scores below this fraction of the read length
    pub full: f32,
}

impl Opts {
    pub fn sr() -> Opts {
        Opts { k: 21, sc: Score { a: 2, b: 8, q: 12, e: 2, q2: 24, e2: 1, zdrop: 100, end_bonus: 10, bw: 32 }, max_occ: 5000, max_gap: 100, chain_bw: 100, min_chain: 25.0, max_align: 20, chain_ratio: 0.0, min_score: 40, per_group: 0, max_occ_seq: 50, lookback: 50, stride: 1, win: 0, full: 0.0 }
    }
    /// short reads against a single sequence: shorter seeds find reads of new alleles with dense differences
    pub fn sr_single() -> Opts {
        Opts { k: 15, max_occ: 50, ..Opts::sr() }
    }
    pub fn asm10() -> Opts {
        Opts { k: 19, sc: Score { a: 1, b: 9, q: 16, e: 2, q2: 41, e2: 1, zdrop: 200, end_bonus: 0, bw: 500 }, max_occ: 50, max_gap: 10000, chain_bw: 2000, min_chain: 40.0, max_align: 5, chain_ratio: 0.0, min_score: 40, per_group: 0, max_occ_seq: 1, lookback: 200, stride: 1, win: 0, full: 0.0 }
    }
    pub fn asm5() -> Opts {
        Opts { sc: Score { a: 1, b: 19, q: 39, e: 3, q2: 81, e2: 1, zdrop: 200, end_bonus: 0, bw: 500 }, ..Opts::asm10() }
    }
    pub fn hifi() -> Opts {
        Opts { k: 19, sc: Score { a: 1, b: 4, q: 6, e: 2, q2: 26, e2: 1, zdrop: 400, end_bonus: 0, bw: 100 }, max_occ: 5000, max_gap: 5000, chain_bw: 5000, min_chain: 40.0, max_align: 50, chain_ratio: 0.0, min_score: 40, per_group: 0, max_occ_seq: 2, lookback: 200, stride: 4, win: 0, full: 0.0 }
    }
}

#[derive(Clone, Debug)]
pub struct Hit {
    pub seq: usize,
    pub fwd: bool,
    /// query interval on the forward query
    pub q_st: usize,
    pub q_en: usize,
    pub r_st: usize,
    pub r_en: usize,
    pub score: i32,
    /// matching bases
    pub mlen: usize,
    /// mismatches plus inserted and deleted bases
    pub nm: usize,
    /// (op, length) with op one of = X I D, query in target orientation
    pub ops: Vec<(u8, u32)>,
    pub primary: bool,
}

impl Hit {
    /// difference string in minimap2 short form (":n", "*rq", "+ins", "-del"), lower case
    pub fn cs(&self, q: &[u8], r: &[u8]) -> String {
        let qo: Vec<u8> = if self.fwd { q.to_vec() } else { revcomp(q) };
        let (mut i, mut j) = (if self.fwd { self.q_st } else { q.len() - self.q_en }, self.r_st);
        let mut s = String::new();
        for &(op, l) in &self.ops {
            let l = l as usize;
            match op {
                b'=' => {
                    s.push_str(&format!(":{l}"));
                    i += l;
                    j += l;
                }
                b'X' => {
                    for _ in 0..l {
                        s.push('*');
                        s.push(r[j].to_ascii_lowercase() as char);
                        s.push(qo[i].to_ascii_lowercase() as char);
                        i += 1;
                        j += 1;
                    }
                }
                b'I' => {
                    s.push('+');
                    s.push_str(&String::from_utf8_lossy(&qo[i..i + l]).to_lowercase());
                    i += l;
                }
                _ => {
                    s.push('-');
                    s.push_str(&String::from_utf8_lossy(&r[j..j + l]).to_lowercase());
                    j += l;
                }
            }
        }
        s
    }
}

/// Score of a hit when mismatches at uncertain bases cost less: a mismatch at a base of quality q
/// costs b·min(q, q0)/q0 instead of b. A base the sequencer was unsure about says little about which
/// allele or gene the read comes from; at q ≥ q0 nothing changes. qual is the forward read's.
pub fn qual_score(h: &Hit, qual: &[u8], b: i32, q0: u8) -> i32 {
    let n = qual.len();
    if n == 0 || q0 == 0 || h.q_en > n {
        return h.score;
    }
    // position in the query as aligned (target orientation)
    let mut i = if h.fwd { h.q_st } else { n - h.q_en };
    let mut add = 0;
    for &(op, l) in &h.ops {
        let l = l as usize;
        match op {
            b'X' => {
                for k in i..i + l {
                    let q = qual[if h.fwd { k } else { n - 1 - k }].min(q0) as i32;
                    add += b - b * q / q0 as i32;
                }
                i += l;
            }
            b'=' | b'I' => i += l,
            _ => {}
        }
    }
    h.score + add
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

fn code(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// invertible integer hash (Thomas Wang), orders k-mers for minimizers without favouring low-complexity ones
fn hash64(mut x: u64) -> u64 {
    x = (!x).wrapping_add(x << 21);
    x ^= x >> 24;
    x = x.wrapping_add(x << 3).wrapping_add(x << 8);
    x ^= x >> 14;
    x = x.wrapping_add(x << 2).wrapping_add(x << 4);
    x ^= x >> 28;
    x.wrapping_add(x << 31)
}

/// k-mers ending at each position (value, end position); positions with N are skipped
fn kmers(s: &[u8], k: usize, mut f: impl FnMut(u64, usize)) {
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut v, mut l) = (0u64, 0usize);
    for (i, &c) in s.iter().enumerate() {
        let x = code(c);
        if x > 3 {
            l = 0;
            continue;
        }
        v = ((v << 2) | x as u64) & mask;
        l += 1;
        if l >= k {
            f(v, i);
        }
    }
}

pub struct Index {
    pub opts: Opts,
    pub names: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
    map: AHashMap<u64, (u32, u32)>,
    pos: Vec<(u32, u32)>,
    group: Vec<u32>,
}

impl Index {
    pub fn new(names: Vec<String>, seqs: Vec<Vec<u8>>, opts: Opts) -> Index {
        let seqs: Vec<Vec<u8>> = seqs.into_iter().map(|s| s.to_ascii_uppercase()).collect();
        let mut all: Vec<(u64, u32, u32)> = Vec::with_capacity(seqs.iter().map(|s| s.len()).sum());
        for (si, s) in seqs.iter().enumerate() {
            kmers(s, opts.k, |v, i| all.push((v, si as u32, i as u32)));
        }
        all.sort_unstable();
        let mut map = AHashMap::with_capacity(all.len() / 4);
        let mut pos = Vec::with_capacity(all.len());
        let mut i = 0;
        while i < all.len() {
            let mut j = i;
            while j < all.len() && all[j].0 == all[i].0 {
                j += 1;
            }
            map.insert(all[i].0, (pos.len() as u32, (j - i) as u32));
            pos.extend(all[i..j].iter().map(|x| (x.1, x.2)));
            i = j;
        }
        let mut gid: AHashMap<&str, u32> = AHashMap::new();
        let group = names.iter().map(|n| { let g = n.split('*').next().unwrap_or(""); let l = gid.len() as u32; *gid.entry(g).or_insert(l) }).collect();
        Index { opts, names, seqs, map, pos, group }
    }

    pub fn single(seq: &[u8], opts: Opts) -> Index {
        Index::new(vec![String::new()], vec![seq.to_vec()], opts)
    }


    /// chains of anchors of the query against all sequences, best first: (score, seq, fwd, anchors (q end, r end))
    fn chains(&self, q: &[u8], win: usize) -> Vec<(f32, usize, bool, Vec<(u32, u32)>)> {
        let k = self.opts.k;
        let ns = self.seqs.len();
        let mut out = vec![];
        let mut tm = if *STAT_ON { Some(std::time::Instant::now()) } else { None };
        let rc = revcomp(q);
        // k-mer hits per strand: (query end, start in pos, count)
        let mut hl: [Vec<(u32, u32, u32)>; 2] = [vec![], vec![]];
        for (si, qs) in [q, rc.as_slice()].into_iter().enumerate() {
            let stride = self.opts.stride.max(1);
            let mut look = |v: u64, i: usize| {
                if let Some(&(st, n)) = self.map.get(&v) {
                    if (n as usize) <= self.opts.max_occ {
                        hl[si].push((i as u32, st, n));
                    }
                }
            };
            if win > 1 {
                let mut km: Vec<(u64, u64, usize)> = vec![];
                kmers(qs, k, |v, i| km.push((hash64(v), v, i)));
                // the first smallest hash of each window of `win` k-mers, with a monotone queue
                // (the earliest of equal hashes stays in front, as a scan from the left finds it)
                let w = win;
                let nwin = km.len().saturating_sub(w - 1).max(1.min(km.len()));
                let mut dq: std::collections::VecDeque<usize> = std::collections::VecDeque::with_capacity(w + 1);
                let mut last = usize::MAX;
                let mut next = 0;
                for s in 0..nwin {
                    let e = (s + w).min(km.len());
                    while next < e {
                        while dq.back().is_some_and(|&b| km[b].0 > km[next].0) {
                            dq.pop_back();
                        }
                        dq.push_back(next);
                        next += 1;
                    }
                    while dq.front().is_some_and(|&b| b < s) {
                        dq.pop_front();
                    }
                    let b = *dq.front().unwrap();
                    if b != last {
                        look(km[b].1, km[b].2);
                        last = b;
                    }
                }
            } else {
                kmers(qs, k, |v, i| {
                    if i % stride == 0 {
                        look(v, i);
                    }
                });
            }
        }
        part(0, &mut tm);
        // with many sequences per group only the sequences with the most anchors per group are chained
        let chosen: Option<Vec<bool>> = if self.opts.per_group > 0 && ns > 1 {
            let mut cnt = vec![0u32; ns * 2];
            for si in 0..2 {
                for &(_, st, n) in &hl[si] {
                    for &(s, _) in &self.pos[st as usize..(st + n) as usize] {
                        cnt[s as usize * 2 + si] += 1;
                    }
                }
            }
            let top = *cnt.iter().max().unwrap_or(&0);
            let min = ((top as f32) * self.opts.chain_ratio).max(1.0) as u32;
            let mut cand: Vec<usize> = (0..ns * 2).filter(|&x| cnt[x] >= min).collect();
            cand.sort_by(|&a, &b| cnt[b].cmp(&cnt[a]).then(a.cmp(&b)));
            // per group the sequences with the most anchors; all that tie with the best of the group are
            // kept as well (up to four times the limit), since equal counts say nothing about which aligns best
            let mut per: AHashMap<u32, (usize, u32)> = AHashMap::new();
            let mut ch = vec![false; ns * 2];
            for x in cand {
                let e = per.entry(self.group[x / 2]).or_insert((0, cnt[x]));
                if e.0 < self.opts.per_group || (cnt[x] + 1 >= e.1 && e.0 < 4 * self.opts.per_group) {
                    e.0 += 1;
                    ch[x] = true;
                }
            }
            Some(ch)
        } else {
            None
        };
        part(1, &mut tm);
        let mut bid = vec![u32::MAX; ns];
        for si in 0..2 {
            // anchors per sequence: (r end, q end); k-mers come in query order, so each bucket is nearly sorted
            let mut buckets: Vec<(u32, Vec<(u32, u32)>)> = vec![];
            for &(i, st, n) in &hl[si] {
                let ps = &self.pos[st as usize..(st + n) as usize];
                let mut a = 0;
                while a < ps.len() {
                    let mut b = a + 1;
                    while b < ps.len() && ps[b].0 == ps[a].0 {
                        b += 1;
                    }
                    let s = ps[a].0;
                    if b - a <= self.opts.max_occ_seq && chosen.as_ref().map_or(true, |c| c[s as usize * 2 + si]) {
                        let mut x = bid[s as usize];
                        if x == u32::MAX {
                            x = buckets.len() as u32;
                            bid[s as usize] = x;
                            buckets.push((s, Vec::with_capacity(hl[si].len())));
                        }
                        buckets[x as usize].1.extend(ps[a..b].iter().map(|&(_, p)| (p, i)));
                    }
                    a = b;
                }
            }
            for (s, _) in &buckets {
                bid[*s as usize] = u32::MAX;
            }
            buckets.sort_unstable_by_key(|b| b.0);
            part(2, &mut tm);
            // alleles of one gene give the same anchors, shifted by where the stretch starts in each
            // sequence: chaining depends only on distances, so the chains are computed once per shape
            let mut memo: AHashMap<Vec<(u32, u32)>, usize> = AHashMap::new();
            let mut shapes: Vec<Vec<(f32, Vec<(u32, u32)>)>> = vec![];
            for (s, mut v) in buckets {
                v.sort_unstable();
                let r0 = v[0].0;
                let g: Vec<(u32, u32)> = v.iter().map(|x| (x.1, x.0 - r0)).collect();
                let id = match memo.get(&g) {
                    Some(&id) => id,
                    None => {
                        shapes.push(self.chain(&g));
                        memo.insert(g, shapes.len() - 1);
                        shapes.len() - 1
                    }
                };
                for (sc, c) in &shapes[id] {
                    out.push((*sc, s as usize, si == 0, c.iter().map(|&(a, b)| (a, b + r0)).collect()));
                }
            }
            part(3, &mut tm);
        }
        out.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)).then(b.2.cmp(&a.2)));
        out
    }

    /// anchors (q end, r end) sorted by r then q; chains in the manner of minimap2
    fn chain(&self, a: &[(u32, u32)]) -> Vec<(f32, Vec<(u32, u32)>)> {
        let k = self.opts.k as f32;
        let n = a.len();
        let mut f = vec![0f32; n];
        let mut p = vec![usize::MAX; n];
        for i in 0..n {
            f[i] = k;
            let (qi, ri) = (a[i].0 as i64, a[i].1 as i64);
            // the previous anchor on the same diagonal right before this one is the best predecessor
            if i > 0 {
                let (qj, rj) = (a[i - 1].0 as i64, a[i - 1].1 as i64);
                let (dq, dr) = (qi - qj, ri - rj);
                if dq == dr && dq > 0 && dq <= self.opts.k as i64 {
                    f[i] = f[i - 1] + dq as f32;
                    p[i] = i - 1;
                    continue;
                }
            }
            let lo = i.saturating_sub(self.opts.lookback);
            for j in (lo..i).rev() {
                let (qj, rj) = (a[j].0 as i64, a[j].1 as i64);
                let (dq, dr) = (qi - qj, ri - rj);
                if dr > self.opts.max_gap {
                    break;
                }
                if dq <= 0 || dr <= 0 || dq > self.opts.max_gap {
                    continue;
                }
                let dd = (dr - dq).abs();
                if dd > self.opts.chain_bw {
                    continue;
                }
                let gain = dq.min(dr).min(self.opts.k as i64) as f32;
                let cost = if dd == 0 { 0.0 } else { 0.01 * k * dd as f32 + 0.5 * (dd as f32).log2() };
                let s = f[j] + gain - cost;
                if s > f[i] {
                    f[i] = s;
                    p[i] = j;
                }
            }
        }
        let mut ord: Vec<usize> = (0..n).collect();
        ord.sort_by(|&x, &y| f[y].partial_cmp(&f[x]).unwrap().then(x.cmp(&y)));
        let mut used = vec![false; n];
        let mut out = vec![];
        for &e in &ord {
            if used[e] {
                continue;
            }
            let mut c = vec![];
            let mut i = e;
            let mut stop = 0f32;
            loop {
                if used[i] {
                    stop = f[i];
                    break;
                }
                used[i] = true;
                c.push(a[i]);
                if p[i] == usize::MAX {
                    break;
                }
                i = p[i];
            }
            let sc = f[e] - stop;
            if sc >= self.opts.min_chain {
                c.reverse();
                out.push((sc, c));
            }
            if out.len() >= 4 {
                break;
            }
        }
        out
    }

    /// all alignments of the query, best score first; the first is marked primary
    pub fn map(&self, q: &[u8]) -> Vec<Hit> {
        let q = q.to_ascii_uppercase();
        let t = if *STAT_ON { Some(std::time::Instant::now()) } else { None };
        let mut ch = self.chains(&q, self.opts.win);
        // minimizers find a read that matches a sequence almost exactly; a read whose best chain covers
        // less than `full` of it (variants, a new allele, a gene boundary) is seeded again with every k-mer
        if self.opts.win > 1 && ch.first().map_or(true, |c| c.0 < self.opts.full * q.len() as f32) {
            ch = self.chains(&q, 0);
        }
        if let Some(t) = t {
            CHAIN_NS.fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        }
        let best = ch.first().map(|c| c.0).unwrap_or(0.0);
        let rc = revcomp(&q);
        let mut per: AHashMap<u32, usize> = AHashMap::new();
        let mut cache: AHashMap<Vec<u8>, (usize, Option<Hit>)> = AHashMap::new();
        let mut hits: Vec<Hit> = ch
            .iter()
            .filter(|c| c.0 >= best * self.opts.chain_ratio)
            .filter(|c| {
                if self.opts.per_group == 0 {
                    return true;
                }
                let n = per.entry(self.group[c.1]).or_insert(0);
                *n += 1;
                *n <= 4 * self.opts.per_group
            })
            .take(self.opts.max_align)
            .filter_map(|(_, s, fwd, c)| {
                let qs = if *fwd { q.as_slice() } else { rc.as_slice() };
                let r = &self.seqs[*s];
                // the alignment only depends on the stretch of r it can reach and on the anchors in it:
                // sequences that agree there (alleles of one gene) are aligned once
                let k = self.opts.k;
                let (q0, r0) = (c[0].0 as usize + 1 - k, c[0].1 as usize + 1 - k);
                let (ql, rl) = (c[c.len() - 1].0 as usize + 1, c[c.len() - 1].1 as usize + 1);
                let ws = r0 - r0.min(q0 + self.opts.sc.bw);
                let we = (rl + (qs.len() - ql + self.opts.sc.bw)).min(r.len());
                let mut key: Vec<u8> = Vec::with_capacity(we - ws + 8 * c.len() + 1);
                key.push(*fwd as u8);
                for &(a, b) in c.iter() {
                    key.extend_from_slice(&a.to_le_bytes());
                    key.extend_from_slice(&(b - ws as u32).to_le_bytes());
                }
                key.extend_from_slice(&r[ws..we]);
                let mut h = match cache.get(&key) {
                    Some((w0, h0)) => {
                        let mut h = h0.clone()?;
                        h.r_st = h.r_st - w0 + ws;
                        h.r_en = h.r_en - w0 + ws;
                        h
                    }
                    None => {
                        let h = align_chain(qs, r, c, k, &self.opts.sc);
                        cache.insert(key, (ws, h.clone()));
                        h?
                    }
                };
                if h.score < self.opts.min_score {
                    return None;
                }
                h.seq = *s;
                h.fwd = *fwd;
                if !*fwd {
                    let (a, b) = (q.len() - h.q_en, q.len() - h.q_st);
                    h.q_st = a;
                    h.q_en = b;
                }
                Some(h)
            })
            .collect();
        hits.sort_by(|a, b| b.score.cmp(&a.score).then(a.seq.cmp(&b.seq)).then(b.fwd.cmp(&a.fwd)).then(a.r_st.cmp(&b.r_st)));
        // the same alignment reached by two chains is reported once
        hits.dedup_by(|a, b| a.seq == b.seq && a.fwd == b.fwd && a.r_st == b.r_st && a.r_en == b.r_en && a.q_st == b.q_st);
        if let Some(h) = hits.first_mut() {
            h.primary = true;
        }
        hits
    }
}

/// Best local alignment of q (in target orientation) anywhere in r, two-piece affine gaps, score only:
/// (score, q start, q end, r start, r end). Used to rescue a mate next to its partner.
pub fn local(q: &[u8], r: &[u8], sc: &Score) -> (i32, usize, usize, usize, usize) {
    let fwd = |q: &[u8], r: &[u8]| -> (i32, usize, usize) {
        let m = r.len();
        let (oe1, oe2) = (sc.q + sc.e, sc.q2 + sc.e2);
        let mut hp = vec![0i32; m + 1];
        let mut hc = vec![0i32; m + 1];
        let mut f1 = vec![NEG; m + 1];
        let mut f2 = vec![NEG; m + 1];
        let (mut best, mut bi, mut bj) = (0, 0, 0);
        let (a, b, e_1, e_2) = (sc.a, sc.b, sc.e, sc.e2);
        for i in 1..=q.len() {
            let (mut e1, mut e2) = (NEG, NEG);
            hc[0] = 0;
            let qi = q[i - 1];
            let qn = qi == b'N';
            // row best, first cell when several tie (as a scan of the cells would find it)
            let (mut rb, mut rj) = (0, 0);
            let mut hleft = 0;
            let (hpr, hcr, f1r, f2r) = (&hp[..=m], &mut hc[..=m], &mut f1[..=m], &mut f2[..=m]);
            for j in 1..=m {
                let hpj = hpr[j];
                let n1 = (hpj - oe1).max(f1r[j] - e_1);
                let n2 = (hpj - oe2).max(f2r[j] - e_2);
                f1r[j] = n1;
                f2r[j] = n2;
                e1 = (hleft - oe1).max(e1 - e_1);
                e2 = (hleft - oe2).max(e2 - e_2);
                let x = r[j - 1];
                let d = hpr[j - 1] + if qn || x == b'N' { -1 } else if qi == x { a } else { -b };
                let h = d.max(e1).max(e2).max(n1).max(n2).max(0);
                hcr[j] = h;
                hleft = h;
                if h > rb {
                    rb = h;
                    rj = j;
                }
            }
            if rb > best {
                best = rb;
                bi = i;
                bj = rj;
            }
            std::mem::swap(&mut hp, &mut hc);
        }
        (best, bi, bj)
    };
    let (s, qe, re) = fwd(q, r);
    if s == 0 {
        return (0, 0, 0, 0, 0);
    }
    let qr: Vec<u8> = q[..qe].iter().rev().copied().collect();
    let rr: Vec<u8> = r[..re].iter().rev().copied().collect();
    let (_, ql, rl) = fwd(&qr, &rr);
    (s, qe - ql, qe, re - rl, re)
}

/// `local` restricted to the cells with dlo <= j - i <= dhi (i in q, j in r, both from 1): an
/// alignment along a known diagonal with indels up to the band's margin. Cells outside the band count
/// as the start of a new alignment, which is exact for the banded problem.
pub fn local_band(q: &[u8], r: &[u8], dlo: i64, dhi: i64, sc: &Score) -> (i32, usize, usize, usize, usize) {
    let fwd = |q: &[u8], r: &[u8], dlo: i64, dhi: i64| -> (i32, usize, usize) {
        let (n, m) = (q.len(), r.len());
        let (oe1, oe2) = (sc.q + sc.e, sc.q2 + sc.e2);
        let (a, b, e_1, e_2) = (sc.a, sc.b, sc.e, sc.e2);
        let mut hp = vec![0i32; m + 2];
        let mut hc = vec![0i32; m + 2];
        let mut f1 = vec![NEG; m + 2];
        let mut f2 = vec![NEG; m + 2];
        let (mut best, mut bi, mut bj) = (0, 0, 0);
        for i in 1..=n {
            let jl = (i as i64 + dlo).max(1);
            let jh = (i as i64 + dhi).min(m as i64);
            if jl > jh {
                if jl > m as i64 {
                    break;
                }
                continue;
            }
            let (jl, jh) = (jl as usize, jh as usize);
            let qi = q[i - 1];
            let qn = qi == b'N';
            let (mut e1, mut e2) = (NEG, NEG);
            let mut hleft = 0;
            let (mut rb, mut rj) = (0, 0);
            for j in jl..=jh {
                let hpj = hp[j];
                let n1 = (hpj - oe1).max(f1[j] - e_1);
                let n2 = (hpj - oe2).max(f2[j] - e_2);
                f1[j] = n1;
                f2[j] = n2;
                e1 = (hleft - oe1).max(e1 - e_1);
                e2 = (hleft - oe2).max(e2 - e_2);
                let x = r[j - 1];
                let d = hp[j - 1] + if qn || x == b'N' { -1 } else if qi == x { a } else { -b };
                let h = d.max(e1).max(e2).max(n1).max(n2).max(0);
                hc[j] = h;
                hleft = h;
                if h > rb {
                    rb = h;
                    rj = j;
                }
            }
            // the cell after the band end is read by the next row as its upper neighbour: outside
            // the band it is a fresh start
            hc[jh + 1] = 0;
            if jl >= 1 {
                hc[jl - 1] = 0;
            }
            if rb > best {
                best = rb;
                bi = i;
                bj = rj;
            }
            std::mem::swap(&mut hp, &mut hc);
        }
        (best, bi, bj)
    };
    let fwd = |q: &[u8], r: &[u8], dlo: i64, dhi: i64| if *RESCUE_AD && dhi > dlo { local_band_ad(q, r, dlo, dhi, sc) } else { fwd(q, r, dlo, dhi) };
    let (s, qe, re) = fwd(q, r, dlo, dhi);
    if s == 0 {
        return (0, 0, 0, 0, 0);
    }
    let qr: Vec<u8> = q[..qe].iter().rev().copied().collect();
    let rr: Vec<u8> = r[..re].iter().rev().copied().collect();
    let off = re as i64 - qe as i64;
    let (_, ql, rl) = fwd(&qr, &rr, off - dhi, off - dlo);
    (s, qe - ql, qe, re - rl, re)
}

static RESCUE_AD: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var("HLA_RESCUE_AD").map_or(true, |v| v != "0"));

/// The forward pass of `local_band` computed by anti-diagonals: the cells of one anti-diagonal do
/// not depend on each other, so the inner loop has no carried dependency and the compiler turns it
/// into vector instructions (SSE/AVX on x86, NEON on ARM). Same result as the row-wise pass,
/// including which of several equally good end cells is reported (smallest i, then smallest j).
fn local_band_ad(q: &[u8], r: &[u8], dlo: i64, dhi: i64, sc: &Score) -> (i32, usize, usize) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU has AVX2
        return unsafe { local_band_ad_avx2(q, r, dlo, dhi, sc) };
    }
    local_band_ad_impl(q, r, dlo, dhi, sc)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn local_band_ad_avx2(q: &[u8], r: &[u8], dlo: i64, dhi: i64, sc: &Score) -> (i32, usize, usize) {
    local_band_ad_impl(q, r, dlo, dhi, sc)
}

#[inline(always)]
fn local_band_ad_impl(q: &[u8], r: &[u8], dlo: i64, dhi: i64, sc: &Score) -> (i32, usize, usize) {
    let (n, m) = (q.len(), r.len());
    if n == 0 || m == 0 {
        return (0, 0, 0);
    }
    let (oe1, oe2, e_1, e_2) = (sc.q + sc.e, sc.q2 + sc.e2, sc.e, sc.e2);
    let (sa, sb) = (sc.a, -sc.b);
    let rrev: Vec<u8> = r.iter().rev().copied().collect();
    // H of anti-diagonals t-2, t-1, t; E in place; F of t-1 and t; all indexed by i
    let mut h0 = vec![0i32; n + 2];
    let mut h1 = vec![0i32; n + 2];
    let mut h2 = vec![0i32; n + 2];
    let mut ea = vec![NEG; n + 2];
    let mut eb = vec![NEG; n + 2];
    let mut fa_p = vec![NEG; n + 2];
    let mut fb_p = vec![NEG; n + 2];
    let mut fa_c = vec![NEG; n + 2];
    let mut fb_c = vec![NEG; n + 2];
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    for t in 2..=(n + m) as i64 {
        let ilo = 1i64.max(t - m as i64).max((t - dhi + 1).div_euclid(2));
        let ihi = (n as i64).min(t - 1).min((t - dlo).div_euclid(2));
        if ilo > ihi {
            if ilo > n as i64 || t - ihi.max(1) > m as i64 + n as i64 {
                break;
            }
            continue;
        }
        let (ilo, ihi) = (ilo as usize, ihi as usize);
        let k0 = (m as i64 - t) as isize;
        let len = ihi - ilo + 1;
        {
            let qs = &q[ilo - 1..ihi];
            let rs = &rrev[(k0 + ilo as isize) as usize..(k0 + ilo as isize) as usize + len];
            let hd = &h0[ilo - 1..ihi];
            let hu = &h1[ilo - 1..ihi]; // (i-1, j): upper neighbour, index i-1 on t-1
            let hl = &h1[ilo..=ihi]; // (i, j-1): left neighbour, index i on t-1
            let fap = &fa_p[ilo - 1..ihi];
            let fbp = &fb_p[ilo - 1..ihi];
            let ear = &mut ea[ilo..=ihi];
            let ebr = &mut eb[ilo..=ihi];
            let fac = &mut fa_c[ilo..=ihi];
            let fbc = &mut fb_c[ilo..=ihi];
            let ho = &mut h2[ilo..=ihi];
            for x in 0..len {
                let (a, b) = (qs[x], rs[x]);
                let sv = if a == b'N' || b == b'N' { -1 } else if a == b { sa } else { sb };
                let e1 = (hl[x] - oe1).max(ear[x] - e_1);
                let e2 = (hl[x] - oe2).max(ebr[x] - e_2);
                let f1 = (hu[x] - oe1).max(fap[x] - e_1);
                let f2 = (hu[x] - oe2).max(fbp[x] - e_2);
                ear[x] = e1;
                ebr[x] = e2;
                fac[x] = f1;
                fbc[x] = f2;
                ho[x] = (hd[x] + sv).max(e1).max(e2).max(f1).max(f2).max(0);
            }
        }
        let row = &h2[ilo..=ihi];
        let mx = row.iter().copied().max().unwrap_or(0);
        if mx > 0 && mx >= best {
            let x = row.iter().position(|&v| v == mx).unwrap();
            let (ci, cj) = (ilo + x, t as usize - ilo - x);
            if mx > best || (ci, cj) < (bi, bj) {
                best = mx;
                bi = ci;
                bj = cj;
            }
        }
        // the cell before this anti-diagonal's start is read by the next two as a neighbour outside the band
        h2[ilo - 1] = 0;
        fa_c[ilo - 1] = NEG;
        fb_c[ilo - 1] = NEG;
        std::mem::swap(&mut h0, &mut h1);
        std::mem::swap(&mut h1, &mut h2);
        std::mem::swap(&mut fa_p, &mut fa_c);
        std::mem::swap(&mut fb_p, &mut fb_c);
    }
    (best, bi, bj)
}

/// `local`, restricted to the stretches of r around the two diagonals with the most shared 11-mers
pub fn local_seeded(q: &[u8], r: &[u8], sc: &Score) -> (i32, usize, usize, usize, usize) {
    let t = if *STAT_ON { Some(std::time::Instant::now()) } else { None };
    let x = local_seeded_in(q, r, sc);
    if let Some(t) = t {
        RESCUE[0].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        RESCUE[1].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
    }
    x
}
static RESCUE_VOTES: std::sync::LazyLock<u32> = std::sync::LazyLock::new(|| std::env::var("HLA_RESCUE_V").ok().and_then(|v| v.parse().ok()).unwrap_or(3));
/// indel margin of the banded rescue alignment (HLA_RESCUE_BAND; 0 = the whole window)
static RESCUE_BAND: std::sync::LazyLock<i64> = std::sync::LazyLock::new(|| std::env::var("HLA_RESCUE_BAND").ok().and_then(|v| v.parse().ok()).unwrap_or(24));
static RESCUE: [std::sync::atomic::AtomicU64; 2] = [const { std::sync::atomic::AtomicU64::new(0) }; 2];

fn local_seeded_in(q: &[u8], r: &[u8], sc: &Score) -> (i32, usize, usize, usize, usize) {
    const K: usize = 11;
    // the read's 11-mers are few: index them (sorted, with a 4096-bit filter in front) and scan r,
    // instead of indexing the longer r
    let mut qk: Vec<(u64, u32)> = Vec::with_capacity(q.len());
    kmers(q, K, |v, i| qk.push((v, i as u32)));
    qk.sort_unstable();
    let fbit = |v: u64| (v.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 52) as usize;
    let mut filt = [0u64; 64];
    for &(v, _) in &qk {
        let b = fbit(v);
        filt[b / 64] |= 1 << (b % 64);
    }
    let qrange = |v: u64| -> &[(u64, u32)] {
        let a = qk.partition_point(|x| x.0 < v);
        let b = a + qk[a..].partition_point(|x| x.0 == v);
        &qk[a..b]
    };
    let mut hit: Vec<(u64, u32)> = vec![];
    kmers(r, K, |v, p| {
        let b = fbit(v);
        if filt[b / 64] >> (b % 64) & 1 == 1 && !qrange(v).is_empty() {
            hit.push((v, p as u32));
        }
    });
    hit.sort_unstable();
    // per diagonal bucket of 8: votes and the range of exact diagonals (r end - q end)
    let base = q.len() as i64 / 8 + 2;
    let nb = (r.len() + q.len()) / 8 + 4;
    let mut votes = vec![0u32; nb];
    let mut span = vec![(i64::MAX, i64::MIN); nb];
    let mut a = 0;
    while a < hit.len() {
        let mut b = a + 1;
        while b < hit.len() && hit[b].0 == hit[a].0 {
            b += 1;
        }
        // 11-mers seen more than 8 times in r (repeats) do not vote
        if b - a <= 8 {
            let qs = qrange(hit[a].0);
            for &(_, p) in &hit[a..b] {
                for &(_, i) in qs {
                    let dg = p as i64 - i as i64;
                    let x = (dg.div_euclid(8) + base) as usize;
                    votes[x] += 1;
                    span[x].0 = span[x].0.min(dg);
                    span[x].1 = span[x].1.max(dg);
                }
            }
        }
        a = b;
    }
    let mut v: Vec<(i64, u32)> = votes.iter().enumerate().filter(|x| *x.1 > 0).map(|(x, &c)| (x as i64 - base, c)).collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    // a stretch that shares fewer than three 11-mers with the read does not hold a usable alignment
    // of it (on HG00438: 0 of 282 such rescues reached the score a hit needs); skip the alignment
    if v.first().map_or(true, |x| x.1 < *RESCUE_VOTES) {
        return (0, 0, 0, 0, 0);
    }
    let mut best = (0, 0, 0, 0, 0);
    // the two best diagonal buckets; a bucket next to an earlier one is part of the same diagonal
    // (indels and the bucket border split it) and is aligned together with it
    let mut done: Vec<i64> = vec![];
    for &(d, nv) in v.iter() {
        if done.len() >= 2 {
            break;
        }
        if done.iter().any(|&x| (x - d).abs() <= 1) {
            continue;
        }
        done.push(d);
        // a perfect local alignment cannot be beaten by the second window
        if best.0 >= sc.a * q.len() as i32 || nv < *RESCUE_VOTES {
            break;
        }
        let x = if *RESCUE_BAND > 0 {
            // the diagonals of the shared 11-mers of this bucket and its neighbours, widened by the
            // indel margin, in the cell coordinates of the window (j - i with both from 1)
            let (mut a0, mut a1) = (i64::MAX, i64::MIN);
            for e in d - 1..=d + 1 {
                if let Some(&(x0, x1)) = span.get((e + base).max(0) as usize) {
                    if x0 <= x1 {
                        a0 = a0.min(x0);
                        a1 = a1.max(x1);
                    }
                }
            }
            let lo = (a0 - *RESCUE_BAND).max(0) as usize;
            let hi = ((a1 + q.len() as i64 + *RESCUE_BAND).max(0) as usize).min(r.len());
            if lo >= hi {
                continue;
            }
            let x = local_band(q, &r[lo..hi], a0 - lo as i64 - *RESCUE_BAND, a1 - lo as i64 + *RESCUE_BAND, sc);
            (x.0, x.1, x.2, x.3 + lo, x.4 + lo)
        } else {
            let lo = (d * 8 - 24).max(0) as usize;
            let hi = ((d * 8 + 32 + q.len() as i64 + 24).max(0) as usize).min(r.len());
            if lo >= hi {
                continue;
            }
            let x = local(q, &r[lo..hi], sc);
            (x.0, x.1, x.2, x.3 + lo, x.4 + lo)
        };
        if x.0 > best.0 {
            best = x;
        }
    }
    best
}

/// primary hit and supplementary hits: best first, then every hit that covers a part of the query not yet covered
pub fn primary_set(hits: &[Hit]) -> Vec<usize> {
    let mut keep: Vec<usize> = vec![];
    for (i, h) in hits.iter().enumerate() {
        let len = (h.q_en - h.q_st).max(1);
        let ov = keep.iter().map(|&j| (h.q_en.min(hits[j].q_en) as i64 - h.q_st.max(hits[j].q_st) as i64).max(0) as usize).max().unwrap_or(0);
        if keep.is_empty() || ov * 2 < len {
            keep.push(i);
        }
    }
    keep
}

fn push(ops: &mut Vec<(u8, u32)>, op: u8, l: u32) {
    if l == 0 {
        return;
    }
    if let Some(x) = ops.last_mut() {
        if x.0 == op {
            x.1 += l;
            return;
        }
    }
    ops.push((op, l));
}

/// alignment of q (in target orientation) along a chain of anchors (q end, r end) of length k
fn align_chain(q: &[u8], r: &[u8], c: &[(u32, u32)], k: usize, sc: &Score) -> Option<Hit> {
    // exact blocks [qs, qe) / [rs, re) from the anchors
    let mut blocks: Vec<(usize, usize, usize)> = vec![]; // q start, r start, length
    for &(qe, re) in c {
        let (qe, re) = (qe as usize + 1, re as usize + 1);
        let (qs, rs) = (qe - k, re - k);
        if let Some(b) = blocks.last_mut() {
            let (bq, br, bl) = *b;
            if qs as i64 - rs as i64 == bq as i64 - br as i64 && qs <= bq + bl {
                b.2 = qe - bq;
                continue;
            }
            // a block that overlaps the previous one is trimmed at its start
            let sh = (bq + bl).saturating_sub(qs).max((br + bl).saturating_sub(rs));
            if sh >= k {
                continue;
            }
            blocks.push((qs + sh, rs + sh, k - sh));
        } else {
            blocks.push((qs, rs, k));
        }
    }
    let (q0, r0, _) = blocks[0];
    let mut ops: Vec<(u8, u32)> = vec![];
    let mut score = 0;
    // left extension on the reversed sequences
    let rl = r0.min(q0 + sc.bw);
    let ql: Vec<u8> = q[..q0].iter().rev().copied().collect();
    let rlv: Vec<u8> = r[r0 - rl..r0].iter().rev().copied().collect();
    let (s, qn, rn, lops) = dp(&ql, &rlv, sc, true);
    score += s;
    for &(op, l) in lops.iter().rev() {
        push(&mut ops, op, l);
    }
    let (q_st, r_st) = (q0 - qn, r0 - rn);
    for (bi, &(bq, br, bl)) in blocks.iter().enumerate() {
        if bi > 0 {
            let (pq, pr, pl) = blocks[bi - 1];
            let (gq, gr) = (&q[pq + pl..bq], &r[pr + pl..br]);
            let (s, gops) = fill(gq, gr, sc, 0);
            score += s;
            for &(op, l) in &gops {
                push(&mut ops, op, l);
            }
        }
        score += sc.a * bl as i32;
        push(&mut ops, b'=', bl as u32);
    }
    let (lq, lr, ll) = *blocks.last().unwrap();
    let (qe, re) = (lq + ll, lr + ll);
    let rr = (r.len() - re).min(q.len() - qe + sc.bw);
    let (s, qn, rn, rops) = dp(&q[qe..], &r[re..re + rr], sc, true);
    score += s;
    for &(op, l) in &rops {
        push(&mut ops, op, l);
    }
    let ops = left_align(&ops, q, r, q_st, r_st);
    let (mut mlen, mut nm) = (0, 0);
    for &(op, l) in &ops {
        match op {
            b'=' => mlen += l as usize,
            _ => nm += l as usize,
        }
    }
    Some(Hit { seq: 0, fwd: true, q_st, q_en: qe + qn, r_st, r_en: re + rn, score, mlen, nm, ops, primary: false })
}

/// Global alignment of a stretch between two anchors. A large stretch is split at exact 13-mers that
/// occur once on both sides and lie in order, so only the pieces between them need the full alignment.
fn fill(q: &[u8], r: &[u8], sc: &Score, depth: u32) -> (i32, Vec<(u8, u32)>) {
    const K: usize = 13;
    if (q.len() as u64) * (r.len() as u64) <= 250_000 || depth >= 2 || q.len() < 2 * K || r.len() < 2 * K {
        let (s, _, _, o) = dp(q, r, sc, false);
        return (s, o);
    }
    let count = |x: &[u8]| {
        let mut m: AHashMap<u64, (u32, u32)> = AHashMap::new();
        kmers(x, K, |v, i| {
            let e = m.entry(v).or_insert((0, i as u32));
            e.0 += 1;
        });
        m
    };
    let (mq, mr) = (count(q), count(r));
    let mut an: Vec<(u32, u32)> = mq.iter().filter(|(_, x)| x.0 == 1).filter_map(|(v, x)| mr.get(v).filter(|y| y.0 == 1).map(|y| (x.1, y.1))).collect();
    an.sort_unstable();
    // longest chain increasing in both, with a diagonal not too far from the stretch's own
    let n = an.len();
    let mut f = vec![1u32; n];
    let mut p = vec![usize::MAX; n];
    for i in 0..n {
        for j in i.saturating_sub(64)..i {
            if an[j].0 < an[i].0 && an[j].1 < an[i].1 && f[j] + 1 > f[i] {
                f[i] = f[j] + 1;
                p[i] = j;
            }
        }
    }
    let Some(mut e) = (0..n).max_by_key(|&i| (f[i], std::cmp::Reverse(i))) else {
        let (s, _, _, o) = dp(q, r, sc, false);
        return (s, o);
    };
    let mut ch = vec![];
    loop {
        ch.push(an[e]);
        if p[e] == usize::MAX {
            break;
        }
        e = p[e];
    }
    ch.reverse();
    // pieces: before the first 13-mer, between them, after the last; 13-mers that overlap the previous are skipped
    let (mut qi, mut rj) = (0usize, 0usize);
    let (mut score, mut ops) = (0, vec![]);
    for &(a, b) in &ch {
        let (a, b) = (a as usize + 1 - K, b as usize + 1 - K);
        if a < qi || b < rj {
            continue;
        }
        let (s, o) = fill(&q[qi..a], &r[rj..b], sc, depth + 1);
        score += s + sc.a * K as i32;
        for (op, l) in o {
            push(&mut ops, op, l);
        }
        push(&mut ops, b'=', K as u32);
        qi = a + K;
        rj = b + K;
    }
    let (s, o) = fill(&q[qi..], &r[rj..], sc, depth + 1);
    score += s;
    for (op, l) in o {
        push(&mut ops, op, l);
    }
    (score, ops)
}

/// gaps in repeats are moved to their leftmost equivalent place, so that the same event gets the
/// same position in every read
fn left_align(ops: &[(u8, u32)], q: &[u8], r: &[u8], q0: usize, r0: usize) -> Vec<(u8, u32)> {
    if !ops.iter().any(|o| o.0 == b'I' || o.0 == b'D') {
        return ops.to_vec();
    }
    let mut ex: Vec<u8> = ops.iter().flat_map(|&(o, l)| std::iter::repeat(o).take(l as usize)).collect();
    let (mut qi, mut rj) = (Vec::with_capacity(ex.len()), Vec::with_capacity(ex.len()));
    let (mut i, mut j) = (q0, r0);
    for &o in &ex {
        qi.push(i);
        rj.push(j);
        if o != b'D' { i += 1; }
        if o != b'I' { j += 1; }
    }
    let mut k = 0;
    while k < ex.len() {
        let o = ex[k];
        if o != b'I' && o != b'D' {
            k += 1;
            continue;
        }
        let (mut s, mut e) = (k, k);
        while e < ex.len() && ex[e] == o {
            e += 1;
        }
        let next = e;
        while s > 0 && ex[s - 1] == b'=' {
            let same = if o == b'D' { r[rj[s - 1]] == r[rj[e - 1]] } else { q[qi[s - 1]] == q[qi[e - 1]] };
            if !same {
                break;
            }
            ex[s - 1] = o;
            ex[e - 1] = b'=';
            s -= 1;
            e -= 1;
        }
        k = next;
    }
    let mut out = vec![];
    for o in ex {
        push(&mut out, o, 1);
    }
    out
}

static STAT_ON: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var("HLA_ALNSTAT").is_ok());
static STAT: [[std::sync::atomic::AtomicU64; 4]; 3] = [const { [const { std::sync::atomic::AtomicU64::new(0) }; 4] }; 3];
static CHAIN_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
/// time in chains(): seeding, candidate counting, anchor collection and sorting, chaining (ns)
static PART_NS: [std::sync::atomic::AtomicU64; 4] = [const { std::sync::atomic::AtomicU64::new(0) }; 4];
fn part(k: usize, t: &mut Option<std::time::Instant>) {
    if let Some(t0) = t {
        let now = std::time::Instant::now();
        PART_NS[k].fetch_add((now - *t0).as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        *t0 = now;
    }
}

/// dp statistics (env HLA_ALNSTAT): calls, cells, seconds, calls over a million cells, for extensions and gap fills
pub fn stat() -> String {
    let g = |k: usize, i: usize| STAT[k][i].load(std::sync::atomic::Ordering::Relaxed);
    format!("long ext {} calls {} Mcells {:.1}s {} big; ext {} calls {} Mcells {:.1}s {} big; fill {} calls {} Mcells {:.1}s {} big; chaining {:.1}s (seeds {:.1}s, count {:.1}s, anchors {:.1}s, chain {:.1}s); rescue {} calls {:.1}s", g(2, 0), g(2, 1) / 1_000_000, g(2, 2) as f64 / 1e9, g(2, 3), g(0, 0), g(0, 1) / 1_000_000, g(0, 2) as f64 / 1e9, g(0, 3), g(1, 0), g(1, 1) / 1_000_000, g(1, 2) as f64 / 1e9, g(1, 3), CHAIN_NS.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9, PART_NS[0].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9, PART_NS[1].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9, PART_NS[2].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9, PART_NS[3].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9, RESCUE[0].load(std::sync::atomic::Ordering::Relaxed), RESCUE[1].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9)
}

const NEG: i32 = i32::MIN / 4;

/// Alignment of q against r with two-piece affine gaps. Global (ext = false) aligns both completely;
/// an extension (ext = true) starts at the beginning of both and ends where the score is best,
/// with a bonus for reaching the end of q and a z-drop stop. Returns (score, q used, r used, ops).
type DpOut = (i32, usize, usize, Vec<(u8, u32)>);

/// results of large alignments by content: alleles of one gene share long stretches, so the same pair
/// of stretches comes up again and again when alleles are compared with each other
static DP_MEMO: std::sync::LazyLock<std::sync::Mutex<AHashMap<u128, DpOut>>> = std::sync::LazyLock::new(|| std::sync::Mutex::new(AHashMap::new()));

thread_local! {
    /// the same for the many small alignments of one thread: the end of a read is extended against
    /// the same stretch of every allele that agrees there
    static DP_SMALL: std::cell::RefCell<AHashMap<u128, DpOut>> = std::cell::RefCell::new(AHashMap::new());
}

pub fn dp(q: &[u8], r: &[u8], sc: &Score, ext: bool) -> DpOut {
    let cells = q.len() * r.len();
    if cells < 1000 {
        return dp_raw(q, r, sc, ext);
    }
    use std::hash::{Hash, Hasher};
    let h = |salt: u8| {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        (salt, q, r, ext).hash(&mut h);
        (sc.a, sc.b, sc.q, sc.e, sc.q2, sc.e2, sc.zdrop, sc.end_bonus, sc.bw).hash(&mut h);
        h.finish()
    };
    let key = ((h(1) as u128) << 64) | h(2) as u128;
    if cells < 200_000 {
        if let Some(v) = DP_SMALL.with(|m| m.borrow().get(&key).cloned()) {
            return v;
        }
        let v = dp_raw(q, r, sc, ext);
        DP_SMALL.with(|m| {
            let mut m = m.borrow_mut();
            if m.len() > 20_000 {
                m.clear();
            }
            m.insert(key, v.clone());
        });
        return v;
    }
    if let Some(v) = DP_MEMO.lock().unwrap().get(&key) {
        return v.clone();
    }
    let v = dp_raw(q, r, sc, ext);
    let mut m = DP_MEMO.lock().unwrap();
    if m.len() > 50_000 {
        m.clear();
    }
    m.insert(key, v.clone());
    v
}

fn dp_raw(q: &[u8], r: &[u8], sc: &Score, ext: bool) -> DpOut {
    let (n, m) = (q.len(), r.len());
    if n == 0 && m == 0 {
        return (0, 0, 0, vec![]);
    }
    if !ext && (n == 0 || m == 0) {
        let (op, l) = if n == 0 { (b'D', m) } else { (b'I', n) };
        let l32 = l as i32;
        return (-(sc.q + l32 * sc.e).min(sc.q2 + l32 * sc.e2), n, m, vec![(op, l as u32)]);
    }
    if ext && n == 0 {
        return (0, 0, 0, vec![]);
    }
    // equal lengths with few mismatches: any gapped path needs two gaps, which cost more
    if !ext && n == m {
        let mm = (0..n).filter(|&i| q[i] != r[i]).count() as i32;
        if mm * (sc.a + sc.b) <= 2 * (sc.q + sc.e).min(sc.q2 + sc.e2) {
            let mut ops = vec![];
            let mut s = 0;
            for i in 0..n {
                let eq = q[i] == r[i] && q[i] != b'N';
                s += if eq { sc.a } else if q[i] == b'N' || r[i] == b'N' { -1 } else { -sc.b };
                push(&mut ops, if eq { b'=' } else { b'X' }, 1);
            }
            return (s, n, n, ops);
        }
    }
    // band: global around the line from (0,0) to (n,m); extension around the main diagonal
    let w = if ext { sc.bw as i64 } else { (n as i64 - m as i64).abs() + 64 };
    let full = if ext { m <= sc.bw } else { (n as u64) * (m as u64) <= 100_000 };
    // extension: an ungapped run with at most two mismatches to the end of q is taken as it is when no
    // path with a gap can beat it: such a path scores at most n·a minus one gap of one base plus the end
    // bonus (near the end of q an insertion and the bonus can pay for themselves)
    if ext && m >= n && !q.contains(&b'N') && !r[..n].contains(&b'N') {
        let mm = (0..n).filter(|&i| q[i] != r[i]).take(3).count();
        if mm <= 2 {
            // best stopping point: before any mismatch, or the end with its bonus
            let (mut sc_i, mut best, mut bn, mut raw) = (0i32, 0i32, 0usize, 0i32);
            for i in 0..n {
                sc_i += if q[i] == r[i] { sc.a } else { -sc.b };
                raw = raw.max(sc_i);
                let v = if i + 1 == n { sc_i + sc.end_bonus } else { sc_i };
                if v > best || (i + 1 == n && v == best) {
                    best = v;
                    bn = i + 1;
                }
            }
            let mut ops = vec![];
            for i in 0..bn {
                push(&mut ops, if q[i] == r[i] { b'=' } else { b'X' }, 1);
            }
            if best >= n as i32 * sc.a - (sc.q + sc.e).min(sc.q2 + sc.e2) + sc.end_bonus {
                // the score is that of the best local end, the path may run on to the end of q for the bonus
                let s = (if bn == n { best - sc.end_bonus } else { best }).max(raw);
                return (s, bn, bn, ops);
            }
        }
    }
    if *STAT_ON {
        let t = std::time::Instant::now();
        let r = SCRATCH.with(|c| {
            let mut c = c.borrow_mut();
            let Scratch { hp, hc, f1, f2, tb, lo_v, off } = &mut *c;
            dp_in(q, r, sc, ext, n, m, full, w, hp, hc, f1, f2, tb, lo_v, off)
        });
        let k = if ext { if sc.bw > 100 { 2 } else { 0 } } else { 1 };
        STAT[k][0].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        STAT[k][1].fetch_add((n * m) as u64, std::sync::atomic::Ordering::Relaxed);
        STAT[k][2].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        if n * m > 1_000_000 { STAT[k][3].fetch_add(1, std::sync::atomic::Ordering::Relaxed); if std::env::var("HLA_ALNBIG").is_ok() { eprintln!("big dp ext={ext} n={n} m={m} full={full} w={w} us={} used={:?}", t.elapsed().as_micros(), (r.1, r.2)); } }
        return r;
    }
    SCRATCH.with(|c| {
        let mut c = c.borrow_mut();
        let Scratch { hp, hc, f1, f2, tb, lo_v, off } = &mut *c;
        dp_in(q, r, sc, ext, n, m, full, w, hp, hc, f1, f2, tb, lo_v, off)
    })
}

struct Scratch {
    hp: Vec<i32>,
    hc: Vec<i32>,
    f1: Vec<i32>,
    f2: Vec<i32>,
    tb: Vec<u8>,
    lo_v: Vec<usize>,
    off: Vec<usize>,
}

thread_local! {
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch { hp: vec![], hc: vec![], f1: vec![], f2: vec![], tb: vec![], lo_v: vec![], off: vec![] });
}

thread_local! {
    /// tests: the row-wise extension, to compare the anti-diagonal one with
    pub(crate) static FORCE_ROWS: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}
static EXT_AD: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var("HLA_EXT_AD").map_or(true, |v| v != "0"));

/// The extension of `dp_in` computed by anti-diagonals, so that the inner loop vectorises. The cells,
/// their values and traceback bits are those of the row-wise version; rows are summarised as they
/// complete, so the best end, the score and the z-drop stop are exactly those of the row-wise scan.
fn dp_ext_ad(q: &[u8], r: &[u8], sc: &Score, n: usize, m: usize, full: bool, w: i64, ext: bool) -> (i32, usize, usize, Vec<(u8, u32)>) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU has AVX2
        return unsafe { dp_ext_ad_avx2(q, r, sc, n, m, full, w, ext) };
    }
    dp_ext_ad_impl(q, r, sc, n, m, full, w, ext)
}

/// the same code compiled for AVX2 (x86 builds default to SSE2, which lacks a 32-bit maximum)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn dp_ext_ad_avx2(q: &[u8], r: &[u8], sc: &Score, n: usize, m: usize, full: bool, w: i64, ext: bool) -> (i32, usize, usize, Vec<(u8, u32)>) {
    dp_ext_ad_impl(q, r, sc, n, m, full, w, ext)
}

#[inline(always)]
fn dp_ext_ad_impl(q: &[u8], r: &[u8], sc: &Score, n: usize, m: usize, full: bool, w: i64, ext: bool) -> (i32, usize, usize, Vec<(u8, u32)>) {
    // rows that exist and their band, as in dp_in
    let nrow = if ext { n.min(m + w as usize) } else { n };
    let (rlo, rhi): (Vec<usize>, Vec<usize>) = (0..=nrow)
        .map(|i| {
            if full {
                return (0, m);
            }
            let c = if ext { i as i64 } else { (i as i64 * m as i64) / n.max(1) as i64 };
            ((c - w).max(0) as usize, ((c + w).min(m as i64)) as usize)
        })
        .unzip();
    let row_hi = |i: usize| rhi[i];
    // the rows of anti-diagonal t are those with i + lo(i) <= t <= i + hi(i); both bounds grow with t
    let (mut pa, mut pb) = (0usize, 0usize);
    let (oe1, oe2, e_1, e_2) = (sc.q + sc.e, sc.q2 + sc.e2, sc.e, sc.e2);
    let (sa, sb) = (sc.a, -sc.b);
    let rrev: Vec<u8> = r.iter().rev().copied().collect();
    // arrays indexed by i + 1 (index 0 is the neighbour of row 0)
    let sz = nrow + 3;
    let mut h0 = vec![NEG; sz];
    let mut h1 = vec![NEG; sz];
    let mut h2 = vec![NEG; sz];
    let mut ea = vec![NEG; sz];
    let mut eb = vec![NEG; sz];
    let mut fa_p = vec![NEG; sz];
    let mut fb_p = vec![NEG; sz];
    let mut fa_c = vec![NEG; sz];
    let mut fb_c = vec![NEG; sz];
    // per row: best h and its first column
    let mut rmax = vec![NEG; nrow + 1];
    let mut rarg = vec![0usize; nrow + 1];
    let mut tb: Vec<u8> = Vec::with_capacity((nrow + 1) * (2 * w.max(0) as usize + 2).min(m + 1));
    let mut toff: Vec<usize> = vec![];
    let mut tlo: Vec<usize> = vec![];
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    let mut rbest = 0i32;
    let mut next_row = 0usize; // first row not yet summarised
    let mut stop = false;
    let tmax = nrow + m;
    for t in 0..=tmax {
        while pa <= nrow && pa + rhi[pa] < t {
            pa += 1;
        }
        while pb < nrow && pb + 1 + rlo[pb + 1] <= t {
            pb += 1;
        }
        let (ilo, ihi) = (pa, pb);
        toff.push(tb.len());
        tlo.push(ilo);
        if ilo <= ihi {
            let base = tb.len();
            tb.resize(base + ihi - ilo + 1, 0);
            // cells with i >= 1 and j >= 1: one vectorisable loop
            let vlo = ilo.max(1);
            let vhi = ihi.min(t.saturating_sub(1));
            if vlo <= vhi && t >= 2 {
                let len = vhi - vlo + 1;
                let k0 = m + vlo - t; // r[t - i - 1] = rrev[m - t + i]
                let qs = &q[vlo - 1..vhi];
                let rs = &rrev[k0..k0 + len];
                // index i + 1 in the arrays; neighbours: up (i - 1) and diagonal (i - 1) at index i
                let hd = &h0[vlo..vlo + len];
                let hu = &h1[vlo..vlo + len];
                let hl = &h1[vlo + 1..vlo + 1 + len];
                let fap = &fa_p[vlo..vlo + len];
                let fbp = &fb_p[vlo..vlo + len];
                let ear = &mut ea[vlo + 1..vlo + 1 + len];
                let ebr = &mut eb[vlo + 1..vlo + 1 + len];
                let fac = &mut fa_c[vlo + 1..vlo + 1 + len];
                let fbc = &mut fb_c[vlo + 1..vlo + 1 + len];
                let ho = &mut h2[vlo + 1..vlo + 1 + len];
                let tbo = &mut tb[base + vlo - ilo..base + vlo - ilo + len];
                for x in 0..len {
                    let (o1, x1) = (hu[x] - oe1, fap[x] - e_1);
                    let (o2, x2) = (hu[x] - oe2, fbp[x] - e_2);
                    let mut tt = ((x1 > o1) as u8) << 5 | ((x2 > o2) as u8) << 6;
                    let (nf1, nf2) = (o1.max(x1), o2.max(x2));
                    let (p1, y1) = (hl[x] - oe1, ear[x] - e_1);
                    let (p2, y2) = (hl[x] - oe2, ebr[x] - e_2);
                    tt |= ((y1 > p1) as u8) << 3 | ((y2 > p2) as u8) << 4;
                    let (e1, e2) = (p1.max(y1), p2.max(y2));
                    let (a, b) = (qs[x], rs[x]);
                    let sv = if a == b'N' || b == b'N' { -1 } else if a == b { sa } else { sb };
                    let mut h = hd[x] + sv;
                    let mut src = 0u8;
                    if e1 > h { h = e1; src = 1; }
                    if e2 > h { h = e2; src = 2; }
                    if nf1 > h { h = nf1; src = 3; }
                    if nf2 > h { h = nf2; src = 4; }
                    ear[x] = e1;
                    ebr[x] = e2;
                    fac[x] = nf1;
                    fbc[x] = nf2;
                    ho[x] = h;
                    tbo[x] = tt | src;
                }
                let rm = &mut rmax[vlo..=vhi];
                let ra = &mut rarg[vlo..=vhi];
                for x in 0..len {
                    if ho[x] > rm[x] {
                        rm[x] = ho[x];
                        ra[x] = t - vlo - x;
                    }
                }
            }
            // row 0 (i = 0, j = t) and column 0 (j = 0, i = t), as in the row-wise code
            if ilo == 0 {
                let j = t;
                let (mut tt, h, e1, e2);
                if j == 0 {
                    h = 0;
                    tt = 7u8;
                    e1 = NEG;
                    e2 = NEG;
                } else {
                    let hl = h1[1];
                    let (o1, x1) = (hl - oe1, ea[1] - e_1);
                    let (o2, x2) = (hl - oe2, eb[1] - e_2);
                    tt = ((x1 > o1) as u8) << 3 | ((x2 > o2) as u8) << 4;
                    e1 = o1.max(x1);
                    e2 = o2.max(x2);
                    let mut hh = NEG;
                    let mut src = 0u8;
                    if e1 > hh { hh = e1; src = 1; }
                    if e2 > hh { hh = e2; src = 2; }
                    tt |= src;
                    h = hh;
                }
                ea[1] = e1;
                eb[1] = e2;
                fa_c[1] = NEG;
                fb_c[1] = NEG;
                h2[1] = h;
                tb[base] = tt;
                if h > rmax[0] {
                    rmax[0] = h;
                    rarg[0] = j;
                }
            }
            if t >= 1 && ihi == t && t <= nrow {
                let i = t;
                let hprev = h1[i];
                let (o1, x1) = (hprev - oe1, fa_p[i] - e_1);
                let (o2, x2) = (hprev - oe2, fb_p[i] - e_2);
                let mut tt = ((x1 > o1) as u8) << 5 | ((x2 > o2) as u8) << 6;
                let (nf1, nf2) = (o1.max(x1), o2.max(x2));
                let mut h = NEG;
                let mut src = 0u8;
                // horizontal gaps do not exist in column 0; the diagonal neither
                if NEG > h { h = NEG; }
                if nf1 > h { h = nf1; src = 3; }
                if nf2 > h { h = nf2; src = 4; }
                tt |= src;
                ea[i + 1] = NEG;
                eb[i + 1] = NEG;
                fa_c[i + 1] = nf1;
                fb_c[i + 1] = nf2;
                h2[i + 1] = h;
                tb[base + i - ilo] = tt;
                if h > rmax[i] {
                    rmax[i] = h;
                    rarg[i] = 0;
                }
            }
            // the cell before this anti-diagonal's start is read by the next two as a neighbour outside the band
            h2[ilo] = NEG;
            fa_c[ilo] = NEG;
            fb_c[ilo] = NEG;
        }
        std::mem::swap(&mut h0, &mut h1);
        std::mem::swap(&mut h1, &mut h2);
        std::mem::swap(&mut fa_p, &mut fa_c);
        std::mem::swap(&mut fb_p, &mut fb_c);
        if !ext {
            continue;
        }
        // rows whose last cell lies on this anti-diagonal are complete: the row-wise scan's bookkeeping
        while next_row <= nrow && next_row + row_hi(next_row) <= t {
            let i = next_row;
            let hb = if i == n { rmax[i] + sc.end_bonus } else { rmax[i] };
            if hb > best || (hb == best && i == n && bi != n) {
                best = hb;
                bi = i;
                bj = rarg[i];
            }
            rbest = rbest.max(rmax[i]);
            next_row += 1;
            if rmax[i] < best - sc.zdrop {
                stop = true;
                break;
            }
        }
        if stop || next_row > nrow {
            break;
        }
    }
    let score = if ext {
        (if bi == n { best - sc.end_bonus } else { best }).max(rbest)
    } else {
        // the last anti-diagonal holds only (n, m); its H ended in h1 after the swaps
        bi = n;
        bj = m;
        h1[n + 1]
    };
    // traceback over the anti-diagonal layout
    let mut ops: Vec<(u8, u32)> = vec![];
    let (mut i, mut j) = (bi, bj);
    let mut state = 0u8;
    while i > 0 || j > 0 {
        let tt = tb[toff[i + j] + i - tlo[i + j]];
        match state {
            0 => {
                let src = tt & 7;
                if src == 0 {
                    let op = if q[i - 1] == r[j - 1] { b'=' } else { b'X' };
                    ops.push((op, 1));
                    i -= 1;
                    j -= 1;
                } else if src == 7 {
                    break;
                } else {
                    state = src;
                }
            }
            1 | 2 => {
                ops.push((b'D', 1));
                let extd = tt & (1 << (state + 2)) != 0;
                j -= 1;
                if !extd { state = 0; }
            }
            _ => {
                ops.push((b'I', 1));
                let extd = tt & (1 << (state + 2)) != 0;
                i -= 1;
                if !extd { state = 0; }
            }
        }
    }
    ops.reverse();
    let mut out: Vec<(u8, u32)> = vec![];
    for (op, l) in ops {
        push(&mut out, op, l);
    }
    (score, bi, bj, out)
}

#[allow(clippy::too_many_arguments)]
fn dp_in(q: &[u8], r: &[u8], sc: &Score, ext: bool, n: usize, m: usize, full: bool, w: i64, hp: &mut Vec<i32>, hc: &mut Vec<i32>, f1: &mut Vec<i32>, f2: &mut Vec<i32>, tb: &mut Vec<u8>, lo_v: &mut Vec<usize>, off: &mut Vec<usize>) -> (i32, usize, usize, Vec<(u8, u32)>) {
    if *EXT_AD && !FORCE_ROWS.with(|c| c.get()) {
        return dp_ext_ad(q, r, sc, n, m, full, w, ext);
    }
    let band = |i: usize| -> (usize, usize) {
        if full {
            return (0, m);
        }
        let c = if ext { i as i64 } else { (i as i64 * m as i64) / n.max(1) as i64 };
        ((c - w).max(0) as usize, ((c + w).min(m as i64)) as usize)
    };
    let s = |a: u8, b: u8| -> i32 {
        if a == b'N' || b == b'N' {
            -1
        } else if a == b {
            sc.a
        } else {
            -sc.b
        }
    };
    let (oe1, oe2) = (sc.q + sc.e, sc.q2 + sc.e2);
    for v in [&mut *hp, &mut *hc, &mut *f1, &mut *f2] {
        v.clear();
        v.resize(m + 1, NEG);
    }
    tb.clear();
    lo_v.clear();
    off.clear();
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    let mut rbest = 0i32;
    for i in 0..=n {
        let (lo, hi) = band(i);
        // an extension does not go deeper than the band into q beyond the end of r
        if ext && i > m + w as usize {
            break;
        }
        let base = tb.len();
        tb.resize(base + hi - lo + 1, 0);
        let (mut e1, mut e2) = (NEG, NEG);
        let mut rmax = NEG;
        // best cell of this row, first one when several tie
        let mut rarg = lo;
        // the first cell of the row, and every cell of row 0, with all the edge cases
        let first_hi = if i == 0 { hi } else { lo };
        for j in lo..=first_hi {
            let mut t = 0u8;
            // the previous row's band starts at or before lo, and its cells beyond its end were never
            // written (NEG); the one cell before its start is set to NEG after that row
            let hprev_j = hp[j];
            // vertical gaps (insertions): consume q
            let (nf1, nf2) = if i > 0 {
                let (fp1, fp2) = (f1[j], f2[j]);
                let (o1, x1) = (hprev_j - oe1, fp1 - sc.e);
                let (o2, x2) = (hprev_j - oe2, fp2 - sc.e2);
                if x1 > o1 { t |= 1 << 5; }
                if x2 > o2 { t |= 1 << 6; }
                (o1.max(x1), o2.max(x2))
            } else {
                (NEG, NEG)
            };
            // horizontal gaps (deletions): consume r
            if j > lo || (j > 0 && j == lo && lo > 0) {
                let hl = if j > lo { hc[j - 1] } else { NEG };
                let (o1, x1) = (hl - oe1, e1 - sc.e);
                let (o2, x2) = (hl - oe2, e2 - sc.e2);
                if x1 > o1 { t |= 1 << 3; }
                if x2 > o2 { t |= 1 << 4; }
                e1 = o1.max(x1);
                e2 = o2.max(x2);
            } else {
                e1 = NEG;
                e2 = NEG;
            }
            let mut h;
            if i == 0 && j == 0 {
                h = 0;
                t |= 7;
            } else {
                let d = if i > 0 && j > 0 { hp[j - 1] + s(q[i - 1], r[j - 1]) } else { NEG };
                h = d;
                let mut src = 0u8;
                if e1 > h { h = e1; src = 1; }
                if e2 > h { h = e2; src = 2; }
                if nf1 > h { h = nf1; src = 3; }
                if nf2 > h { h = nf2; src = 4; }
                t |= src;
            }
            hc[j] = h;
            f1[j] = nf1;
            f2[j] = nf2;
            tb[base + j - lo] = t;
            if h > rmax {
                rmax = h;
                rarg = j;
            }
        }
        // the rest of the row (i > 0, lo < j <= hi): the same recurrences without edge cases, on
        // slices of the band so that the loop runs without bounds checks
        if i > 0 && hi > lo {
            let qi = q[i - 1];
            let qn = qi == b'N';
            let (a, b) = (sc.a, sc.b);
            let (e_1, e_2) = (sc.e, sc.e2);
            let hpr = &hp[lo..=hi];
            let hcr = &mut hc[lo..=hi];
            let f1r = &mut f1[lo..=hi];
            let f2r = &mut f2[lo..=hi];
            let rr = &r[lo..hi];
            let tbr = &mut tb[base..base + hi - lo + 1];
            let mut hleft = hcr[0];
            for k in 1..hpr.len() {
                let hprev = hpr[k];
                let (o1, x1) = (hprev - oe1, f1r[k] - e_1);
                let (o2, x2) = (hprev - oe2, f2r[k] - e_2);
                let mut t = ((x1 > o1) as u8) << 5 | ((x2 > o2) as u8) << 6;
                let (nf1, nf2) = (o1.max(x1), o2.max(x2));
                let (p1, y1) = (hleft - oe1, e1 - e_1);
                let (p2, y2) = (hleft - oe2, e2 - e_2);
                t |= ((y1 > p1) as u8) << 3 | ((y2 > p2) as u8) << 4;
                e1 = p1.max(y1);
                e2 = p2.max(y2);
                let rj = rr[k - 1];
                let sv = if qn || rj == b'N' { -1 } else if rj == qi { a } else { -b };
                let mut h = hpr[k - 1] + sv;
                let mut src = 0u8;
                if e1 > h { h = e1; src = 1; }
                if e2 > h { h = e2; src = 2; }
                if nf1 > h { h = nf1; src = 3; }
                if nf2 > h { h = nf2; src = 4; }
                hcr[k] = h;
                f1r[k] = nf1;
                f2r[k] = nf2;
                tbr[k] = t | src;
                hleft = h;
                if h > rmax {
                    rmax = h;
                    rarg = lo + k;
                }
            }
        }
        if ext {
            // as a scan of the row cell by cell: a better cell replaces the best; on the last row of q
            // (with its bonus) the first cell that reaches the best takes it
            let hb = if i == n { rmax + sc.end_bonus } else { rmax };
            if hb > best || (hb == best && i == n && bi != n) {
                best = hb;
                bi = i;
                bj = rarg;
            }
            rbest = rbest.max(rmax);
        }
        // bands only move right: the next row reads this row from lo - 1 on, so that cell must not be
        // a stale value of an older row in this buffer
        if lo > 0 {
            hc[lo - 1] = NEG;
        }
        off.push(base);
        lo_v.push(lo);
        std::mem::swap(hp, hc);
        if ext && rmax < best - sc.zdrop {
            break;
        }
    }
    let (ei, ej) = if ext { (bi, bj) } else { (n, m) };
    let score = if ext {
        (if bi == n { best - sc.end_bonus } else { best }).max(rbest)
    } else {
        hp[m]
    };
    // traceback
    let mut ops: Vec<(u8, u32)> = vec![];
    let (mut i, mut j) = (ei, ej);
    let mut state = 0u8; // 0 H, 1 E1, 2 E2, 3 F1, 4 F2
    while i > 0 || j > 0 {
        let t = tb[off[i] + j - lo_v[i]];
        match state {
            0 => {
                let src = t & 7;
                if src == 0 {
                    let op = if q[i - 1] == r[j - 1] { b'=' } else { b'X' };
                    ops.push((op, 1));
                    i -= 1;
                    j -= 1;
                } else if src == 7 {
                    break;
                } else {
                    state = src;
                }
            }
            1 | 2 => {
                ops.push((b'D', 1));
                let extd = t & (1 << (state + 2)) != 0;
                j -= 1;
                if !extd { state = 0; }
            }
            _ => {
                ops.push((b'I', 1));
                let extd = t & (1 << (state + 2)) != 0;
                i -= 1;
                if !extd { state = 0; }
            }
        }
    }
    ops.reverse();
    let mut out: Vec<(u8, u32)> = vec![];
    for (op, l) in ops {
        push(&mut out, op, l);
    }
    (score, ei, ej, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn global_and_extension() {
        let sc = Opts::sr().sc;
        let (s, _, _, o) = dp(b"ACGTACGT", b"ACGTACGT", &sc, false);
        assert_eq!((s, o), (16, vec![(b'=', 8)]));
        let (s, _, _, o) = dp(b"ACGTTACGT", b"ACGTACGT", &sc, false);
        assert_eq!(s, 16 - 14);
        assert_eq!(o.iter().filter(|x| x.0 == b'I').count(), 1);
        // an extension stops before a run of mismatches near the end unless the end bonus pays for it
        let (_, qn, _, _) = dp(b"ACGTACGTAA", b"ACGTACGTTT", &sc, true);
        assert_eq!(qn, 8);
    }

    #[test]
    fn map_read_both_strands() {
        let mut x = 12345u64;
        let r: Vec<u8> = (0..400)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 62) as usize]
            })
            .collect();
        let idx = Index::single(&r, Opts::sr_single());
        let mut q = r[100..250].to_vec();
        q[70] = if q[70] == b'A' { b'C' } else { b'A' };
        let h = &idx.map(&q)[0];
        assert!(h.fwd && h.r_st == 100 && h.r_en == 250 && h.nm == 1 && h.score == 149 * 2 - 8);
        let h = &idx.map(&revcomp(&q))[0];
        assert!(!h.fwd && h.r_st == 100 && h.q_st == 0 && h.q_en == 150);
        let mut d = r[100..170].to_vec();
        d.extend_from_slice(&r[173..253]);
        let h = &idx.map(&d)[0];
        assert_eq!(h.cs(&d, &r), format!(":70-{}:80", String::from_utf8_lossy(&r[170..173]).to_lowercase()));
    }
}

#[cfg(test)]
#[path = "aln_props.rs"]
mod props;
