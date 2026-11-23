//! Local de-novo assembly of a candidate window, after GATK HaplotypeCaller: a de Bruijn graph of
//! the reference window and the read k-mers (bases >= Q10), edges seen in fewer than 2 reads
//! removed, the best paths from the first to the last reference k-mer as haplotypes, each aligned
//! back to the reference (Smith-Waterman parameters of GATK) to give its variant events. It finds
//! complex alleles that the read aligner placed differently. Research use only.
use std::collections::{HashMap, HashSet};

pub struct Settings {
    pub kmers: &'static [usize],
    pub min_base_q: u8,
    pub min_support: u32,
    pub max_haplotypes: usize,
}

pub const DEFAULT: Settings = Settings { kmers: &[25, 35, 45, 55], min_base_q: 10, min_support: 2, max_haplotypes: 32 };

/// haplotype sequences of the window (the reference first), or only the reference when no k gives
/// an acyclic graph with unique reference k-mers
pub fn assemble<'a>(refw: &'a [u8], reads: &[(&'a [u8], &'a [u8])], st: &Settings) -> Vec<Vec<u8>> {
    for &k in st.kmers {
        if refw.len() < k + 2 {
            break;
        }
        if let Some(h) = assemble_k(refw, reads, k, st) {
            return h;
        }
    }
    vec![refw.to_vec()]
}

fn node<'a>(s: &'a [u8], id: &mut HashMap<&'a [u8], usize>, kmer: &mut Vec<&'a [u8]>) -> usize {
    *id.entry(s).or_insert_with(|| {
        kmer.push(s);
        kmer.len() - 1
    })
}

fn assemble_k<'a>(refw: &'a [u8], reads: &[(&'a [u8], &'a [u8])], k: usize, st: &Settings) -> Option<Vec<Vec<u8>>> {
    let mut id: HashMap<&'a [u8], usize> = HashMap::new();
    let mut kmer: Vec<&'a [u8]> = vec![];
    // reference path; its k-mers must be unique
    let mut ref_nodes = vec![];
    for i in 0..=refw.len() - k {
        let n = node(&refw[i..i + k], &mut id, &mut kmer);
        if ref_nodes.contains(&n) {
            return None;
        }
        ref_nodes.push(n);
    }
    let mut edges: HashMap<(usize, usize), u32> = HashMap::new();
    let ref_edges: HashSet<(usize, usize)> = ref_nodes.windows(2).map(|w| (w[0], w[1])).collect();
    for e in &ref_edges {
        edges.insert(*e, 0);
    }
    for (seq, qual) in reads {
        if seq.len() < k + 1 {
            continue;
        }
        let mut seen: HashSet<(usize, usize)> = HashSet::new();
        let (seq, qual): (&'a [u8], &'a [u8]) = (seq, qual);
        let mut prev: Option<usize> = None;
        let mut bad_until = 0usize; // first k-mer start that has no low-quality base or N
        for i in 0..seq.len() {
            if qual.get(i).is_some_and(|&q| q < st.min_base_q) || !matches!(seq[i], b'A' | b'C' | b'G' | b'T') {
                bad_until = i + 1;
            }
            if i + 1 < k {
                continue;
            }
            let s = i + 1 - k;
            if s < bad_until {
                prev = None;
                continue;
            }
            let n = node(&seq[s..=i], &mut id, &mut kmer);
            if let Some(p) = prev {
                if seen.insert((p, n)) {
                    *edges.entry((p, n)).or_insert(0) += 1;
                }
            }
            prev = Some(n);
        }
    }
    // prune weak non-reference edges
    edges.retain(|e, c| ref_edges.contains(e) || *c >= st.min_support);
    let nn = kmer.len();
    let mut out: Vec<Vec<(usize, u32)>> = vec![vec![]; nn];
    let mut inc: Vec<Vec<usize>> = vec![vec![]; nn];
    for (&(a, b), &c) in &edges {
        out[a].push((b, c));
        inc[b].push(a);
    }
    let (src, snk) = (ref_nodes[0], *ref_nodes.last().unwrap());
    // keep nodes reachable from the source and reaching the sink
    let reach = |start: usize, adj: &dyn Fn(usize) -> Vec<usize>| -> Vec<bool> {
        let mut seen = vec![false; nn];
        let mut st = vec![start];
        seen[start] = true;
        while let Some(x) = st.pop() {
            for y in adj(x) {
                if !seen[y] {
                    seen[y] = true;
                    st.push(y);
                }
            }
        }
        seen
    };
    let fwd = reach(src, &|x| out[x].iter().map(|e| e.0).collect());
    let bwd = reach(snk, &|x| inc[x].clone());
    let keep: Vec<bool> = (0..nn).map(|i| fwd[i] && bwd[i]).collect();
    for x in 0..nn {
        out[x].retain(|e| keep[e.0]);
    }
    // cycles among the kept nodes: try a larger k
    let mut state = vec![0u8; nn];
    fn cyclic(x: usize, out: &[Vec<(usize, u32)>], keep: &[bool], state: &mut [u8]) -> bool {
        state[x] = 1;
        for &(y, _) in &out[x] {
            if !keep[y] {
                continue;
            }
            if state[y] == 1 || (state[y] == 0 && cyclic(y, out, keep, state)) {
                return true;
            }
        }
        state[x] = 2;
        false
    }
    if cyclic(src, &out, &keep, &mut state) {
        return None;
    }
    // best paths: score = sum of log(edge count / total outgoing count), reference edges count 1
    let w = |c: u32| (c.max(1)) as f64;
    let mut paths: Vec<(f64, Vec<usize>)> = vec![];
    let mut stack: Vec<(usize, Vec<usize>, f64)> = vec![(src, vec![src], 0.0)];
    let mut steps = 0usize;
    while let Some((x, path, score)) = stack.pop() {
        steps += 1;
        if steps > 200_000 {
            break;
        }
        if x == snk {
            paths.push((score, path));
            continue;
        }
        let tot: f64 = out[x].iter().map(|e| w(e.1)).sum();
        for &(y, c) in &out[x] {
            let mut p = path.clone();
            p.push(y);
            stack.push((y, p, score + (w(c) / tot).ln()));
        }
    }
    paths.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let mut haps: Vec<Vec<u8>> = vec![refw.to_vec()];
    for (_, p) in paths {
        let mut s = kmer[p[0]].to_vec();
        for &n in &p[1..] {
            s.push(kmer[n][k - 1]);
        }
        if !haps.contains(&s) {
            haps.push(s);
        }
        if haps.len() > st.max_haplotypes {
            break;
        }
    }
    Some(haps)
}

/// events of a haplotype against the reference window: (offset in the window, REF, ALT) with an
/// anchor base for insertions and deletions; global alignment with affine gaps, GATK's
/// Smith-Waterman parameters (match 200, mismatch -150, gap open -260, extension -11)
pub fn events(refw: &[u8], hap: &[u8]) -> Vec<(usize, Vec<u8>, Vec<u8>)> {
    let (n, m) = (refw.len(), hap.len());
    const NEG: i32 = i32::MIN / 4;
    let (mat, mis, open, ext) = (200, -150, -260, -11);
    // states: 0 = M (match/mismatch), 1 = D (gap in the haplotype), 2 = I (gap in the reference)
    let w = m + 1;
    let idx = |i: usize, j: usize| i * w + j;
    let mut sc = [vec![NEG; (n + 1) * w], vec![NEG; (n + 1) * w], vec![NEG; (n + 1) * w]];
    let mut bp = [vec![0u8; (n + 1) * w], vec![0u8; (n + 1) * w], vec![0u8; (n + 1) * w]];
    sc[0][idx(0, 0)] = 0;
    let best3 = |v: [i32; 3]| -> (i32, u8) {
        let mut b = 0;
        for s in 1..3 {
            if v[s] > v[b] {
                b = s;
            }
        }
        (v[b], b as u8)
    };
    for i in 0..=n {
        for j in 0..=m {
            if i == 0 && j == 0 {
                continue;
            }
            if i > 0 && j > 0 {
                let s = if refw[i - 1] == hap[j - 1] { mat } else { mis };
                let (v, b) = best3([sc[0][idx(i - 1, j - 1)], sc[1][idx(i - 1, j - 1)], sc[2][idx(i - 1, j - 1)]]);
                sc[0][idx(i, j)] = v + s;
                bp[0][idx(i, j)] = b;
            }
            if i > 0 {
                let (v, b) = best3([sc[0][idx(i - 1, j)] + open, sc[1][idx(i - 1, j)] + ext, sc[2][idx(i - 1, j)] + open]);
                sc[1][idx(i, j)] = v;
                bp[1][idx(i, j)] = b;
            }
            if j > 0 {
                let (v, b) = best3([sc[0][idx(i, j - 1)] + open, sc[1][idx(i, j - 1)] + open, sc[2][idx(i, j - 1)] + ext]);
                sc[2][idx(i, j)] = v;
                bp[2][idx(i, j)] = b;
            }
        }
    }
    let (_, mut state) = best3([sc[0][idx(n, m)], sc[1][idx(n, m)], sc[2][idx(n, m)]]);
    let (mut i, mut j) = (n, m);
    let mut ops: Vec<u8> = vec![];
    while i > 0 || j > 0 {
        let prev = bp[state as usize][idx(i, j)];
        ops.push(state);
        match state {
            0 => {
                i -= 1;
                j -= 1;
            }
            1 => i -= 1,
            _ => j -= 1,
        }
        state = prev;
    }
    ops.reverse();
    let mut ev = vec![];
    let (mut ri, mut hj) = (0usize, 0usize);
    let mut k = 0;
    while k < ops.len() {
        match ops[k] {
            0 => {
                if refw[ri] != hap[hj] {
                    ev.push((ri, vec![refw[ri]], vec![hap[hj]]));
                }
                ri += 1;
                hj += 1;
                k += 1;
            }
            1 => {
                let mut l = 0;
                while k + l < ops.len() && ops[k + l] == 1 {
                    l += 1;
                }
                if ri > 0 {
                    ev.push((ri - 1, refw[ri - 1..ri + l].to_vec(), vec![refw[ri - 1]]));
                }
                ri += l;
                k += l;
            }
            _ => {
                let mut l = 0;
                while k + l < ops.len() && ops[k + l] == 2 {
                    l += 1;
                }
                if ri > 0 {
                    let mut alt = vec![refw[ri - 1]];
                    alt.extend(&hap[hj..hj + l]);
                    ev.push((ri - 1, vec![refw[ri - 1]], alt));
                }
                hj += l;
                k += l;
            }
        }
    }
    ev
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snv_and_deletion_events() {
        let r = b"ACGTACGTTAGCATCGATCGGATCCATGCAAGT";
        let mut h = r.to_vec();
        h[5] = b'T'; // C>T at offset 5
        h.drain(15..18); // delete 3 bases at offsets 15..18 (anchor 14)
        let e = events(r, &h);
        assert_eq!(e.len(), 2, "{e:?}");
        assert_eq!(e[0], (5, b"C".to_vec(), b"T".to_vec()));
        assert_eq!(e[1].0, 14);
        assert_eq!(e[1].1.len(), 4);
        assert_eq!(e[1].2.len(), 1);
    }

    #[test]
    fn assembles_a_read_supported_haplotype() {
        let r: Vec<u8> = b"TTGACCATGCAGTACGATCGATCGTAGCTAGCTAGGATCGATCGTTAGCATCGACTGACTAGCATCGATCAGCTAGGCATCGA".to_vec();
        let mut alt = r.clone();
        alt[40] = if r[40] == b'A' { b'C' } else { b'A' };
        let q = vec![30u8; alt.len()];
        let reads: Vec<(&[u8], &[u8])> = vec![(&alt[..], &q[..]), (&alt[..], &q[..]), (&alt[..], &q[..])];
        let h = assemble(&r, &reads, &DEFAULT);
        assert!(h.len() >= 2);
        assert!(h.iter().any(|x| x == &alt));
    }
}
