// Rust implementations of star: STAR run as a program.

#[derive(Clone, Debug, Default)]
pub struct StarIndexOptsRs {
    pub threads: i64,
}

#[derive(Clone, Debug, Default)]
pub struct StarOptsRs {
    pub threads: i64,
    pub gtf: String,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct StarRunRs {
    pub genomeBam: String,
    pub transcriptomeBam: String,
    pub junctions: String,
    pub log: String,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct StarStatsRs {
    pub inputReads: i64,
    pub uniquelyMapped: i64,
    pub multiMapped: i64,
    pub tooManyLoci: i64,
    pub unmappedTooShort: i64,
    pub unmappedOther: i64,
    pub splices: i64,
    pub mismatchRate: f64,
}

pub fn star_sa_index_bases(len: i64) -> i64 {
    if len <= 0 {
        panic!("saIndexBases: genome length {} is not positive", len);
    }
    let x = (len as f64).log2() / 2.0 - 1.0;
    (x.round_ties_even() as i64).min(14)
}

fn star_threads(n: i64) -> String {
    let k = if n > 0 { n as usize } else { std::thread::available_parallelism().map(|k| k.get()).unwrap_or(1) };
    k.to_string()
}

fn star_genome_length(fasta: &str) -> i64 {
    use std::io::BufRead;
    let f = std::fs::File::open(fasta).unwrap_or_else(|e| panic!("buildStarIndex: cannot open {}: {}", fasta, e));
    let mut n = 0i64;
    for line in std::io::BufReader::new(f).lines() {
        let line = line.unwrap_or_else(|e| panic!("buildStarIndex: cannot read {}: {}", fasta, e));
        if !line.starts_with('>') {
            n += line.trim_end().len() as i64;
        }
    }
    n
}

fn star_run(args: &[String], what: &str, log: &str) {
    let exe = "STAR";
    let out = std::process::Command::new(exe)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("{}: cannot run {}: {}", what, exe, e));
    if !out.status.success() {
        let tail = |b: &[u8]| {
            let s = String::from_utf8_lossy(b);
            let k = s.len().saturating_sub(2000);
            s[k..].to_string()
        };
        let log_tail = std::fs::read(log).map(|b| tail(&b)).unwrap_or_default();
        panic!("{}: {} failed ({}):\n{}{}\n{}", what, exe, out.status, tail(&out.stdout), tail(&out.stderr), log_tail);
    }
}

pub fn star_build_index(o: &StarIndexOptsRs, fasta: &String, gtf: &String, dir: &String) -> String {
    std::fs::create_dir_all(dir).unwrap_or_else(|e| panic!("buildStarIndex: cannot create {}: {}", dir, e));
    let n = star_sa_index_bases(star_genome_length(fasta));
    let args: Vec<String> = vec![
        "--runMode".into(),
        "genomeGenerate".into(),
        "--genomeDir".into(),
        dir.clone(),
        "--genomeFastaFiles".into(),
        fasta.clone(),
        "--sjdbGTFfile".into(),
        gtf.clone(),
        "--runThreadN".into(),
        star_threads(o.threads),
        "--genomeSAindexNbases".into(),
        n.to_string(),
        "--outFileNamePrefix".into(),
        format!("{}/", dir),
    ];
    star_run(&args, "buildStarIndex", &format!("{}/Log.out", dir));
    dir.clone()
}

fn star_files_command(files: &[String]) -> Option<&'static str> {
    let gz = files.iter().filter(|f| f.ends_with(".gz")).count();
    if gz != 0 && gz != files.len() {
        panic!("alignStar: some read files are gzipped and some are not");
    }
    if gz > 0 { Some("zcat") } else { None }
}

fn star_align(o: &StarOptsRs, index: &str, mates: &[&[String]], prefix: &str) -> StarRunRs {
    if let Some(d) = std::path::Path::new(prefix).parent() {
        if !d.as_os_str().is_empty() {
            std::fs::create_dir_all(d).unwrap_or_else(|e| panic!("alignStar: cannot create {}: {}", d.display(), e));
        }
    }
    let all: Vec<String> = mates.iter().flat_map(|m| m.iter().cloned()).collect();
    if mates.iter().any(|m| m.is_empty()) {
        panic!("alignStar: no read files");
    }
    if mates.len() == 2 && mates[0].len() != mates[1].len() {
        panic!("alignStar: {} first-mate files but {} second-mate files", mates[0].len(), mates[1].len());
    }
    let sample = std::path::Path::new(prefix).file_name().map(|s| s.to_string_lossy().trim_end_matches('.').to_string()).unwrap_or_default();
    let mut args: Vec<String> = vec!["--genomeDir".into(), index.to_string(), "--readFilesIn".into()];
    for m in mates {
        args.push(m.join(","));
    }
    for a in [
        "--runThreadN",
        &star_threads(o.threads),
        "--outFileNamePrefix",
        prefix,
        "--sjdbGTFfile",
        &o.gtf,
        "--quantMode",
        "TranscriptomeSAM",
        "--outSAMtype",
        "BAM",
        "Unsorted",
        "--outSAMattributes",
        "NH",
        "HI",
        "AS",
        "NM",
        "MD",
        "--twopassMode",
        "Basic",
        "--runRNGseed",
        "0",
        "--outFilterMultimapNmax",
        "20",
        "--alignSJDBoverhangMin",
        "1",
        "--outSAMstrandField",
        "intronMotif",
        "--quantTranscriptomeSAMoutput",
        "BanSingleEnd",
        "--outSAMattrRGline",
    ] {
        args.push(a.to_string());
    }
    args.push(format!("ID:{}", sample));
    args.push(format!("SM:{}", sample));
    if let Some(c) = star_files_command(&all) {
        args.push("--readFilesCommand".into());
        args.push(c.into());
    }
    star_run(&args, "alignStar", &format!("{}Log.out", prefix));
    StarRunRs {
        genomeBam: format!("{}Aligned.out.bam", prefix),
        transcriptomeBam: format!("{}Aligned.toTranscriptome.out.bam", prefix),
        junctions: format!("{}SJ.out.tab", prefix),
        log: format!("{}Log.final.out", prefix),
    }
}

pub fn star_align_pairs(o: &StarOptsRs, index: &String, r1: &Vec<String>, r2: &Vec<String>, prefix: &String) -> StarRunRs {
    star_align(o, index, &[r1.as_slice(), r2.as_slice()], prefix)
}

pub fn star_align_reads(o: &StarOptsRs, index: &String, r: &Vec<String>, prefix: &String) -> StarRunRs {
    star_align(o, index, &[r.as_slice()], prefix)
}

pub fn star_stats(text: &String) -> StarStatsRs {
    let mut kv = std::collections::HashMap::new();
    for line in text.lines() {
        if let Some((k, v)) = line.split_once('|') {
            kv.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    let int = |k: &str| -> i64 {
        kv.get(k).unwrap_or_else(|| panic!("starStats: no '{}' in Log.final.out", k)).parse().unwrap_or_else(|_| panic!("starStats: '{}' is not a number", k))
    };
    let pct = |k: &str| -> f64 {
        let v = kv.get(k).unwrap_or_else(|| panic!("starStats: no '{}' in Log.final.out", k));
        v.trim_end_matches('%').parse().unwrap_or_else(|_| panic!("starStats: '{}' is not a percentage", k))
    };
    StarStatsRs {
        inputReads: int("Number of input reads"),
        uniquelyMapped: int("Uniquely mapped reads number"),
        multiMapped: int("Number of reads mapped to multiple loci"),
        tooManyLoci: int("Number of reads mapped to too many loci"),
        unmappedTooShort: int("Number of reads unmapped: too short"),
        unmappedOther: int("Number of reads unmapped: other"),
        splices: int("Number of splices: Total"),
        mismatchRate: pct("Mismatch rate per base, %"),
    }
}
