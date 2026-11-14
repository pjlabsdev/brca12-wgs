#!/usr/bin/env nextflow
/*
 * brca12-wgs · BRCA1/BRCA2 germline analysis from whole-genome sequencing (GRCh38/hg38 only)
 *
 * FASTQ (R1/R2) or an existing BAM/CRAM
 *   -> BWA-MEM -> samtools fixmate/sort/markdup/index
 *   -> GATK HaplotypeCaller (or DeepVariant)    SNVs/small indels
 *   -> Manta                                    structural variants
 *   -> CNVkit (WGS, flat reference)             deletions/duplications
 *   -> extraction of BRCA1/BRCA2 via brca12_regions.bed
 *   -> Ensembl VEP (+ ClinVar)                  annotation
 *   -> mosdepth                                 coverage per exon
 *   -> brca12_report.html / .md / brca12_results.json
 *
 * Research use only (RUO). Results must not be interpreted medically without clinically
 * validated confirmation (accredited laboratory, e.g. Sanger sequencing or MLPA).
 */

nextflow.enable.dsl = 2

params.sample          = 'sample'
params.reads_r1        = null      // sample_R1.fastq.gz
params.reads_r2        = null      // sample_R2.fastq.gz
params.bam             = null      // alternative to FASTQ: coordinate-sorted BAM/CRAM (GRCh38)
params.reference       = null      // GRCh38 FASTA with .fai, .dict and BWA index (.bwt/.pac/.ann/.amb/.sa)
params.regions_bed     = "${projectDir}/assets/brca12_regions.bed"
params.exons_bed       = "${projectDir}/assets/brca12_exons.bed"
params.vep_cache       = null      // Ensembl VEP cache directory (GRCh38)
params.vep_cache_version = null    // optional, e.g. 113
params.clinvar_vcf     = null      // clinvar.vcf.gz (+ .tbi), GRCh38
// params.outdir: see nextflow.config (default brca12_results)
// Optional evidence for the ENIGMA v1.2 hints. All are bgzip-compressed with a tabix index and
// are supplied by the user (SpliceAI and dbNSFP/BayesDel: non-commercial research licence).
params.gnomad_exomes_vcf  = null   // gnomAD v2.1 exomes (GRCh38 liftover) with non_cancer FAF fields
params.gnomad_genomes_vcf = null   // gnomAD v3.1 genomes with non_cancer FAF fields
params.spliceai_vcf       = null   // spliceai_scores.raw.snv.hg38.vcf.gz
params.bayesdel_tsv       = null   // chrom, pos, ref, alt, BayesDel_noAF score
params.threads         = 8
params.snv_caller      = 'haplotypecaller'   // or 'deepvariant'
params.calling_padding = 200       // bp around BRCA1/BRCA2 for variant calling
params.sv_padding      = 1000000   // bp search window for Manta breakpoints
params.min_exon_depth  = 20        // exon flagged uncertain if its mean depth is below this
params.min_exon_pct20  = 95        // or if fewer than N % of bases reach 20x
params.low_depth       = 10        // intervals below this depth are listed as uncertain
params.max_mapq0       = 0.10      // fraction of reads with MAPQ < 20 above which an exon counts as mapping-uncertain

// Executables can be overridden (own tool installations instead of PATH).
params.bwa       = 'bwa'
params.samtools  = 'samtools'
params.bcftools  = 'bcftools'
params.tabix     = 'tabix'
params.bgzip     = 'bgzip'
params.gatk      = 'gatk'
params.deepvariant = 'run_deepvariant'
params.manta_config = 'configManta.py'
params.cnvkit    = 'cnvkit.py'
params.vep       = 'vep'
params.mosdepth  = 'mosdepth'
params.python    = 'python3'
params.report_lang = 'en'      // language of brca12_report.html/.md: en, de or ru

def requireParam(name) {
    if (!params[name]) { error "Parameter --${name} is missing (see README.md)." }
}

process VALIDATE_INPUTS {
    tag "${params.sample}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path reference_fai
    path regions_bed
    path exons_bed

    output:
    path 'reference_check.json', emit: check
    path 'regions.calling.bed',  emit: calling_bed
    path 'regions.sv.bed',       emit: sv_bed
    path 'regions.exact.bed',    emit: exact_bed
    path 'exons.bed',            emit: exons
    path 'validate.command.sh'

    script:
    """
    ${params.python} ${projectDir}/bin/brca12_prepare.py \\
        --fai ${reference_fai} --regions ${regions_bed} --exons ${exons_bed} \\
        --calling-padding ${params.calling_padding} --sv-padding ${params.sv_padding} \\
        --out-check reference_check.json --out-calling regions.calling.bed \\
        --out-sv regions.sv.bed --out-exact regions.exact.bed --out-exons exons.bed
    cp .command.sh validate.command.sh
    """
}

process ALIGN_BWA_MEM {
    tag "${params.sample}"
    cpus params.threads
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '*.{bam,bai,txt,sh}'

    input:
    tuple path(r1), path(r2)
    path reference
    path reference_index_files

    output:
    tuple path("${params.sample}.markdup.bam"), path("${params.sample}.markdup.bam.bai"), emit: bam
    path "${params.sample}.markdup.stats.txt"
    path 'align.command.sh'

    script:
    def rg = "@RG\\tID:${params.sample}\\tSM:${params.sample}\\tLB:${params.sample}\\tPL:ILLUMINA"
    """
    set -euo pipefail
    ${params.bwa} mem -t ${task.cpus} -K 100000000 -Y -R '${rg}' ${reference} ${r1} ${r2} \\
      | ${params.samtools} fixmate -m -@ ${task.cpus} - - \\
      | ${params.samtools} sort -@ ${task.cpus} -m 1G -o ${params.sample}.sorted.bam -
    ${params.samtools} markdup -@ ${task.cpus} -s -f ${params.sample}.markdup.stats.txt \\
      ${params.sample}.sorted.bam ${params.sample}.markdup.bam
    ${params.samtools} index -@ ${task.cpus} ${params.sample}.markdup.bam
    rm -f ${params.sample}.sorted.bam
    cp .command.sh align.command.sh
    """
}

process CALL_SMALL_VARIANTS {
    tag "${params.sample}:${params.snv_caller}"
    cpus params.threads
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path calling_bed

    output:
    tuple path("${params.sample}.small.vcf.gz"), path("${params.sample}.small.vcf.gz.tbi"), emit: vcf
    path 'small_variants.command.sh'

    script:
    if (params.snv_caller == 'deepvariant')
        """
        set -euo pipefail
        ${params.deepvariant} --model_type=WGS --ref=${reference} --reads=${bam} \\
          --regions=${calling_bed} --output_vcf=${params.sample}.small.vcf.gz \\
          --num_shards=${task.cpus}
        ${params.tabix} -f -p vcf ${params.sample}.small.vcf.gz
        cp .command.sh small_variants.command.sh
        """
    else
        """
        set -euo pipefail
        ${params.gatk} --java-options "-Xmx4g" HaplotypeCaller \\
          -R ${reference} -I ${bam} -L ${calling_bed} \\
          --native-pair-hmm-threads ${task.cpus} \\
          -O ${params.sample}.small.vcf.gz
        cp .command.sh small_variants.command.sh
        """
}

process CALL_SV_MANTA {
    tag "${params.sample}"
    cpus params.threads
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path sv_bed

    output:
    tuple path("${params.sample}.manta.diploidSV.vcf.gz"), path("${params.sample}.manta.diploidSV.vcf.gz.tbi"), emit: vcf
    path 'manta.command.sh'

    script:
    """
    set -euo pipefail
    sort -k1,1 -k2,2n ${sv_bed} | ${params.bgzip} -c > call_regions.bed.gz
    ${params.tabix} -f -p bed call_regions.bed.gz
    ${params.manta_config} --bam ${bam} --referenceFasta ${reference} \\
      --callRegions call_regions.bed.gz --runDir manta_run
    manta_run/runWorkflow.py -j ${task.cpus}
    cp manta_run/results/variants/diploidSV.vcf.gz ${params.sample}.manta.diploidSV.vcf.gz
    cp manta_run/results/variants/diploidSV.vcf.gz.tbi ${params.sample}.manta.diploidSV.vcf.gz.tbi
    cp .command.sh manta.command.sh
    """
}

process CALL_CNV_CNVKIT {
    tag "${params.sample}"
    cpus params.threads
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path reference_check

    output:
    path "cnvkit/${params.sample}.call.cns", emit: cns
    path "cnvkit/${params.sample}.cnr", emit: cnr
    path 'cnv.command.sh'

    script:
    // CNVkit WGS mode with a flat reference on chr13 and chr17 (normalised across both
    // chromosomes). A single sample without a cohort yields candidates only; GATK gCNV
    // would need a trained cohort model. --drop-low-coverage keeps bins with (almost) no reads
    // out of the segmentation; without it such coverage holes dominate the CBS statistics and
    // heterozygous deletions up to whole-gene size were missed in the spike-in validation.
    """
    set -euo pipefail
    ${params.cnvkit} access ${reference} -o access.bed
    ${params.python} ${projectDir}/bin/brca12_prepare.py --cnv-targets ${reference_check} \\
      --cnv-access access.bed --out-cnv-targets wgs_regions.bed
    # In WGS mode CNVkit requires the same BED for --targets and --access.
    ${params.cnvkit} batch ${bam} --method wgs --normal --fasta ${reference} \\
      --targets wgs_regions.bed --access wgs_regions.bed --target-avg-size 1000 \\
      --drop-low-coverage \\
      --output-reference flat_reference.cnn --output-dir cnvkit -p ${task.cpus}
    # CNVkit names its output after the BAM; rename unless the BAM is already called like the sample.
    if [ "${bam.baseName}" != "${params.sample}" ]; then
      mv cnvkit/${bam.baseName}.call.cns cnvkit/${params.sample}.call.cns
      mv cnvkit/${bam.baseName}.cnr cnvkit/${params.sample}.cnr
    fi
    cp .command.sh cnv.command.sh
    """
}

process EXON_COVERAGE {
    tag "${params.sample}"
    cpus 4
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path exons
    path exact_bed

    output:
    path 'exons_q20.regions.bed.gz', emit: q20
    path 'exons_q20.thresholds.bed.gz', emit: thresholds
    path 'exons_q0.regions.bed.gz', emit: q0
    path 'gene_q20.per-base.bed.gz', emit: perbase
    path 'chrom_median.txt', emit: chrom_median
    path 'coverage.command.sh'
    path 'insertion_signals.json', emit: insertion

    script:
    """
    set -euo pipefail
    ${params.mosdepth} -t ${task.cpus} -f ${reference} -Q 20 -n --by ${exons} --thresholds 10,20,30 exons_q20 ${bam}
    ${params.mosdepth} -t ${task.cpus} -f ${reference} -Q 0  -n --by ${exons} exons_q0 ${bam}
    # Per-base depth inside BRCA1/BRCA2 only (MAPQ >= 20, no duplicates) for uncertain intervals.
    ${params.samtools} depth -a -Q 20 -b ${exact_bed} --reference ${reference} ${bam} \\
      | awk 'BEGIN{OFS="\\t"} {print \$1, \$2-1, \$2, \$3}' | gzip -c > gene_q20.per-base.bed.gz
    # Chromosome median depth (chr13/chr17) as baseline for per-exon dosage.
    for C in \$(cut -f1 ${exact_bed} | sort -u); do
      ${params.mosdepth} -t ${task.cpus} -f ${reference} -Q 20 -n --fast-mode --by 100000 -c \$C chrom_\$C ${bam}
    done
    ${params.python} ${projectDir}/bin/brca12_prepare.py --chrom-median chrom_*.regions.bed.gz --out-chrom-median chrom_median.txt
    # Insertion signatures (e.g. Alu) from soft-clipped read ends inside the genes.
    ${params.samtools} view -L ${exact_bed} --reference ${reference} ${bam} \\
      | ${params.python} ${projectDir}/bin/brca12_insertion_signals.py --regions ${exact_bed} --out insertion_signals.json
    cp .command.sh coverage.command.sh
    """
}

process EXTRACT_AND_ANNOTATE {
    tag "${params.sample}"
    cpus 4
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    path exact_bed
    path reference
    path reference_index_files
    path vep_cache
    tuple path(clinvar), path(clinvar_tbi)

    output:
    path 'BRCA1.snv_indel.vep.vcf.gz', emit: brca1
    path 'BRCA2.snv_indel.vep.vcf.gz', emit: brca2
    path 'annotate.command.sh'

    script:
    def cacheVersion = params.vep_cache_version ? "--cache_version ${params.vep_cache_version}" : ''
    """
    set -euo pipefail
    ${params.bcftools} norm -m -any -f ${reference} ${vcf} -Oz -o norm.vcf.gz
    ${params.tabix} -f -p vcf norm.vcf.gz
    for GENE in BRCA1 BRCA2; do
      awk -v g=\$GENE 'BEGIN{OFS="\\t"} \$4==g {print \$1,\$2,\$3,\$4}' ${exact_bed} > \$GENE.bed
      ${params.bcftools} view -R \$GENE.bed norm.vcf.gz -Oz -o \$GENE.snv_indel.vcf.gz
      ${params.tabix} -f -p vcf \$GENE.snv_indel.vcf.gz
      ${params.vep} --offline --cache --dir_cache ${vep_cache} ${cacheVersion} --assembly GRCh38 \\
        --fasta ${reference} --format vcf --vcf --compress_output bgzip --force_overwrite \\
        --symbol --hgvs --canonical --mane --variant_class --numbers --biotype \\
        --af_gnomade --af_gnomadg --max_af --sift b --polyphen b --flag_pick \\
        --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \\
        --custom file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \\
        --fork ${task.cpus} \\
        -i \$GENE.snv_indel.vcf.gz -o \$GENE.snv_indel.vep.vcf.gz
    done
    cp .command.sh annotate.command.sh
    """
}

process TOOL_VERSIONS {
    output:
    path 'tool_versions.tsv', emit: versions

    script:
    """
    {
      echo -e "nextflow\\t${workflow.nextflow.version}"
      echo -e "bwa\\t\$( (${params.bwa} 2>&1 || true) | awk '/^Version/{print \$2}')"
      echo -e "samtools\\t\$(${params.samtools} --version 2>/dev/null | head -1 | awk '{print \$2}')"
      echo -e "bcftools\\t\$(${params.bcftools} --version 2>/dev/null | head -1 | awk '{print \$2}')"
      if [ "${params.snv_caller}" = "deepvariant" ]; then
        echo -e "deepvariant\\t\$(${params.deepvariant} --version 2>&1 | tail -1 | awk '{print \$NF}')"
      else
        echo -e "gatk\\t\$(${params.gatk} --version 2>/dev/null | awk '/The Genome Analysis Toolkit/{print \$NF}')"
      fi
      echo -e "manta\\t\$(${params.manta_config} --version 2>/dev/null | tail -1)"
      echo -e "cnvkit\\t\$(${params.cnvkit} version 2>/dev/null | tail -1)"
      echo -e "ensembl-vep\\t\$(${params.vep} --help 2>/dev/null | awk '/ensembl-vep/{print \$NF; exit}')"
      echo -e "mosdepth\\t\$(${params.mosdepth} --version 2>/dev/null | awk '{print \$2}')"
      echo -e "python\\t\$(${params.python} --version 2>&1 | awk '{print \$2}')"
    } > tool_versions.tsv
    """
}

process REPORT {
    tag "${params.sample}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path reference_check
    path exons
    path exact_bed
    path brca1_vcf
    path brca2_vcf
    path manta_vcf
    path cnv_cns
    path cnv_cnr
    path cov_q20
    path cov_thresholds
    path cov_q0
    path cov_perbase
    path chrom_median
    path insertion_signals
    path versions
    path command_files
    path evidence_files

    output:
    path 'brca12_report.html'
    path 'brca12_report.md'
    path 'brca12_results.json'
    path 'brca12_commands.sh'

    script:
    def input_mode = params.bam ? 'bam' : 'fastq'
    def inputs = params.bam ? params.bam : "${params.reads_r1},${params.reads_r2}"
    def staged = { p -> p ? file(p).name : null }
    def evidence = [['--gnomad-exomes', params.gnomad_exomes_vcf], ['--gnomad-genomes', params.gnomad_genomes_vcf],
                    ['--spliceai', params.spliceai_vcf], ['--bayesdel', params.bayesdel_tsv]]
        .findAll { it[1] }.collect { "${it[0]} ${staged(it[1])}" }.join(' ')
    """
    cat ${command_files} > brca12_commands.sh
    ${params.python} ${projectDir}/bin/brca12_report.py \\
      --sample ${params.sample} --input-mode ${input_mode} --inputs '${inputs}' \\
      --reference ${params.reference} --reference-check ${reference_check} \\
      --exons ${exons} --regions ${exact_bed} \\
      --small BRCA1=${brca1_vcf} --small BRCA2=${brca2_vcf} \\
      --sv ${manta_vcf} --cnv ${cnv_cns} --cnr ${cnv_cnr} \\
      --cov-q20 ${cov_q20} --cov-thresholds ${cov_thresholds} --cov-q0 ${cov_q0} \\
      --cov-perbase ${cov_perbase} --chrom-median ${chrom_median} \\
      --insertion-signals ${insertion_signals} \\
      --versions ${versions} --commands brca12_commands.sh \\
      --snv-caller ${params.snv_caller} \\
      --min-exon-depth ${params.min_exon_depth} --min-exon-pct20 ${params.min_exon_pct20} \\
      --low-depth ${params.low_depth} --max-mapq0 ${params.max_mapq0} \\
      --nextflow-version ${workflow.nextflow.version} --session-id ${workflow.sessionId} \\
      --lang ${params.report_lang} --tabix ${params.tabix} ${evidence} \\
      --out-json brca12_results.json --out-md brca12_report.md --out-html brca12_report.html
    """
}

process REPORT_PDF {
    // Optional PDF of the HTML report (--report_pdf), printed with headless Chrome/Chromium
    // (--chrome); skipped with a note when none is available.
    tag "${params.sample}"
    publishDir "${params.outdir}", mode: 'copy'

    when:
    params.report_pdf

    input:
    path report_html

    output:
    path 'brca12_report.pdf', optional: true

    script:
    """
    bash ${projectDir}/bin/brca12_pdf.sh ${report_html} brca12_report.pdf '${params.chrome ?: ''}'
    """
}

workflow {
    requireParam('reference')
    requireParam('vep_cache')
    requireParam('clinvar_vcf')
    if (!params.bam && !(params.reads_r1 && params.reads_r2)) {
        error 'Provide either --reads_r1/--reads_r2 (FASTQ) or --bam.'
    }

    reference = file(params.reference, checkIfExists: true)
    fai = file("${params.reference}.fai", checkIfExists: true)
    ref_files = Channel.fromPath([
        "${params.reference}.fai",
        "${params.reference}.{amb,ann,bwt,pac,sa}",
        "${params.reference.replaceAll(/\.(fa|fasta)(\.gz)?$/, '')}.dict"
    ].flatten(), checkIfExists: false).collect()

    VALIDATE_INPUTS(fai, file(params.regions_bed, checkIfExists: true), file(params.exons_bed, checkIfExists: true))

    if (params.bam) {
        def bai = file("${params.bam}.bai").exists() ? file("${params.bam}.bai")
                : file("${params.bam}".replaceAll(/\.bam$/, '.bai')).exists() ? file("${params.bam}".replaceAll(/\.bam$/, '.bai'))
                : file("${params.bam}.crai", checkIfExists: true)
        bam_ch = Channel.value(tuple(file(params.bam, checkIfExists: true), bai))
    } else {
        reads = Channel.value(tuple(file(params.reads_r1, checkIfExists: true), file(params.reads_r2, checkIfExists: true)))
        bam_ch = ALIGN_BWA_MEM(reads, reference, ref_files).bam
    }

    small = CALL_SMALL_VARIANTS(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.calling_bed)
    sv = CALL_SV_MANTA(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.sv_bed)
    cnv = CALL_CNV_CNVKIT(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.check)
    cov = EXON_COVERAGE(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.exons, VALIDATE_INPUTS.out.exact_bed)
    clinvar = Channel.value(tuple(file(params.clinvar_vcf, checkIfExists: true), file("${params.clinvar_vcf}.tbi", checkIfExists: true)))
    ann = EXTRACT_AND_ANNOTATE(small.vcf, VALIDATE_INPUTS.out.exact_bed, reference, ref_files,
                               file(params.vep_cache, checkIfExists: true), clinvar)
    versions = TOOL_VERSIONS()

    commands = Channel.empty()
        .mix(VALIDATE_INPUTS.out[5])
        .mix(params.bam ? Channel.empty() : ALIGN_BWA_MEM.out[2])
        .mix(CALL_SMALL_VARIANTS.out[1], CALL_SV_MANTA.out[1], CALL_CNV_CNVKIT.out[2],
             EXON_COVERAGE.out[5], EXTRACT_AND_ANNOTATE.out[2])
        .collect()

    // Optional evidence files and their tabix indexes; an empty placeholder keeps the input non-empty.
    evidence_files = Channel.fromPath(
        [params.gnomad_exomes_vcf, params.gnomad_genomes_vcf, params.spliceai_vcf, params.bayesdel_tsv]
            .findAll { it }.collectMany { [it, "${it}.tbi"] } + ["${projectDir}/assets/NO_EVIDENCE"],
        checkIfExists: true).collect()

    REPORT(VALIDATE_INPUTS.out.check, VALIDATE_INPUTS.out.exons, VALIDATE_INPUTS.out.exact_bed,
           ann.brca1, ann.brca2, sv.vcf.map { it[0] }, cnv.cns, cnv.cnr,
           cov.q20, cov.thresholds, cov.q0, cov.perbase, cov.chrom_median, cov.insertion,
           versions.versions, commands, evidence_files)
    REPORT_PDF(REPORT.out[0])
}
