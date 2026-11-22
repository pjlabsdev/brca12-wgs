#!/usr/bin/env nextflow
/*
 * brca12-wgs · BRCA1/BRCA2 germline analysis from whole-genome sequencing (GRCh38/hg38 only)
 *
 * FASTQ (R1/R2) or an existing BAM/CRAM; WGS (default), exome, gene panel or amplicon panel (--assay)
 *   -> BWA-MEM -> samtools fixmate/sort/markdup/index (no duplicate marking for amplicons)
 *   -> samtools ampliconclip                    primer removal (amplicon, --primers_bed)
 *   -> samtools stats/depth                     read QC, target coverage, uniformity, fold-80
 *   -> GATK HaplotypeCaller (or DeepVariant)    SNVs/small indels
 *   -> GATK Mutect2 (optional, --mosaic)        low-allele-fraction (mosaic) candidates
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
params.sv_caller       = 'delly'   // structural variants: 'delly' (default, native on Apple Silicon) or 'manta'
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
params.delly     = 'delly'
params.cnvkit    = 'cnvkit.py'
params.rscript   = null       // Rscript with Bioconductor DNAcopy for CNVkit; default: next to cnvkit, else PATH
params.vep       = 'vep'
params.mosdepth  = 'mosdepth'
params.python    = 'python3'
params.report_lang = 'en'      // first language of brca12_report.html/.md (both carry English and German): en or de

// Assay. 'wgs' (default) or targeted data: 'exome', 'panel' (hybrid capture) or 'amplicon'.
params.assay            = 'wgs'
params.targets_bed      = null     // target intervals; required for exome, panel and amplicon
params.primers_bed      = null     // amplicon primers (BED, optional strand in column 6) for samtools ampliconclip
params.min_target_depth = 30       // targets with a lower mean depth are listed as weakly covered
params.cnv_reference    = null     // CNVkit reference (.cnn) built from normal samples of the same assay
// Exon copy number against a reference cohort. WGS default: assets/exon_cnv_panel_1kg30x.json (80 1000 Genomes
// 30x samples); other assays need a panel built from normal samples of the same assay; 'none' switches it off.
params.exon_cnv_panel   = null
// Low-allele-fraction (mosaic) candidates with GATK Mutect2 in tumour-only mode.
params.mosaic           = false
params.mosaic_min_af    = 0.02
params.mosaic_max_af    = 0.30
params.mosaic_min_alt_reads = 3
params.mosaic_min_depth = 20
// Sequencing platform: 'illumina' (default) or 'iontorrent' (read group PL; indels flagged, mosaic SNVs only).
params.platform         = 'illumina'

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
    path targets_in
    path primers_in

    output:
    path 'reference_check.json', emit: check
    path 'regions.calling.bed',  emit: calling_bed
    path 'regions.sv.bed',       emit: sv_bed
    path 'regions.exact.bed',    emit: exact_bed
    path 'exons.bed',            emit: exons
    path 'validate.command.sh'
    path 'targets.brca.bed',     emit: targets, optional: true
    path 'primers.brca.bed',     emit: primers, optional: true

    script:
    def targets = targets_in.name != 'NO_TARGETS' ? "--targets ${targets_in} --out-targets targets.brca.bed" : ''
    def primers = primers_in.name != 'NO_PRIMERS' ? "--primers ${primers_in} --out-primers primers.brca.bed" : ''
    """
    ${params.python} "${projectDir}/bin/brca12_prepare.py" \\
        --fai ${reference_fai} --regions ${regions_bed} --exons ${exons_bed} ${targets} ${primers} \\
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
    def rg = "@RG\\tID:${params.sample}\\tSM:${params.sample}\\tLB:${params.sample}\\tPL:${params.platform == 'iontorrent' ? 'IONTORRENT' : 'ILLUMINA'}"
    // Single-end reads (e.g. Ion Torrent or R1 only): no mate information, fixmate is skipped.
    def single = r2.name == 'NO_R2'
    """
    set -euo pipefail
    if [ "${single}" = "true" ]; then
      ${params.bwa} mem -t ${task.cpus} -K 100000000 -Y -R '${rg}' ${reference} ${r1} \\
        | ${params.samtools} sort -@ ${task.cpus} -m 1G -o ${params.sample}.sorted.bam -
    else
      ${params.bwa} mem -t ${task.cpus} -K 100000000 -Y -R '${rg}' ${reference} ${r1} ${r2} \\
        | ${params.samtools} fixmate -m -@ ${task.cpus} - - \\
        | ${params.samtools} sort -@ ${task.cpus} -m 1G -o ${params.sample}.sorted.bam -
    fi
    if [ "${params.assay}" = "amplicon" ]; then
      # Amplicon reads share their start positions by design; duplicate marking would remove almost all.
      mv ${params.sample}.sorted.bam ${params.sample}.markdup.bam
      echo "duplicate marking skipped (assay: amplicon)" > ${params.sample}.markdup.stats.txt
    else
      ${params.samtools} markdup -@ ${task.cpus} -s -f ${params.sample}.markdup.stats.txt \\
        ${params.sample}.sorted.bam ${params.sample}.markdup.bam
    fi
    ${params.samtools} index -@ ${task.cpus} ${params.sample}.markdup.bam
    rm -f ${params.sample}.sorted.bam
    cp .command.sh align.command.sh
    """
}

process PRIMER_CLIP {
    tag "${params.sample}"
    cpus params.threads
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '*.{bam,bai,txt,sh}'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path calling_bed
    path primers

    output:
    tuple path("${params.sample}.primerclip.bam"), path("${params.sample}.primerclip.bam.bai"), emit: bam
    path 'primer_clip.stats.txt', emit: stats
    path 'primer_clip.command.sh'

    script:
    // Strand-aware clipping when the primer BED has a strand column.
    """
    set -euo pipefail
    STRAND=\$(awk -F'\\t' 'NF>=6 && (\$6=="+" || \$6=="-"){n++} END{print (n==NR && NR>0) ? "--strand" : ""}' ${primers})
    ${params.samtools} view -u -M -L ${calling_bed} --reference ${reference} ${bam} \\
      | ${params.samtools} ampliconclip --both-ends \$STRAND --filter-len 30 -b ${primers} \\
          -f primer_clip.stats.txt --primer-counts primer_clip.counts.txt -u - \\
      | ${params.samtools} sort -@ ${task.cpus} -o ${params.sample}.primerclip.bam -
    ${params.samtools} index ${params.sample}.primerclip.bam
    cp .command.sh primer_clip.command.sh
    """
}

process READ_QC {
    tag "${params.sample}"
    cpus 4
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path calling_bed
    path extras

    output:
    path 'read_qc.json', emit: json
    path 'region.samtools.stats.txt'
    path 'read_qc.command.sh'

    script:
    // extras: targets.brca.bed and primer_clip.stats.txt when present (plus a placeholder).
    """
    set -euo pipefail
    REGIONS=\$(awk '{printf "%s:%d-%d ", \$1, \$2 + 1, \$3}' ${calling_bed})
    ${params.samtools} stats -@ ${task.cpus} --reference ${reference} ${bam} \$REGIONS > region.samtools.stats.txt
    ${params.samtools} idxstats ${bam} > idxstats.txt || true
    OPTS=""
    if [ -s targets.brca.bed ]; then
      ${params.samtools} depth -a -Q 20 -b targets.brca.bed --reference ${reference} ${bam} | gzip -c > target_depth.tsv.gz
      ${params.samtools} view -c -F 0x904 -M -L targets.brca.bed --reference ${reference} ${bam} > on_target_reads.txt
      OPTS="--targets targets.brca.bed --target-depth target_depth.tsv.gz --on-target-reads on_target_reads.txt"
    fi
    if [ -s primer_clip.stats.txt ]; then OPTS="\$OPTS --clip-stats primer_clip.stats.txt"; fi
    if [ -s idxstats.txt ]; then OPTS="\$OPTS --idxstats idxstats.txt"; fi
    # header (read-group platform, duplicate-marking step) and read names of the first calling region
    # (read-name format of the sequencing platform); both only feed QC notes
    ${params.samtools} view -H ${bam} > bam_header.txt 2>/dev/null || true
    FIRST=\$(awk 'NR == 1 {printf "%s:%d-%d", \$1, \$2 + 1, \$3}' ${calling_bed})
    ( ${params.samtools} view --reference "${reference}" ${bam} \$FIRST 2>/dev/null || true ) | head -1000 | cut -f1 > read_names.txt || true
    if [ -s bam_header.txt ]; then OPTS="\$OPTS --header bam_header.txt"; fi
    if [ -s read_names.txt ]; then OPTS="\$OPTS --read-names read_names.txt"; fi
    ${params.python} "${projectDir}/bin/brca12_read_qc.py" --assay ${params.assay} --stats region.samtools.stats.txt \\
      --min-target-depth ${params.min_target_depth} \$OPTS --out read_qc.json
    cp .command.sh read_qc.command.sh
    """
}

process CALL_MOSAIC {
    tag "${params.sample}"
    cpus params.threads
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path calling_bed

    output:
    tuple path("${params.sample}.mosaic.vcf.gz"), path("${params.sample}.mosaic.vcf.gz.tbi"), emit: vcf
    path 'mosaic.command.sh'

    script:
    // GATK reads a BAM of the calling regions written by samtools, for two reasons: htsjdk fails on some
    // CRAMs when it converts the .crai ("Attempt to create a bai entry for an unmapped slice", e.g. the
    // 1000 Genomes 30x CRAMs), and soft-clipped poly-G read-through tails (two-colour chemistry) are
    // hard-clipped first, because GATK assembles soft-clipped bases and calls false variants from them.
    // Tumour-only Mutect2 reports variants at any allele fraction; germline heterozygous calls are
    // removed later by allele fraction (--mosaic_max_af) and by comparison with the germline calls.
    def amplicon = params.assay == 'amplicon' ? '--max-reads-per-alignment-start 0' : ''
    """
    set -euo pipefail
    ${params.samtools} view -h -M -T ${reference} -L ${calling_bed} ${bam} \\
      | ${params.python} "${projectDir}/bin/brca12_clip_polyg.py" \\
      | ${params.samtools} view -b -o calling_regions.bam -
    ${params.samtools} index calling_regions.bam
    ${params.gatk} --java-options "-Xmx4g" Mutect2 -R ${reference} -I calling_regions.bam -L ${calling_bed} \\
      --native-pair-hmm-threads ${task.cpus} ${amplicon} -O mosaic.unfiltered.vcf.gz
    ${params.gatk} --java-options "-Xmx4g" FilterMutectCalls -R ${reference} -V mosaic.unfiltered.vcf.gz \\
      -O ${params.sample}.mosaic.vcf.gz
    cp .command.sh mosaic.command.sh
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
    // HaplotypeCaller reads a region BAM with hard-clipped poly-G tails, see CALL_MOSAIC.
    if (params.snv_caller == 'deepvariant')
        """
        set -euo pipefail
        ${params.deepvariant} --model_type=${params.assay == 'wgs' ? 'WGS' : 'WES'} --ref=${reference} --reads=${bam} \\
          --regions=${calling_bed} --output_vcf=${params.sample}.small.vcf.gz \\
          --num_shards=${task.cpus}
        ${params.tabix} -f -p vcf ${params.sample}.small.vcf.gz
        cp .command.sh small_variants.command.sh
        """
    else
        """
        set -euo pipefail
        ${params.samtools} view -h -M -T ${reference} -L ${calling_bed} ${bam} \\
          | ${params.python} "${projectDir}/bin/brca12_clip_polyg.py" \\
          | ${params.samtools} view -b -o calling_regions.bam -
        ${params.samtools} index calling_regions.bam
        ${params.gatk} --java-options "-Xmx4g" HaplotypeCaller \\
          -R ${reference} -I calling_regions.bam -L ${calling_bed} \\
          --native-pair-hmm-threads ${task.cpus} ${params.assay == 'amplicon' ? '--max-reads-per-alignment-start 0' : ''} \\
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
    if (params.reads_r1 && !params.reads_r2)
        // Manta needs read pairs; with single-end reads an empty SV call set is written.
        """
        set -euo pipefail
        printf '##fileformat=VCFv4.1\\n##brca12=Manta skipped: single-end reads\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${params.sample}\\n' \\
          | ${params.bgzip} -c > ${params.sample}.manta.diploidSV.vcf.gz
        ${params.tabix} -f -p vcf ${params.sample}.manta.diploidSV.vcf.gz
        cp .command.sh manta.command.sh
        """
    else
    """
    set -euo pipefail
    sort -k1,1 -k2,2n ${sv_bed} | ${params.bgzip} -c > call_regions.bed.gz
    ${params.tabix} -f -p bed call_regions.bed.gz
    ${params.manta_config} --bam ${bam} --referenceFasta ${reference} ${params.assay == 'wgs' ? '' : '--exome'} \\
      --callRegions call_regions.bed.gz --runDir manta_run
    manta_run/runWorkflow.py -j ${task.cpus}
    cp manta_run/results/variants/diploidSV.vcf.gz ${params.sample}.manta.diploidSV.vcf.gz
    cp manta_run/results/variants/diploidSV.vcf.gz.tbi ${params.sample}.manta.diploidSV.vcf.gz.tbi
    cp .command.sh manta.command.sh
    """
}

process CALL_SV_DELLY {
    tag "${params.sample}"
    cpus 2
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path reference
    path reference_index_files
    path sv_bed

    output:
    tuple path("${params.sample}.delly.sv.vcf.gz"), path("${params.sample}.delly.sv.vcf.gz.tbi"), emit: vcf
    path 'delly.command.sh'

    script:
    if (params.reads_r1 && !params.reads_r2)
        // Delly needs read pairs; with single-end reads an empty SV call set is written.
        """
        set -euo pipefail
        printf '##fileformat=VCFv4.2\\n##brca12=Delly skipped: single-end reads\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${params.sample}\\n' \\
          | ${params.bgzip} -c > ${params.sample}.delly.sv.vcf.gz
        ${params.tabix} -f -p vcf ${params.sample}.delly.sv.vcf.gz
        cp .command.sh delly.command.sh
        """
    else
    """
    set -euo pipefail
    # Only the chromosomes of the gene windows are analysed; all other contigs are excluded.
    cut -f1 ${sv_bed} | sort -u > keep_contigs.txt
    awk 'NR==FNR{keep[\$1]=1; next} !(\$1 in keep){print \$1"\\t0\\t"\$2}' keep_contigs.txt ${reference}.fai > exclude.bed
    exclude=""
    if [ -s exclude.bed ]; then exclude="-x exclude.bed"; fi   # Delly rejects an empty exclude file
    ${params.delly} call -g ${reference} \$exclude -o delly.bcf ${bam}
    # Germline calls of one sample: PASS, non-reference genotype, DEL/DUP/INV/BND (insertions are handled
    # by the insertion signatures and the known-allele realignment); then, as Manta's --callRegions, both
    # ends of DEL/DUP/INV inside the search windows and at least 50 bp (shorter events are small variants).
    ${params.bcftools} view -f PASS \\
      -i 'GT!="0/0" && GT!="./." && (SVTYPE=="DEL" || SVTYPE=="DUP" || SVTYPE=="INV" || SVTYPE=="BND")' \\
      delly.bcf \\
      | ${params.python} "${projectDir}/bin/brca12_filter_sv.py" --windows ${sv_bed} --min-length 50 \\
      | ${params.bgzip} -c > ${params.sample}.delly.sv.vcf.gz
    ${params.tabix} -f -p vcf ${params.sample}.delly.sv.vcf.gz
    cp .command.sh delly.command.sh
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
    path targets
    path cnv_reference

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
    // Targeted assays use CNVkit hybrid (exome, panel) or amplicon mode on the target BED. With
    // --cnv_reference a reference built from normal samples of the same assay replaces the flat one.
    def method = params.assay == 'wgs' ? 'wgs' : params.assay == 'amplicon' ? 'amplicon' : 'hybrid'
    def useReference = cnv_reference.name != 'NO_CNV_REFERENCE'
    // CNVkit segments with R (DNAcopy). Another Rscript earlier on PATH (for example a conda base
    // without DNAcopy) made batch silently skip segmentation: with -p it does not pass on errors of
    // its worker processes. Rscript is therefore resolved explicitly and checked before the run, and
    // the outputs are checked afterwards.
    def rscript = params.rscript ?: ''
    """
    set -euo pipefail
    RSCRIPT="${rscript}"
    if [ -z "\$RSCRIPT" ]; then
      CNVKIT_BIN=\$(dirname "\$(command -v ${params.cnvkit} || true)")
      if [ -x "\$CNVKIT_BIN/Rscript" ]; then RSCRIPT="\$CNVKIT_BIN/Rscript"; else RSCRIPT=Rscript; fi
    fi
    if ! "\$RSCRIPT" --no-restore --no-environ -e 'suppressMessages(library(DNAcopy))' > rscript_check.log 2>&1; then
      cat rscript_check.log >&2
      echo "ERROR: CNVkit needs R with the Bioconductor package DNAcopy; '\$RSCRIPT' cannot load it. Install DNAcopy in that R or pass --rscript /path/to/Rscript." >&2
      exit 1
    fi
    # Accessible regions only on the BRCA chromosomes (chr13/chr17): the same intervals as on the
    # whole genome, since access works per contig, but about four times faster.
    BRCA_CHROMS=\$(${params.python} -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["chromLengths"]))' ${reference_check})
    ${params.samtools} faidx ${reference} \$BRCA_CHROMS > access_ref.fa
    ${params.cnvkit} access access_ref.fa -o access.bed
    rm -f access_ref.fa access_ref.fa.fai
    ${params.python} "${projectDir}/bin/brca12_prepare.py" --cnv-targets ${reference_check} \\
      --cnv-access access.bed --out-cnv-targets wgs_regions.bed
    if [ "${useReference}" = "true" ]; then
      ${params.cnvkit} batch ${bam} --reference ${cnv_reference} --method ${method} \\
        --drop-low-coverage --output-dir cnvkit -p ${task.cpus} --rscript-path "\$RSCRIPT"
    elif [ "${method}" = "wgs" ]; then
      # In WGS mode CNVkit requires the same BED for --targets and --access.
      ${params.cnvkit} batch ${bam} --method wgs --normal --fasta ${reference} \\
        --targets wgs_regions.bed --access wgs_regions.bed --target-avg-size 1000 \\
        --drop-low-coverage \\
        --output-reference flat_reference.cnn --output-dir cnvkit -p ${task.cpus} --rscript-path "\$RSCRIPT"
    elif [ "${method}" = "hybrid" ]; then
      ${params.cnvkit} batch ${bam} --method hybrid --normal --fasta ${reference} \\
        --targets ${targets} --access wgs_regions.bed --drop-low-coverage \\
        --output-reference flat_reference.cnn --output-dir cnvkit -p ${task.cpus} --rscript-path "\$RSCRIPT"
    else
      ${params.cnvkit} batch ${bam} --method amplicon --normal --fasta ${reference} \\
        --targets ${targets} --drop-low-coverage \\
        --output-reference flat_reference.cnn --output-dir cnvkit -p ${task.cpus} --rscript-path "\$RSCRIPT"
    fi
    for f in call.cns cnr; do
      if [ ! -s "cnvkit/${bam.baseName}.\$f" ]; then
        echo "ERROR: CNVkit batch finished without cnvkit/${bam.baseName}.\$f. With -p CNVkit does not report errors of its worker processes; run '${params.cnvkit} segment' on the .cnr to see the cause." >&2
        exit 1
      fi
    done
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
    path exon_cnv_panel

    output:
    path 'exons_q20.regions.bed.gz', emit: q20
    path 'exons_q20.thresholds.bed.gz', emit: thresholds
    path 'exons_q0.regions.bed.gz', emit: q0
    path 'gene_q20.per-base.bed.gz', emit: perbase
    path 'chrom_median.txt', emit: chrom_median
    path 'coverage.command.sh'
    path 'insertion_signals.json', emit: insertion
    path 'exon_cnv.json', emit: exon_cnv

    script:
    """
    set -euo pipefail
    # All reads that overlap an analysed interval are fetched once through the index (-M) into a
    # small BAM. Without this, mosdepth --by and samtools view -L read the whole alignment file
    # (hours for a 70 GB genome) although only the gene regions are needed. Depths are unchanged:
    # every read overlapping the intervals is in the small BAM.
    ${params.python} "${projectDir}/bin/brca12_exon_cnv.py" --panel ${exon_cnv_panel} --write-bed exon_cnv_regions.bed
    cat ${exons} ${exact_bed} exon_cnv_regions.bed | grep -v -e '^#' -e '^track' -e '^browser' \\
      | cut -f1-3 | sort -k1,1 -k2,2n > coverage_regions.bed
    ${params.samtools} view -b -M -L coverage_regions.bed -T ${reference} -o coverage_regions.bam ${bam}
    ${params.samtools} index coverage_regions.bam
    ${params.mosdepth} -t ${task.cpus} -f ${reference} -Q 20 -n --by ${exons} --thresholds 10,20,30 exons_q20 coverage_regions.bam
    ${params.mosdepth} -t ${task.cpus} -f ${reference} -Q 0  -n --by ${exons} exons_q0 coverage_regions.bam
    # Per-base depth inside BRCA1/BRCA2 only (MAPQ >= 20, no duplicates) for uncertain intervals.
    ${params.samtools} depth -a -Q 20 -b ${exact_bed} --reference ${reference} coverage_regions.bam \\
      | awk 'BEGIN{OFS="\\t"} {print \$1, \$2-1, \$2, \$3}' | gzip -c > gene_q20.per-base.bed.gz
    if [ "${params.assay}" = "wgs" ]; then
      # Chromosome median depth (chr13/chr17) as baseline for per-exon dosage.
      for C in \$(cut -f1 ${exact_bed} | sort -u); do
        ${params.mosdepth} -t ${task.cpus} -f ${reference} -Q 20 -n --fast-mode --by 100000 -c \$C chrom_\$C ${bam}
      done
      ${params.python} "${projectDir}/bin/brca12_prepare.py" --chrom-median chrom_*.regions.bed.gz --out-chrom-median chrom_median.txt
    else
      # Targeted data: median exon depth over both genes as baseline (a whole-gene loss stays visible).
      ${params.python} "${projectDir}/bin/brca12_prepare.py" --exon-baseline exons_q20.regions.bed.gz --out-chrom-median chrom_median.txt
    fi
    # Insertion signatures (Alu, extra tandem repeat copies) from soft-clipped read ends inside the genes.
    ${params.samtools} faidx ${reference} \$(awk '{print \$1":"\$2+1"-"\$3}' ${exact_bed}) > genes.fa
    ${params.samtools} view -L ${exact_bed} --reference ${reference} coverage_regions.bam \\
      | ${params.python} "${projectDir}/bin/brca12_insertion_signals.py" --regions ${exact_bed} --gene-fasta genes.fa --out insertion_signals.json
    # Exon copy number against a reference cohort: depth over exons +/- pad and control windows (MAPQ >= 20).
    if [ -s exon_cnv_regions.bed ]; then
      ${params.samtools} depth -a -Q 20 -b exon_cnv_regions.bed --reference ${reference} coverage_regions.bam \\
        | ${params.python} "${projectDir}/bin/brca12_exon_cnv.py" --panel ${exon_cnv_panel} --out exon_cnv.json
    else
      ${params.python} "${projectDir}/bin/brca12_exon_cnv.py" --panel ${exon_cnv_panel} --out exon_cnv.json
    fi
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
    path mosaic_files

    output:
    path 'BRCA1.snv_indel.vep.vcf.gz', emit: brca1
    path 'BRCA2.snv_indel.vep.vcf.gz', emit: brca2
    path 'BRCA*.mosaic.vep.vcf.gz', emit: mosaic, optional: true
    path 'annotate.command.sh'

    script:
    def cacheVersion = params.vep_cache_version ? "--cache_version ${params.vep_cache_version}" : ''
    """
    set -euo pipefail
    ${params.bcftools} norm -m -any -f ${reference} ${vcf} -Oz -o norm.vcf.gz
    ${params.tabix} -f -p vcf norm.vcf.gz
    SETS="snv_indel:norm.vcf.gz"
    # Optional mosaic candidates (Mutect2): only PASS calls, normalised the same way.
    MOSAIC=\$(ls *.mosaic.vcf.gz 2>/dev/null | head -1 || true)
    if [ -n "\$MOSAIC" ]; then
      ${params.bcftools} view -f PASS \$MOSAIC | ${params.bcftools} norm -m -any -f ${reference} -Oz -o mosaic.norm.vcf.gz -
      ${params.tabix} -f -p vcf mosaic.norm.vcf.gz
      SETS="\$SETS mosaic:mosaic.norm.vcf.gz"
    fi
    for SET in \$SETS; do
    KIND=\${SET%%:*}; SRC=\${SET#*:}
    for GENE in BRCA1 BRCA2; do
      awk -v g=\$GENE 'BEGIN{OFS="\\t"} \$4==g {print \$1,\$2,\$3,\$4}' ${exact_bed} > \$GENE.bed
      ${params.bcftools} view -R \$GENE.bed \$SRC -Oz -o \$GENE.\$KIND.vcf.gz
      ${params.tabix} -f -p vcf \$GENE.\$KIND.vcf.gz
      ${params.vep} --offline --cache --dir_cache ${vep_cache} ${cacheVersion} --assembly GRCh38 \\
        --fasta ${reference} --format vcf --vcf --compress_output bgzip --force_overwrite \\
        --symbol --hgvs --canonical --mane --variant_class --numbers --biotype \\
        --af_gnomade --af_gnomadg --max_af --sift b --polyphen b --flag_pick \\
        --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \\
        --custom file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \\
        --fork ${task.cpus} \\
        -i \$GENE.\$KIND.vcf.gz -o \$GENE.\$KIND.vep.vcf.gz
    done
    done
    cp .command.sh annotate.command.sh
    """
}

process TOOL_VERSIONS {
    input:
    tuple path(clinvar), path(clinvar_tbi)

    output:
    path 'tool_versions.tsv', emit: versions

    script:
    """
    # Each tool starts its own runtime (Java for GATK, Perl for VEP, R): queried in parallel, written in order.
    echo -e "nextflow\\t${workflow.nextflow.version}" > tv_00
    ( echo -e "bwa\\t\$( (${params.bwa} 2>&1 || true) | awk '/^Version/{print \$2}')" ) > tv_01 &
    ( echo -e "samtools\\t\$(${params.samtools} --version 2>/dev/null | head -1 | awk '{print \$2}')" ) > tv_02 &
    ( echo -e "bcftools\\t\$(${params.bcftools} --version 2>/dev/null | head -1 | awk '{print \$2}')" ) > tv_03 &
    if [ "${params.snv_caller}" = "deepvariant" ]; then
      ( echo -e "deepvariant\\t\$(${params.deepvariant} --version 2>&1 | tail -1 | awk '{print \$NF}')" ) > tv_04 &
    else
      ( echo -e "gatk\\t\$(${params.gatk} --version 2>/dev/null | awk '/The Genome Analysis Toolkit/{print \$NF}')" ) > tv_04 &
    fi
    ( echo -e "manta\\t\$(${params.manta_config} --version 2>/dev/null | tail -1)" ) > tv_05 &
    ( echo -e "delly\\t\$(${params.delly} --version 2>&1 | awk '/Delly version/{print \$NF; exit}')" ) > tv_06 &
    ( echo -e "cnvkit\\t\$(${params.cnvkit} version 2>/dev/null | tail -1)" ) > tv_07 &
    (
      R_FOR_CNVKIT="${params.rscript ?: ''}"; [ -n "\$R_FOR_CNVKIT" ] || { d=\$(dirname "\$(command -v ${params.cnvkit} || true)"); [ -x "\$d/Rscript" ] && R_FOR_CNVKIT="\$d/Rscript" || R_FOR_CNVKIT=Rscript; }
      echo -e "dnacopy\\t\$("\$R_FOR_CNVKIT" --no-restore --no-environ -e 'cat(as.character(packageVersion("DNAcopy")))' 2>/dev/null || true)"
    ) > tv_08 &
    ( echo -e "ensembl-vep\\t\$(${params.vep} --help 2>/dev/null | awk '/ensembl-vep/{print \$NF; exit}')" ) > tv_09 &
    ( echo -e "mosdepth\\t\$(${params.mosdepth} --version 2>/dev/null | awk '{print \$2}')" ) > tv_10 &
    ( echo -e "python\\t\$(${params.python} --version 2>&1 | awk '{print \$2}')" ) > tv_11 &
    # ClinVar release: the fileDate line of the ClinVar VCF header (for example 20251104)
    ( echo -e "clinvar\\t\$(${params.bcftools} view -h ${clinvar} 2>/dev/null | awk -F= '/^##fileDate=/{print \$2; exit}')" ) > tv_12 &
    wait
    cat tv_?? > tool_versions.tsv
    rm -f tv_??
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
    path exon_cnv
    path versions
    path command_files
    path evidence_files
    path read_qc
    path mosaic_files, stageAs: 'mosaic/*'

    output:
    path 'brca12_report.html'
    path 'brca12_report.md'
    path 'brca12_results.json'
    path 'brca12_commands.sh'

    script:
    def input_mode = params.bam ? 'bam' : 'fastq'
    def inputs = params.bam ? params.bam : (params.reads_r2 ? "${params.reads_r1},${params.reads_r2}" : "${params.reads_r1}")
    def staged = { p -> p ? file(p).name : null }
    def evidence = [['--gnomad-exomes', params.gnomad_exomes_vcf], ['--gnomad-genomes', params.gnomad_genomes_vcf],
                    ['--spliceai', params.spliceai_vcf], ['--bayesdel', params.bayesdel_tsv]]
        .findAll { it[1] }.collect { "${it[0]} ${staged(it[1])}" }.join(' ')
    def mosaicNames = (mosaic_files instanceof List ? mosaic_files : [mosaic_files]).collect { it.name }
    def mosaicArgs = ''
    if (params.mosaic) {
        def genes = ['BRCA1', 'BRCA2'].findAll { g -> mosaicNames.any { n -> n.endsWith("${g}.mosaic.vep.vcf.gz".toString()) } }
        mosaicArgs = genes.collect { g -> "--mosaic ${g}=mosaic/${g}.mosaic.vep.vcf.gz" }.join(' ') +
            " --mosaic-min-af ${params.mosaic_min_af} --mosaic-max-af ${params.mosaic_max_af}" +
            " --mosaic-min-alt-reads ${params.mosaic_min_alt_reads} --mosaic-min-depth ${params.mosaic_min_depth}"
    }
    """
    cat ${command_files} > brca12_commands.sh
    ${params.python} "${projectDir}/bin/brca12_report.py" \\
      --sample ${params.sample} --input-mode ${input_mode} --inputs '${inputs}' \\
      --reference '${params.reference}' --reference-check ${reference_check} \\
      --exons ${exons} --regions ${exact_bed} \\
      --small BRCA1=${brca1_vcf} --small BRCA2=${brca2_vcf} \\
      --sv ${manta_vcf} --cnv ${cnv_cns} --cnr ${cnv_cnr} \\
      --cov-q20 ${cov_q20} --cov-thresholds ${cov_thresholds} --cov-q0 ${cov_q0} \\
      --cov-perbase ${cov_perbase} --chrom-median ${chrom_median} \\
      --insertion-signals ${insertion_signals} --exon-cnv ${exon_cnv} \\
      --assay ${params.assay} --platform ${params.platform} --read-qc ${read_qc} ${mosaicArgs} \\
      --versions ${versions} --commands brca12_commands.sh \\
      --snv-caller ${params.snv_caller} \\
      --min-exon-depth ${params.min_exon_depth} --min-exon-pct20 ${params.min_exon_pct20} \\
      --low-depth ${params.low_depth} --max-mapq0 ${params.max_mapq0} \\
      --nextflow-version ${workflow.nextflow.version} --session-id ${workflow.sessionId} \\
      --lang ${params.report_lang} --tabix ${params.tabix} ${evidence} \\
      --out-json brca12_results.json --out-md brca12_report.md --out-html brca12_report.html
    """
}

process XLSX_EXPORT {
    // Excel export of the result (bin/brca12_xlsx.py, no third-party packages).
    tag "${params.sample}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path results_json

    output:
    path 'brca12_report.xlsx'

    script:
    """
    ${params.python} "${projectDir}/bin/brca12_xlsx.py" ${results_json} brca12_report.xlsx
    """
}

workflow {
    requireParam('reference')
    requireParam('vep_cache')
    requireParam('clinvar_vcf')
    // the sample name becomes part of file names and of the read group
    if (!(params.sample.toString() ==~ /[A-Za-z0-9._-]+/)) {
        error "--sample may only contain letters, digits, dot, underscore and hyphen (got '${params.sample}')."
    }
    if (!params.bam && !params.reads_r1) {
        error 'Provide either --reads_r1 (and --reads_r2 for paired-end reads) or --bam.'
    }

    reference = file(params.reference, checkIfExists: true)
    fai = file("${params.reference}.fai", checkIfExists: true)
    ref_files = Channel.fromPath([
        "${params.reference}.fai",
        "${params.reference}.{amb,ann,bwt,pac,sa}",
        "${params.reference.replaceAll(/\.(fa|fasta)(\.gz)?$/, '')}.dict"
    ].flatten(), checkIfExists: false).collect()

    if (!(params.assay in ['wgs', 'exome', 'panel', 'amplicon'])) {
        error "--assay must be wgs, exome, panel or amplicon (got '${params.assay}')."
    }
    if (!(params.sv_caller in ['delly', 'manta'])) {
        error "--sv_caller must be delly or manta (got '${params.sv_caller}')."
    }
    if (params.assay != 'wgs' && !params.targets_bed) {
        error "--assay ${params.assay} needs --targets_bed (the capture or amplicon intervals)."
    }
    if (params.primers_bed && params.assay != 'amplicon') {
        error '--primers_bed is only used with --assay amplicon.'
    }
    targets_in = file(params.targets_bed ?: "${projectDir}/assets/NO_TARGETS", checkIfExists: true)
    primers_in = file(params.primers_bed ?: "${projectDir}/assets/NO_PRIMERS", checkIfExists: true)
    VALIDATE_INPUTS(fai, file(params.regions_bed, checkIfExists: true), file(params.exons_bed, checkIfExists: true),
                    targets_in, primers_in)
    targets_ch = VALIDATE_INPUTS.out.targets.ifEmpty(file("${projectDir}/assets/NO_TARGETS"))

    if (params.bam) {
        // Index next to the alignment: BAM -> .bam.bai, .bai or .bam.csi; CRAM -> .cram.crai or .crai.
        def stem = "${params.bam}".replaceAll(/\.(bam|cram)$/, '')
        def candidates = "${params.bam}".endsWith('.cram') ? ["${params.bam}.crai", "${stem}.crai"]
                       : ["${params.bam}.bai", "${stem}.bai", "${params.bam}.csi"]
        def found = candidates.find { file(it).exists() }
        if (!found) { error "No index for ${params.bam} (looked for ${candidates.join(', ')})." }
        def bai = file(found)
        def bam_in = file(params.bam, checkIfExists: true)
        // A file name with spaces or other shell characters would break the tool commands: such a BAM or
        // CRAM is linked under a safe name (only letters, digits, dot, underscore, hyphen); data are not copied.
        if (!(bam_in.name ==~ /[A-Za-z0-9._-]+/)) {
            def dir = file("${workflow.workDir}/input-links/${workflow.sessionId}")
            dir.mkdirs()
            def safe = bam_in.name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            def link = dir.resolve(safe)
            def bai_link = dir.resolve(bai.name.startsWith(bam_in.name) ? safe + bai.name.substring(bam_in.name.length())
                                                                         : bai.name.replaceAll(/[^A-Za-z0-9._-]/, '_'))
            bam_in.mklink(link, overwrite: true)
            bai.mklink(bai_link, overwrite: true)
            log.info "Input ${bam_in.name} is linked as ${safe} (file names with spaces or shell characters)."
            bam_in = link
            bai = bai_link
        }
        aligned_ch = Channel.value(tuple(bam_in, bai))
    } else {
        r2 = params.reads_r2 ? file(params.reads_r2, checkIfExists: true) : file("${projectDir}/assets/NO_R2")
        reads = Channel.value(tuple(file(params.reads_r1, checkIfExists: true), r2))
        aligned_ch = ALIGN_BWA_MEM(reads, reference, ref_files).bam
    }
    if (params.primers_bed) {
        PRIMER_CLIP(aligned_ch, reference, ref_files, VALIDATE_INPUTS.out.calling_bed, VALIDATE_INPUTS.out.primers)
        bam_ch = PRIMER_CLIP.out.bam
        clip_stats = PRIMER_CLIP.out.stats
    } else {
        bam_ch = aligned_ch
        clip_stats = Channel.empty()
    }
    qc_extras = Channel.empty().mix(VALIDATE_INPUTS.out.targets, clip_stats)
        .mix(Channel.fromPath("${projectDir}/assets/NO_EVIDENCE")).collect()
    READ_QC(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.calling_bed, qc_extras)

    small = CALL_SMALL_VARIANTS(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.calling_bed)
    if (params.sv_caller == 'manta') {
        sv = CALL_SV_MANTA(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.sv_bed)
        sv_cmd = CALL_SV_MANTA.out[1]
    } else {
        sv = CALL_SV_DELLY(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.sv_bed)
        sv_cmd = CALL_SV_DELLY.out[1]
    }
    cnv = CALL_CNV_CNVKIT(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.check, targets_ch,
                          file(params.cnv_reference ?: "${projectDir}/assets/NO_CNV_REFERENCE", checkIfExists: true))
    exon_cnv_panel_in = file(params.exon_cnv_panel && params.exon_cnv_panel != 'none' ? params.exon_cnv_panel
        : "${projectDir}/assets/" + (params.exon_cnv_panel != 'none' && params.assay == 'wgs' ? 'exon_cnv_panel_1kg30x.json' : 'NO_EXON_CNV_PANEL'),
        checkIfExists: true)
    cov = EXON_COVERAGE(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.exons, VALIDATE_INPUTS.out.exact_bed,
                        exon_cnv_panel_in)
    clinvar = Channel.value(tuple(file(params.clinvar_vcf, checkIfExists: true), file("${params.clinvar_vcf}.tbi", checkIfExists: true)))
    if (params.mosaic) {
        CALL_MOSAIC(bam_ch, reference, ref_files, VALIDATE_INPUTS.out.calling_bed)
        mosaic_in = CALL_MOSAIC.out.vcf.map { it.toList() }
        mosaic_cmd = CALL_MOSAIC.out[1]
    } else {
        mosaic_in = Channel.value([file("${projectDir}/assets/NO_MOSAIC")])
        mosaic_cmd = Channel.empty()
    }
    ann = EXTRACT_AND_ANNOTATE(small.vcf, VALIDATE_INPUTS.out.exact_bed, reference, ref_files,
                               file(params.vep_cache, checkIfExists: true), clinvar, mosaic_in)
    versions = TOOL_VERSIONS(clinvar)

    commands = Channel.empty()
        .mix(VALIDATE_INPUTS.out[5])
        .mix(params.bam ? Channel.empty() : ALIGN_BWA_MEM.out[2])
        .mix(params.primers_bed ? PRIMER_CLIP.out[2] : Channel.empty(), READ_QC.out[2], mosaic_cmd)
        .mix(CALL_SMALL_VARIANTS.out[1], sv_cmd, CALL_CNV_CNVKIT.out[2],
             EXON_COVERAGE.out[5], EXTRACT_AND_ANNOTATE.out[2])
        .collect()

    // Optional evidence files and their tabix indexes; an empty placeholder keeps the input non-empty.
    evidence_files = Channel.fromPath(
        [params.gnomad_exomes_vcf, params.gnomad_genomes_vcf, params.spliceai_vcf, params.bayesdel_tsv]
            .findAll { it }.collectMany { [it, "${it}.tbi"] } + ["${projectDir}/assets/NO_EVIDENCE"],
        checkIfExists: true).collect()

    REPORT(VALIDATE_INPUTS.out.check, VALIDATE_INPUTS.out.exons, VALIDATE_INPUTS.out.exact_bed,
           ann.brca1, ann.brca2, sv.vcf.map { it[0] }, cnv.cns, cnv.cnr,
           cov.q20, cov.thresholds, cov.q0, cov.perbase, cov.chrom_median, cov.insertion, cov.exon_cnv,
           versions.versions, commands, evidence_files, READ_QC.out.json,
           ann.mosaic.flatten().unique { it.name }.collect().ifEmpty([file("${projectDir}/assets/NO_MOSAIC")]))
    XLSX_EXPORT(REPORT.out[2])
}
