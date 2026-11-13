#!/usr/bin/env python3
import sys, os, gzip, shutil
name = os.path.basename(sys.argv[0]); a = sys.argv[1:]
def opt(flag):
    return a[a.index(flag)+1] if flag in a else None
def touch(p): open(p,'a').close()
CSQF="Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|VARIANT_CLASS|CANONICAL|MANE_SELECT|PICK|MAX_AF|gnomADe_AF|gnomADg_AF|ClinVar|ClinVar_CLNSIG|ClinVar_CLNREVSTAT|ClinVar_CLNDN|ClinVar_CLNSIGCONF"
VCF="##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr17\t43093000\t.\tAG\tA\t500\tPASS\tMQ=60\tGT:AD:DP:GQ\t0/1:20,18:38:99\nchr13\t32338000\t.\tA\tG\t900\tPASS\tMQ=60\tGT:AD:DP:GQ\t1/1:0,35:35:99\n"
if name=="bwa":
    if not a: print("Version: 0.7.18-r1243", file=sys.stderr); sys.exit(1)
    print("@HD\tVN:1.6")
elif name=="samtools":
    if a[0]=="--version": print("samtools 1.21")
    elif a[0]=="fixmate": shutil.copyfileobj(sys.stdin, sys.stdout)
    elif a[0]=="sort": open(opt("-o"),"w").write(sys.stdin.read())
    elif a[0]=="markdup": open(a[-1],"w").write("bam"); open(opt("-f"),"w").write("stats")
    elif a[0]=="index": touch(a[-1]+".bai")
    elif a[0]=="depth":
        for p in range(43045001,43045101): print(f"chr17\t{p}\t{5 if p<43045041 else 25}")
elif name=="gatk":
    if a[0]=="--version": print("The Genome Analysis Toolkit (GATK) v4.6.1.0"); sys.exit(0)
    o=opt("-O"); open(o,"w").write(VCF); touch(o+".tbi")
elif name=="bcftools":
    if a[0]=="--version": print("bcftools 1.21"); sys.exit(0)
    o=opt("-o"); src=[x for x in a if x.endswith(".vcf.gz") and x!=o][-1]
    lines=open(src).read().splitlines(True)
    if a[0]=="view":
        bed=open(opt("-R")).read().split()
        chrom=bed[0]
        lines=[l for l in lines if l.startswith("#") or l.startswith(chrom+"\t")]
    open(o,"w").writelines(lines)
elif name=="tabix": touch(a[-1]+".tbi")
elif name=="bgzip": shutil.copyfileobj(sys.stdin, sys.stdout)
elif name=="vep":
    if a[0]=="--help": print("Versions:\n  ensembl-vep          : 113.0"); sys.exit(0)
    src,o=opt("-i"),opt("-o"); out=[]
    for l in open(src):
        if l.startswith("#CHROM"):
            out.append(f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {CSQF}">\n'); out.append(l)
        elif l.startswith("#"): out.append(l)
        else:
            f=l.rstrip("\n").split("\t"); g="BRCA1" if f[0]=="chr17" else "BRCA2"
            c=["-","frameshift_variant","HIGH",g,"E","Transcript","ENST00000357654","protein_coding","10/23","","c.3000del","p.Ser1000fs","deletion","YES","NM_007294.4","1","","","","","Pathogenic","reviewed_by_expert_panel","Hereditary_breast_ovarian_cancer_syndrome",""] if g=="BRCA1" else ["G","synonymous_variant","LOW",g,"E","Transcript","ENST00000380152","protein_coding","11/27","","c.4000A>G","p.=","SNV","YES","NM_000059.4","1","0.3","","","","Benign","reviewed_by_expert_panel","not_specified",""]
            f[7]+=";CSQ="+"|".join(c); out.append("\t".join(f)+"\n")
    open(o,"w").writelines(out)
elif name=="configManta.py":
    if "--version" in a: print("1.6.0"); sys.exit(0)
    rd=opt("--runDir"); os.makedirs(rd,exist_ok=True)
    w=os.path.join(rd,"runWorkflow.py")
    open(w,"w").write("#!/bin/sh\nmkdir -p manta_run/results/variants\nprintf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tS\\nchr17\\t43090000\\tMantaDEL:1\\tN\\t<DEL>\\t200\\tPASS\\tSVTYPE=DEL;END=43096000;SVLEN=-6000\\tGT:PR:SR\\t0/1:20,8:15,6\\n' > manta_run/results/variants/diploidSV.vcf.gz\ntouch manta_run/results/variants/diploidSV.vcf.gz.tbi\n")
    os.chmod(w,0o755)
elif name=="cnvkit.py":
    if a[0]=="version": print("0.9.11"); sys.exit(0)
    if a[0]=="access":
        open(opt("-o"),"w").write("chr1\t0\t1000\nchr13\t10000\t114000000\nchr17\t60000\t83000000\n"); sys.exit(0)
    if opt("--access") != opt("--targets"): sys.exit("cnvkit: WGS mode needs the same --targets and --access")
    d=opt("--output-dir"); os.makedirs(d,exist_ok=True); stem=os.path.basename(a[1])[:-4]
    open(f"{d}/{stem}.call.cns","w").write("chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\nchr17\t40000000\t46000000\t-\t0.0\t2\t30\t6000\t1\nchr13\t30000000\t35000000\t-\t0.01\t2\t30\t5000\t1\n")
    open(f"{d}/{stem}.cnr","w").write("chromosome\tstart\tend\tgene\tdepth\tlog2\tweight\nchr17\t43090000\t43091000\t-\t30\t0.0\t1\n")
elif name=="mosdepth":
    if a[0]=="--version": print("mosdepth 0.3.10"); sys.exit(0)
    by=opt("--by"); prefix=a[-2]
    if by.isdigit():
        c=opt("-c"); 
        with gzip.open(prefix+".regions.bed.gz","wt") as h:
            for i in range(5): h.write(f"{c}\t{i*100000}\t{(i+1)*100000}\t30.00\n")
    else:
        rows=[l.split("\t") for l in open(by) if l.strip()]
        q0 = opt("-Q")=="0"
        with gzip.open(prefix+".regions.bed.gz","wt") as h:
            for r in rows:
                m = 8.0 if r[3]=="BRCA1_ex23" else 31.0
                if q0: m = m*1.5 if r[3]=="BRCA1_ex23" else m
                h.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{m}\n")
        if "--thresholds" in a:
            with gzip.open(prefix+".thresholds.bed.gz","wt") as h:
                h.write("#chrom\tstart\tend\tregion\t10X\t20X\t30X\n")
                for r in rows:
                    L=int(r[2])-int(r[1]); h.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{L//2 if r[3]=='BRCA1_ex23' else L}\t{0 if r[3]=='BRCA1_ex23' else L}\t0\n")
