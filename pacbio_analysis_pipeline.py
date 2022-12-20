
# import dependencies
import os
import sys
import subprocess
import pandas as pd
from pathlib import Path
import shutil
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

import matplotlib
matplotlib.use('Agg')

def getArgs():
    """
    function to obtain the commandline arguments
    """
    
    assay_dir = sys.argv[1]
    assembler = sys.argv[2]
    nproc = sys.argv[3]
    return assay_dir, assembler, nproc

def runCMD(cmd):
    print("==============================================================================")
    print(cmd)
    print("==============================================================================")
    subprocess.run(cmd, shell=True, executable="/bin/bash")


def unzipData(raw_dir, sample):
    """
    function to check if data is unzip
    and unzip
    """

    # obtain samples files
    sampleFiles = [file for file in os.listdir(raw_dir) if ~file.find(sample)]

    # check if data is already cleaned
    filtered = [file for file in sampleFiles if ~file.find("filtered")]
    unzipped = [file for file in sampleFiles if ~file.endswith("gz")]
    zipped = [file for file in sampleFiles if file.endswith("gz")]
    if len(filtered) == 1:
        print(f"raw data already clean, skipping unzipping..")
        return 0
    elif (len(filtered) == 0) and (len(unzipped) > 0):
        print(f"raw data already unzipped, skipping..")
    elif (len(zipped) == 1):
        print(f"unzipping data...")
        cmd = f"gunzip {raw_dir}/{zipped[0]}"
        runCMD(cmd)
        return 1


def createRefDB(assay_dir, refFasta, refGTF):
    """
    function to index fasta file
    and create reference db for snpEff
    """

    cmd = f"mkdir -p {assay_dir}/reference/data/athaliana; "
    cmd += f"echo \"athaliana.genome : athaliana\" > {assay_dir}/reference/snpEff.config; "
    cmd += f"cp {assay_dir}/reference/{refFasta} {assay_dir}/reference/data/athaliana/sequences.fa; "
    cmd += f"cp {assay_dir}/reference/{refGTF} {assay_dir}/reference/data/athaliana/genes.gtf; "
    cmd += f"export PATH=\"/home/nkiel/snippy/binaries/noarch/:$PATH\"; "
    cmd += f"snpEff build -v -c {assay_dir}/reference/snpEff.config -gtf22 athaliana; "
    cmd += f"samtools faidx {assay_dir}/reference/{refFasta}; "
    runCMD(cmd)


def remove2kb(raw_dir, sample, nproc):
    """
    function to check if data is unzip
    and unzip
    """

    # obtain samples files
    sampleFiles = [file for file in os.listdir(raw_dir) if ~file.find(sample)]
    filtered = [file for file in sampleFiles if ~file.find("filtered")]
    rawfasta = [file for file in sampleFiles if ~file.find(".hifi_reads.fastq")][0]
    
    # check if data is already cleaned
    if len(filtered) == 1:
        print(f"raw data already clean, skipping..")
        return 0
    
    elif len(filtered) == 0:
        print(f"raw data not clean, removing 2kb control")

        # path to fasta
        pacbio_control="/netscratch/dep_psl/grp_rgo/flores-uribe/reiso/reiso1_pacbio_bact_genomes/2kb.fasta"

        # output names
        outname = f'{raw_dir}/{sample}_filtered.fasta'
        rawname = f'{raw_dir}/{rawfasta}'

        # define subprocess
        cmd = f'blasr -nproc {nproc} {rawname} {pacbio_control} -unaligned {outname}'
        runCMD(cmd)
        return 1


def assembly(assay_dir, sample, nproc, assembler, haplotypeResolved, filterReads, assemblyDir):
    """
    perform assembly using the specified assembler
    and haplotype options
    """

    print(f"performing assembly using: {assembler}")


    if filterReads:
        rawInput = f"{assay_dir}/raw_data/{sample}_filtered.fasta"
    elif not filterReads:
        rawInput = f"{assay_dir}/raw_data/{sample}.hifi_reads.fastq"


    # define output directories
    assemblyDir = f'{assemblyDir}/{sample}'
    outdir = f'{assemblyDir}/assembly'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    if haplotypeResolved:
        assemblyOut1 = f'{assemblyDir}/{sample}_hap1.fasta'
        assemblyGraph1 = f"{assemblyDir}/{sample}_hap1_graph.gfa"
        assemblyOut2 = f'{assemblyDir}/{sample}_hap2.fasta'
        assemblyGraph2 = f"{assemblyDir}/{sample}_hap2_graph.gfa"
    
    elif not haplotypeResolved:
        assemblyOut = f'{assemblyDir}/{sample}.fasta'
        assemblyGraph = f"{assemblyDir}/{sample}_graph.gfa"
        assemblyImage = f'{assemblyDir}/{sample}.png'


    if (assembler == "flye") and not (os.path.isfile(f'{outdir}/assembly.fasta')):
        cmd = f"source /home/nkiel/activate;"
        cmd += f'conda activate svimasm_env;'
        cmd += f"flye --pacbio-hifi {rawInput} --genome-size 135M --out-dir {outdir} --threads {nproc} --iterations 4"
        runCMD(cmd)

        if not haplotypeResolved:
            # move assembly into sample folder
            shutil.copy(f'{outdir}/assembly.fasta', assemblyOut)
            shutil.copy(f'{outdir}/assembly_graph.gfa', assemblyGraph)

        if haplotypeResolved:
            # recover haplotypes from diploid assembly
            mapped = f"{outdir}/mapped_reads.bam"
            cmd = f"minimap2 -ax map-hifi -t {nproc} {outdir}/assembly.fasta {rawInput} | samtools sort -@ 4 -m 4G > {mapped};"
            cmd += f"samtools index -@ 4 {mapped};"
            cmd += f"hapdup --assembly {outdir}/assembly.fasta --bam {mapped} --out-dir {outdir} -t {nproc} --rtype hifi"
            runCMD(cmd)

            # work from here...
            # hapdup_dual_{1,2}.fasta dual assembly
            # hapdup_phased_{1,2}.fasta - haplotype-resolved assmebly
            shutil.copy(f'{outdir}/hapdup_phased_1.fasta', assemblyOut1)
            shutil.copy(f'{outdir}/hapdup_phased_1.fasta', assemblyOut2)


    if assembler == "hifiasm":
        cmd = f"source /home/nkiel/activate;"
        cmd += f'conda activate svimasm_env;'

        if (not haplotypeResolved) and not (os.path.isfile(assemblyOut)):
            cmd += f"hifiasm -l 0 -z 20 -r 4 -t {nproc} -o {outdir}/{sample} {rawInput};"
            cmd += f'awk \'/^S/{{print \">\"$2;print $3}}\' {outdir}/{sample}.bp.p_ctg.gfa > {assemblyOut}'
            runCMD(cmd)

            # move assembly graphs
            shutil.copy(f'{outdir}/{sample}.bp.p_ctg.gfa', assemblyGraph)

        if haplotypeResolved and not ((os.path.isfile(assemblyOut1)) and (os.path.isfile(assemblyOut2))):
            cmd += f"hifiasm -s 0.65 -t {nproc} -o {outdir}/{sample} {rawInput};" #--hom-cov
            cmd += f'awk \'/^S/{{print \">\"$2;print $3}}\' {outdir}/{sample}.bp.hap1.p_ctg.gfa > {assemblyOut1};'
            cmd += f'awk \'/^S/{{print \">\"$2;print $3}}\' {outdir}/{sample}.bp.hap2.p_ctg.gfa > {assemblyOut2}'
            runCMD(cmd)

            # extract assemblies
            shutil.copy(f'{outdir}/{sample}.bp.hap1.p_ctg.gfa', assemblyGraph1)
            shutil.copy(f'{outdir}/{sample}.bp.hap2.p_ctg.gfa', assemblyGraph2)


    if not haplotypeResolved:
        # create image of assembly graph
        cmd = f"Bandage image {assemblyGraph} {assemblyImage}"
        runCMD(cmd)
        return assemblyOut, assemblyDir, 1

    if haplotypeResolved:
        try:
            # create image of assembly graph
            cmd = f"cd {assemblyDir};"
            cmd += f"Bandage image {sample}_hap1_graph.gfa {sample}_hap1.png;"
            cmd += f"Bandage image {sample}_hap2_graph.gfa {sample}_hap2.png;"
            runCMD(cmd)
        except:
            pass
        return [assemblyOut1, assemblyOut2], assemblyDir, 1



def getAssemblyStats(assemblyPath, sample, assembler, assemblyDir):
    """
    function to obtain assembly statistics and make
    a coverage plot
    """

    print(f'gathering assembly stats...')

    # output paths
    haplo = os.path.basename(assemblyPath).replace(".fasta", "")
    tableOut = f"{assemblyDir}/{haplo}_coverage.csv"
    figOut = f"{assemblyDir}/{haplo}_coverage.pdf"

    # gather assembly statistucts
    cmd = f"source /home/nkiel/activate;"
    cmd += f'conda activate svimasm_env;'
    cmd += f'assembly-stats {assemblyPath} > {assemblyDir}/{haplo}_{assembler}.stats'
    runCMD(cmd)


    # get the commulative coverage
    cmd = f"source /home/nkiel/activate;"
    cmd += f'conda activate svimasm_env;'
    cmd += f"STRAIN1={haplo};TYPE1=contig;"
    cmd += f'LEN1=`bioawk -c fastx \'{{sum+=length($seq)}}END{{print sum}}\' {assemblyPath}`;'
    cmd += f'echo "line,length,type,coverage" > {tableOut};'
    cmd += f"cat {assemblyPath} | bioawk -c fastx -v line=\"$STRAIN1\" \'{{print line\",\"length($seq)\",\"length($seq)}}\' | sort -k3rV -t \",\" | awk -F \",\" -v len=\"$LEN1\" -v type=\"$TYPE1\" \'OFS=\",\"{{ print $1,$2,type,(sum+0)/len; sum+=$3 }}\' >> {tableOut}"
    runCMD(cmd)

    
    # make a plot from the assembly statistics
    coverageTag = pd.read_csv(tableOut)
    coverageTag["length"] = coverageTag["length"] / 1000000
    p = sns.lineplot(data=coverageTag, x="coverage", y = "length", hue="line", drawstyle='steps-post')
    p.set_xlabel("Cumulative coverage")
    p.set_ylabel("Length (Mb)")
    p.set_xlim(0,1)
    p.set_ylim(0,coverageTag["length"].max()*1.05)
    p.axvline(0.5, color='r', linestyle='--', lw=1)
    fig = p.get_figure()
    fig.tight_layout()
    fig.savefig(figOut)
    plt.clf()

    return 1


def runBusco(assemblyPath, sample, assemblyDir, nproc):
    """
    function to run busco to asses assembly
    quality and completeness
    """

    print(f"running busco..")

    # reference lineage
    lineage = "/netscratch/dep_psl/grp_rgo/kieln/references/databases/busco_downloads/lineages/brassicales_odb10/"

    # output directories
    haplo = os.path.basename(assemblyPath).replace(".fasta", "")
    outdir = f'{assemblyDir}/busco/{haplo}'
    Path(outdir).mkdir(parents=True, exist_ok=True)
    busco_tab=f'{outdir}/short_summary.specific..{haplo}.txt'
    busco_out=f'{assemblyDir}/{haplo}_busco.csv'
    figOut = f"{assemblyDir}/{haplo}_busco.pdf"

    if not (os.path.isfile(busco_tab)):
        # command line inputs
        cmd = f"source /home/nkiel/activate;"
        cmd += f'conda activate svimasm_env;'
        cmd += f'cd {assemblyDir}/busco;'
        cmd += f"busco -i {assemblyPath} -l {lineage} -o {haplo} -m genome -c {nproc} -f"
        runCMD(cmd)


    # make table from busco output
    cmd = f"echo \"Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing\" > {busco_out};"
    cmd += f"cat {busco_tab} | grep  \"(S)\" | awk -v strain=\"{haplo}\" \'{{print strain\",\"$1}}\' > {outdir}/complete_single.txt;"
    cmd += f"cat {busco_tab} | grep \"(D)\" | awk \'{{print $1}}\' > {outdir}/complete_duplicated.txt;"
    cmd += f"cat {busco_tab} | grep \"(F)\" | awk \'{{print $1}}\' > {outdir}/fragmented.txt;"
    cmd += f"cat {busco_tab} | grep \"(M)\" | awk \'{{print $1}}\' > {outdir}/missing.txt;"
    cmd += f"paste -d \",\" {outdir}/complete_single.txt {outdir}/complete_duplicated.txt {outdir}/fragmented.txt {outdir}/missing.txt >> {busco_out}"
    runCMD(cmd)


    # plot busco
    buscoTab = pd.read_csv(busco_out).set_index("Strain")
    plot = buscoTab.plot(kind="bar", stacked=True)
    plot.set_ylabel("BUSCO counts")
    plot.set_xlabel("")
    fig = plot.get_figure()
    fig.tight_layout()
    fig.savefig(figOut)
    plt.clf()

    return 1


def scaffolding(assemblyDir, nproc, assemblyPath, sample, ref):
    """
    function to perform scaffling 
    of the assembly using ragtag
    """

    print("performing scaffolding..")

  
    # output directories
    haplo = os.path.basename(assemblyPath).replace(".fasta", "")
    outdir = f'{assemblyDir}/scaffolding/{haplo}'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # perform scaffolding
    cmd = f"source /home/nkiel/activate;"
    cmd += f'conda activate svimasm_env;'
    cmd += f'samtools faidx {assemblyPath};'
    cmd += f"ragtag.py scaffold -t {nproc} -u -o {outdir} {ref} {assemblyPath}"
    runCMD(cmd)


    # remove ragtag-tags from filenames and header
    scaffoldIn = f'{outdir}/ragtag.scaffold.fasta'
    scaffoldOut = f'{assemblyDir}/{haplo}.scaffold.fasta'
    chrom = list(SeqIO.to_dict(SeqIO.parse(ref, "fasta")).keys())

    seqList = []
    with open(scaffoldIn) as handle:
        for seq in SeqIO.parse(handle, "fasta"):
            seqID = seq.id.replace("_RagTag", "")
            
            # append seq list
            if seqID in chrom:
                seqList.append(SeqIO.SeqRecord(seq.seq, id=seqID, description=""))

    SeqIO.write(seqList, scaffoldOut, "fasta")


    return 1, scaffoldOut


def plotSyn(scaffold, ref, sample, nproc, assemblyDir, ref_dir):
    """
    function to plot the synteny 
    of the assembled scaffolds
    """

    print("generating synteny plots..")

    # inputs
    haplo = os.path.basename(scaffold).replace(".scaffold.fasta", "")
    centromere = f"{ref_dir}/Giraut2011_centromeres.bed"

    # output paths
    outDir = f'{assemblyDir}/syri'
    samOut = f"{outDir}/{haplo}_scaffold_alignment.sam"
    bamOut = f"{outDir}/{haplo}_scaffold_alignment.bam"
    Path(outDir).mkdir(parents=True, exist_ok=True)
    genomesFile = f'{outDir}/genomes.txt'


    # variants
    variantsOut = f"{outDir}/{haplo}_syri.vcf"
    variantsFiltered = f"{outDir}/{haplo}_filtered.vcf"
    vcfClean = f"{outDir}/{haplo}_clean.vcf"
    

    # alignment using minimap2
    cmd = f"minimap2 -a -x asm5 --eqx -t {nproc} {ref} {scaffold} > {samOut}; "
    cmd += f"samtools view -b -@ {nproc} -m 8G {samOut} | samtools sort -@ {nproc} -m 8G -o {bamOut}; "
    cmd += f"samtools index {bamOut}; "
    runCMD(cmd)


    # command line inputs for syri
    cmd = f"source /home/nkiel/activate;"
    cmd += f'conda activate syenv;'
    cmd += f'cd {outDir};'
    cmd += f'syri -c {samOut} -r {ref} -q {scaffold} -k -F S --prefix {haplo}_ ; '
    cmd += f"awk \'$1 ~ /^#/ {{print $0;next}} {{if ($4 ~ /A|C|T|G/ && $5 ~ /A|C|T|G/) print $0}}\' {variantsOut} > {variantsFiltered} ;"
    cmd += f"bedtools intersect -v -a {variantsFiltered} -b {centromere} -wa > {vcfClean}"
    runCMD(cmd)


    # plot synteny using plotr
    file = [ref, scaffold]
    name = ["ref", haplo]
    pd.DataFrame({"#file": file, "name": name}).to_csv(genomesFile, sep='\t', index=False)

    # plotr cli
    sr = f'{outDir}/{haplo}_syri.out'
    plotOut = f'{assemblyDir}/{haplo}_synteny.pdf'
    markers = f'{ref_dir}/markers.bed'
    tracks = f'{ref_dir}/tracks.txt'

    cmd = f"source /home/nkiel/activate;"
    cmd += f'conda activate syenv;'
    cmd += f"plotsr --sr {sr} --genomes {genomesFile} -H 8 -o {plotOut} --markers {markers}"
    runCMD(cmd)
        
    return 1


def svCalling(assemblyDir, sample, ref, scaffold, nproc, ref_dir, haplotypeResolved):
    """
    function to perform search for structural variant calling
    """

    # log
    print("performing structural variant calling..")

    # paths
    outDir = f'{assemblyDir}/sv_calling'

    if not haplotypeResolved:
        samOut = f"{outDir}/{sample}_scaffold_alignment.bam"

    if haplotypeResolved:
        samOut1 = f"{outDir}/{sample}_hap1_scaffold_alignment.bam"
        samOut2 = f"{outDir}/{sample}_hap1_scaffold_alignment.bam"

    # output paths
    vcf = f"{outDir}/variants.vcf"
    vcfMapped = f'{assemblyDir}/{sample}_mapped.vcf'
    centromere = f"{ref_dir}/Giraut2011_centromeres.bed"
    vcfClean = f'{assemblyDir}/{sample}_mapped_clean.vcf'
    Path(outDir).mkdir(parents=True, exist_ok=True)

    # cli command
    cmd = f"source /home/nkiel/activate;"
    cmd += f'conda activate svimasm_env;'

    if not haplotypeResolved:
        cmd += f"minimap2 -a -x asm5 --cs -r2k -t {nproc} {ref} {scaffold} | samtools sort -m4G -@10 -O BAM -o {samOut};"
        cmd += f"samtools index {samOut};"
        cmd += f"svim-asm haploid --max_sv_size 2500 {outDir} {samOut} {ref};"

    if haplotypeResolved:
        cmd += f"minimap2 -a -x asm5 --cs -r2k -t {nproc} {ref} {scaffold[0]} | samtools sort -m4G -@10 -O BAM -o {samOut1};"
        cmd += f"minimap2 -a -x asm5 --cs -r2k -t {nproc} {ref} {scaffold[1]} | samtools sort -m4G -@10 -O BAM -o {samOut2};"
        cmd += f"samtools index {samOut1};"
        cmd += f"samtools index {samOut2};"
        cmd += f"svim-asm diploid --max_sv_size 2500 {outDir} {samOut1} {samOut2} {ref};"

    # execute the sv-calling
    cmd += f"export PATH=\"/home/nkiel/snippy/binaries/noarch/:$PATH\";"
    cmd += f"snpEff ann -v -noLog -noStats -no-downstream -no-upstream -c {ref_dir}/snpEff.config -csvStats {outDir}/snpEff.stats.tab athaliana {vcf} > {vcfMapped};"
    cmd += f"bedtools intersect -v -a {vcfMapped} -b {centromere} -wa > {vcfClean}"
    runCMD(cmd)

    return 1, vcfClean


def mapReads(assay_dir, filterReads, sample, nproc, ref):
    """
    function to map the reads against the reference
    and store as sorted bam file for manual inspection
    with IGV
    """

    print(f"mapping reads against reference..")
    
    # select the correct reads
    if filterReads:
        rawInput = f"{assay_dir}/raw_data/{sample}_filtered.fasta"
    elif not filterReads:
        rawInput = f"{assay_dir}/raw_data/{sample}.hifi_reads.fastq"


    # define output directories
    outdir = f'{assay_dir}/mappedReads'
    Path(outdir).mkdir(parents=True, exist_ok=True)
    mappedSam = f'{assay_dir}/mappedReads/{sample}.sam'
    mappedBam = f'{assay_dir}/mappedReads/{sample}.bam'


    # command for mapping
    cmd = f"source $HOME/activate; conda activate svimasm_env; "
    cmd += f"minimap2 -ax map-hifi -t {nproc} {ref} {rawInput} | samtools sort -@ {nproc} -m 8G > {mappedSam};"
    cmd += f"samtools view -bS -@ {nproc} {mappedSam} > {mappedBam};"
    cmd += f"samtools index -@ {nproc} {mappedBam}; "
    runCMD(cmd)

    return 1, mappedSam, mappedBam



def checkOffTargets_sample(mapping_dir, sample, assemblyDir, variantSet = "syri"):
    """
    function to map potential off targets provided 
    by jana and johannes to the syri and svimasm
    calls and check for variants up to 3 bp downstream

    according to jana:
        "3 nts downstream von der pam würde ich erwarten 😃"

        NGG that is found directly downstream of the target sequence in 
        the genomic DNA, on the non-target strand.
        (https://eu.idtdna.com/pages/support/faqs/what-is-a-pam-sequence-and-where-is-it-located)

        PAM is included in the sequences obtained by jana --> no need to scan for that :)
    """

    # get mapping table and reformat
    sgMappingTab = f"{mapping_dir}/CRISPOR_predicted off targts-WRKY30del.xlsx"
    sgMappingTab = pd.read_excel(sgMappingTab)
    newCols = ["relevant_1", "relevant_2", "relevant_3"]
    sgMappingTab[newCols] = sgMappingTab['relevant_for(ID2)'].str.split(',', 3, expand=True)
    idCols = sgMappingTab.columns[~sgMappingTab.columns.isin(newCols)]
    sgMappingTab = sgMappingTab.melt(idCols).dropna()
    sgMappingTab["sample"] = "5603_" + sgMappingTab["value"]
    sgMappingTab = sgMappingTab.loc[sgMappingTab["sample"] == sample,:]
    cols = ["chrom", "start", "end", "strand", "sample"]
    sgMappingTab = sgMappingTab[cols].reset_index(drop=True)
    sgMappingTab["chrom"] = "chr" + sgMappingTab["chrom"].astype(str)


    # load the variant table from syri
    if variantSet == "syri":

        # filtered syri output (centromeric regions removed)
        variantsTab = f"{assemblyDir}/{sample}/syri/{sample}_clean.vcf"

        # use raw syri output
        variantsTab = f"{assemblyDir}/{sample}/syri/{sample}_syri.vcf"

        variantsTab = pd.read_csv(variantsTab, sep="\t", header=None, comment="#")
        variantsTab.columns = ["chrom", "start", "id", "ref", "alt", "qual", "filter", "format"]
        variantsTab["end"] = variantsTab["format"].str.split(";").str[0].str.split("=").str[1]


    # reformat svim-asm predictions
    if variantSet == "svimasm":
        variantsTab = f"{assemblyDir}/{sample}/{sample}_mapped_clean.vcf"
        variantsTab = pd.read_csv(variantsTab, sep="\t")
        variantsTab["start"] = variantsTab["pos"]
        variantsTab["type"] = variantsTab["source"].str.split(".").str[1]
        variantsTab = variantsTab.loc[variantsTab["type"].isin(["DEL", "INS"]),:]
        variantsTab["end"] = variantsTab["format"].str.split(";").str[1].str.split("=").str[1].astype(int)

    # things to consider for off-target predictions:
        # reference assembly used to check for off targets
        # cas9 cuts 3-4 nt unpstream of PAM: https://www.addgene.org/crispr/zhang/faq/#:~:text=The%20Cas9%20cuts%203%2D4bp,target%20DSBs%20using%20wildtype%20Cas9.
        # InDels, or also SNPs --> all of them if feasable
        # go upstream for (-) strand binding sites? --> yes

    offTargetStart = []
    offTargetEnd = []

    for chrom in sgMappingTab["chrom"].unique():
        sgMappingTabChrom = sgMappingTab.loc[sgMappingTab["chrom"] == chrom,:]
        variantsTabChrom = variantsTab.loc[variantsTab["chrom"] == chrom]

        # first, just check for variants that are within the whole sgRNA range
        possibleTargets = [np.array(range(e[0], e[1])) for e in zip(sgMappingTabChrom["start"], sgMappingTabChrom["end"])]
        for target in possibleTargets:
            offTargetStartdf = variantsTabChrom.loc[variantsTabChrom["start"].isin(target),:]
            offTargetEnddf = variantsTabChrom.loc[variantsTabChrom["end"].isin(target),:]

            offTargetStart.append(offTargetStartdf)
            offTargetEnd.append(offTargetEnddf)

    offTargetStart = pd.concat(offTargetStart)
    offTargetEnd = pd.concat(offTargetEnd)


    offTargets = pd.concat([offTargetStart, offTargetEnd])
    offTargets["variantSet"] = variantSet

    return offTargets


def checkOffTargets(mapping_dir, assemblyDir, sampleList):
    """
    function to apply off-target search to all samples
    """
    for sample in sampleList:
            offTargets_syri = checkOffTargets_sample(mapping_dir, sample, assemblyDir, "syri")
            offTargets_svim = checkOffTargets_sample(mapping_dir, sample, assemblyDir, "svimasm")

            offTargets = pd.concat([offTargets_syri, offTargets_svim])
            offTargets_out = f'{assemblyDir}/{sample}_offTargets.tsv'
            offTargets.to_csv(offTargets_out, sep='\t', index=False)


def runSample(sample, assay_dir, nproc, assembler, haplotypeResolved, filterReads, assemblyDir, refFasta, performAsm, mapping_dir):
    """
    function to completly process one sample
    end to end
    """

    # derive paths to data
    raw_dir=f"{assay_dir}/raw_data"
    ref_dir = f'{assay_dir}/reference'
    ref = f"{ref_dir}/{refFasta}"
    ref_gRNA = f"{ref_dir}/GCF_000001735.4_TAIR10.1_genomic_gRNA.fna"


    # print
    print(f"processing sample:  {sample}")
    print(f'raw data from:      {raw_dir}')


    # unzip the data
    unzipLog = unzipData(raw_dir, sample)

    # remove 2kb control
    log2k = remove2kb(raw_dir, sample, nproc)

    # map reads against reference for IGV visualization
    mappingLog, mappedSam, mappedBam = mapReads(assay_dir, filterReads, sample, nproc, ref)

    # perform variant calling using deepVariant
    # deepVariant(assay_dir, filterReads, sample, nproc, ref)


    if performAsm:

        # perform assembly
        assemblyPath, assemblyDir, assemblyLog = assembly(assay_dir, sample, nproc, assembler, haplotypeResolved, filterReads, assemblyDir)

        if not haplotypeResolved:
            # get assembly stats
            statsLog = getAssemblyStats(assemblyPath, sample, assembler, assemblyDir)

            # score assembly using busco
            buscoLog = runBusco(assemblyPath, sample, assemblyDir, nproc)

            # scaffolding using ragtag
            ragtagLog, scaffold = scaffolding(assemblyDir, nproc, assemblyPath, sample, ref)

            # plot synteny
            synLog = plotSyn(scaffold, ref, sample, nproc, assemblyDir, ref_dir)

            # snv calling
            svLog, vcfClean = svCalling(assemblyDir, sample, ref, scaffold, nproc, ref_dir, haplotypeResolved)


        if haplotypeResolved:
            scaffolds = []
            for haplotype in assemblyPath:
                # get assembly stats
                statsLog = getAssemblyStats(haplotype, sample, assembler, assemblyDir)

                # score assembly using busco
                buscoLog = runBusco(haplotype, sample, assemblyDir, nproc)

                # scaffolding using ragtag
                ragtagLog, scaffold = scaffolding(assemblyDir, nproc, haplotype, sample, ref)
                scaffolds.append(scaffold)

                # plot synteny
                synLog = plotSyn(scaffold, ref, sample, nproc, assemblyDir, ref_dir)

            # snv calling
            svLog, vcfClean = svCalling(assemblyDir, sample, ref, scaffolds, nproc, ref_dir, haplotypeResolved)


        # add sample information to cleaned variants
        vcfTab = pd.read_csv(vcfClean, sep='\t', header=None)
        vcfTab.columns = ["chrom", "pos", "source", "ref", "alt", "dot", "filter", "format", "info", "gt"]
        vcfTab["sample"] = sample
        vcfTab.to_csv(vcfClean, sep='\t', index=False)


        # gather the outputs of the logs
        logOut = f"{assemblyDir}/log"
        logs = [unzipLog, log2k, assemblyLog, statsLog, buscoLog, ragtagLog, synLog, svLog]
        steps = ["unzip", "remove2kb", "assembly", "assemblyStats", "busco", "scaffolding", "synteny", "sv-calling"]
        tools = ["-", "blasr", assembler, "-", "busco", "ragtab", "syri", "svim-asm"]
        logDf = pd.DataFrame([steps, tools, logs]).T
        logDf.columns = ["step", "tool", "done"]
        logDf.to_csv(logOut, sep = "\t", index=False)

        return vcfTab




def run():
    """
    main wrapper around the pipeline
    """

    #assay_dir, assembler, nproc = getArgs()

    # debugging args
    assay_dir = "/home/niklas/mount/netscratch/dep_psl/grp_rgo/kieln/collaboration/jana/study_snv_calling_at/assay_snv_calling_at"
    assay_dir = "/netscratch/dep_psl/grp_rgo/kieln/collaboration/jana/study_snv_calling_at/assay_snv_calling_at"
    assembler = "flye" # flye hifiasm
    goi = "WRKY30"
    haplotypeResolved = False
    filterReads = True
    force = True
    nproc = 60
    performAsm = True


    # create reference files
    refFasta = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"  #GCF_000001735.4_TAIR10.1_genomic.fna
    refGTF = "TAIR10_GFF3_genes.gff" # "genomic.gtf"


    # get output name of assay
    reads = "filteredReads" if filterReads else ""
    haplo = "haplo" if haplotypeResolved else ""
    assemblyDir = f'{assay_dir}/assay_{assembler}_{reads}_{haplo}'
    
    if force:
        try:
            print("removing existing dir..")
            shutil.rmtree(assemblyDir)
            createRefDB(assay_dir, refFasta, refGTF)
        except:
            print("dir did not exist")

    # load tables
    mapping_dir = f"{assay_dir}/mapping"
    variants_out = f'{assemblyDir}/variants.tsv'
    mappingTab = pd.read_csv(f'{mapping_dir}/mapping.txt', sep='\t')
    sampleList = list(mappingTab["sample"])

    # run the pipeline
    ret = [runSample(sample, assay_dir, nproc, assembler, haplotypeResolved, filterReads, assemblyDir, refFasta, performAsm, mapping_dir) for sample in sampleList]

    # check for off-targets
    checkOffTargets(mapping_dir, assemblyDir, sampleList)

    # if assembly is performed for snp calling
    if performAsm:
        ret = pd.concat(ret)

        # filter variants and aggregate
        ret = ret.loc[ret["filter"] == "PASS",]
        ret = ret.loc[~ret["ref"].str.contains("N")]
        ret = ret.loc[~ret["alt"].str.contains("N")]

        # drop duplicated variants
        ret = ret.drop_duplicates(subset=["sample", "pos"], keep = False)
        ret["type"] = list(ret["format"].str.split(";").str[0].str.split('=').str[1])
        
        # subset ins and del for separate filtering
        insertions = ret.loc[ret["type"] == "INS"].copy()
        insertions = insertions.drop_duplicates(subset=["sample", "alt"], keep = False)
        deletions = ret.loc[ret["type"] == "DEL"].copy()
        deletions = deletions.drop_duplicates(subset=["sample", "ref"], keep = False)
        ret = pd.concat([deletions, insertions])

        # write final table
        ret.to_csv(variants_out, sep="\t", index=False)

        if goi:
            goiSV = ret.loc[ret["format"].str.contains(goi)]
            goiSV.to_csv(f'{assemblyDir}/{goi}_variants.tsv', sep="\t", index=False)


if __name__ == "__main__":
    run()


