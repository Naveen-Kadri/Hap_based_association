configfile:"config.yaml"

#PROGRAMS
plink = config ['plink']
GCTA =config ['gcta']


#WILDCARDS
chromosomes = range (1,30)
beagle = config ['beagle']
sexes =['male', 'female']
breeds = ['bv', 'fv', 'all']
winsizes = [12,15]
frq_threshs = [0.01,0.05,0.1]


#OUT/INPUT FILES
OUT_DIR = config ['OUT_DIR']
vcf=OUT_DIR + '/REFALT/{breed}/CHR{chr}/clean.vcf.gz'
grr_file=OUT_DIR + "/REFALT/LINKPHASE_RUN2/{breed}/cleaned_GRR.txt"
good_gams=OUT_DIR + '/REFALT/LINKPHASE_RUN2/{breed}/good_gams.txt'


rule all:
    input:
        expand(OUT_DIR + "/REFALT/HAPBASED/{winsize}/{breed}/{sex}/result_{frqthresh}.txt", winsize=winsizes, frqthresh=frq_threshs, breed=breeds,sex=sexes)

        
rule make_phenotypes:
    input:
        grr_file=grr_file,
        good_gams=good_gams
    output:
        outfile_raw=OUT_DIR + '/REFALT/HAPBASED/PHENOTYPES/{breed}/{sex}/phenotypes_raw.txt',
        outfile_mean=OUT_DIR + '/REFALT/HAPBASED/PHENOTYPES/{breed}/{sex}/phenotypes_mean.txt'
    script:
        "make_phenotypes.R"


rule common_markers:
    input:
        infiles = expand ('/cluster/work/pausch/naveen/RECOMBINATION/REFALT/{breed}/CHR{{chr}}/clean.map', breed=breeds)
    output:
        outfile=OUT_DIR + "/REFALT/HAPBASED/COMMON/merged/CHR{chr}/common_markers.txt"
    script:
        "common_markers.py"


rule combine_phenotypes:
    input:
        infiles =expand(OUT_DIR + '/REFALT/HAPBASED/PHENOTYPES/{breed}/{{sex}}/phenotypes_{{type}}.txt',  breed=breeds)
    output:
        outfile=OUT_DIR + '/REFALT/HAPBASED/PHENOTYPES/all/{sex}/phenotypes_{type}.txt'
    shell:
        "cat {input.infiles} >{output.outfile}"
        
rule make_vcf:
    input:
        vcf=vcf,
        phenotyped = OUT_DIR + '/REFALT/HAPBASED/PHENOTYPES/merged/{sex}/phenotypes_mean.txt',
        common_markers = rules.common_markers.output.outfile
    output:
        outfile=temp( OUT_DIR + "/REFALT/HAPBASED/COMMON/{breed}/{sex}/CHR{chr}/cleaned.vcf.gz") 
    params:
        out_prefix = lambda wildcarads, output : output.outfile [ :-7]
    resources:
        mem_mb = 16000,
        walltime="00:30"
    shell:
        '''
        module load plink;
        {plink} --cow
        --vcf  {input.vcf} 
        --keep-allele-order
        --keep {input.phenotyped}
        --extract {input.common_markers}
        --recode vcf-fid bgz 
        --out {params.out_prefix}
        '''

rule rename_vcf:
    ''' Add breed info to the id '''
    input:
        vcf=rules.make_vcf.output.outfile
    output:
        vcf=OUT_DIR + "/REFALT/HAPBASED/COMMON/{breed}/{sex}/CHR{chr}/cleaned_renamed.vcf"
    script:
        "rename_vcf.py"

rule zip_index:
    input:
        vcf=rules.rename_vcf.output.vcf
    output:
        vcf=OUT_DIR + "/REFALT/HAPBASED/COMMON/{breed}/{sex}/CHR{chr}/cleaned_renamed.vcf.gz",
        index=OUT_DIR + "/REFALT/HAPBASED/COMMON/{breed}/{sex}/CHR{chr}/cleaned_renamed.vcf.gz.tbi"
    shell:
        '''
        module load htslib
        bgzip {input.vcf}
        tabix -p vcf {input.vcf}.gz
        '''

rule merge_vcfs:
    input:
        vcfs = expand (OUT_DIR + "/REFALT/HAPBASED/COMMON/{breed}/{{sex}}/CHR{{chr}}/cleaned_renamed.vcf.gz", breed=breeds),
    output:
        vcf = OUT_DIR + "/REFALT/HAPBASED/COMMON/merged/{sex}/CHR{chr}/merged.vcf.gz"
    resources:
        mem_mb=8000,
        walltime="01:00"
    threads:
        4
    shell:
        '''
        module load gcc/8.2.0 bcftools/1.6;
        bcftools merge {input.vcfs} \
        --output {output.vcf} \
        --output-type z\
        --threads {threads}
        '''
        
rule phase_merged:
    ''' For a combined breed hapbased association analysis.. all samples should be phased together'''
    input:
        vcf=rules.merge_vcfs.output.vcf
    output:
        vcf = OUT_DIR + "/REFALT/HAPBASED/COMMON/merged/{sex}/CHR{chr}/phased_cleaned.vcf.gz"
    params:
        window=60,
        overlap=10,
        outfile=lambda wildcards,output : output.vcf [:-7],
        tempdir = '/cluster/work/pausch/temp_scratch/naveen/',
        mem_g = lambda wildcards, resources, threads : int ( (resources.mem_mb/1000) * threads )
    resources:
        mem_mb=8000,
        walltime='24:00'
    threads:
        12
    shell:
        '''
        module load jdk \n"
        java -Djava.io.tmpdir={params.tempdir} -Xss5m -Xmx{params.mem_g}g -jar {beagle} \
        gt={input.vcf} \
        out={params.outfile}  \
        window={params.window} \
        overlap={params.overlap} \
        nthreads={threads}
        '''
        
rule make_bed:
    input:
        vcf = rules.merge_vcfs.output.vcf
    output:
        binaries = expand ( OUT_DIR + "/HAPBASED/COMMON/merged/{{sex}}/CHR{{chr}}/cleaned.{ext}", ext =['bim', 'fam', 'bed'] )
    params:
        out_prefix = lambda wildcards, output : output.binaries[0][:-4]
    resources:
        mem_mb=16000,
        walltime="00:30"
    shell:
        '''
        module load plink;
        {plink} --cow \
        --vcf {input.vcf}  \
        --keep-allele-order \
        --make-bed \
        --out {params.out_prefix}
        '''

rule make_grm:
    input:
        infiles=rules.make_bed.output.binaries
    output:
        outfiles=expand (OUT_DIR + '/REFALT/GRM/{{sex}}/CHR{{chr}}/chrwise.{ext}', ext=['grm.id', 'grm.bin', 'grm.N.bin'])
    params:
        in_prefix = lambda wildcards, input : input.infiles [0] [:-4],
        out_prefix = lambda wildcards, output : output.outfiles [0][ :-7],
        maf_thresh = 0.05
    resources:
        mem_mb=8000,
        walltime="01:00"
    threads:
        10
    shell:
        '''
        {GCTA} --autosome-num 29 \
        --make-grm-bin \
        --bfile {params.in_prefix}\
        --out {params.out_prefix}\
        --maf {params.maf_thresh} \
        --threads {threads}"
        '''

rule make_grm_list:
    input:
        infiles = expand (OUT_DIR + '/REFALT/GRM/{{sex}}/CHR{chr}/chrwise.{ext}', ext=['grm.id', 'grm.bin', 'grm.N.bin'], chr=chromosomes)
    output:
        outfile=OUT_DIR + "/REFALT/GRM/{sex}/grm_list.txt"
    script:
        "make_grm_list.py"

rule merge_grm:
    input:
        infiles = expand (OUT_DIR + '/REFALT/GRM/{{sex}}/CHR{chr}/chrwise.{ext}', chr=chromosomes, ext=['grm.id', 'grm.bin', 'grm.N.bin'] ),
        grm_list=OUT_DIR + "/REFALT/GRM/{sex}/grm_list.txt"
    output:
        outfiles=expand (OUT_DIR + '/REFALT/GRM/{{sex}}/grm.{ext}', ext=['grm.id', 'grm.bin', 'grm.N.bin'] )
    resources:
        mem_mb=64000,
        walltime="10:00"
    threads:
        4
    params:
        out_prefix = lambda wildcards, output : output.outfiles [0] [ :-7]
    shell:
        '''
        {GCTA} --autosome-num 29 \
         --mgrm {input.grm_list} \
         --make-grm-bin \
         --out {params.out_prefix} \
         --threads {threads}
        '''

rule pca:
    input:
        infiles=expand (OUT_DIR + '/REFALT/GRM/{{sex}}/grm.{ext}', ext=['grm.id', 'grm.bin', 'grm.N.bin'] )
    output:
        outfiles=expand (OUT_DIR + '/REFALT/GRM/{{sex}}/genome.{ext}', ext=['eigenval', 'eigenvec'] )
    params:
        in_prefix = lambda wildcards, input : input.infiles [0] [:-7],
        out_prefix = lambda wildcards, output : output.outfiles [0] [:-9],
        num_pcs=20
    resources:
        mem_mb=64000,
        walltime="24:00"
    threads:
        4
    shell:
        '''
        {GCTA} --autosome-num 29
        --pca {params.num_pcs}
        --grm {params.in_prefix}
        --out {params.out_prefix}
        --threads {threads}
        '''

rule plot_pca:
    input:
        infile=OUT_DIR + '/REFALT/GRM/{sex}/genome.eigenvec' 
    output:
        plotfile=OUT_DIR + '/REFALT/GRM/{sex}/pca_{sex}.pdf' 
    script:
        "plot_pca.R"
        
rule association:
    input:
        vcf=rules.phase_merged.output.vcf,
        pheno=OUT_DIR + '/REFALT/HAPBASED/PHENOTYPES/{breed}/{sex}/phenotypes_mean.txt',
        pc=OUT_DIR + '/REFALT/GRM/{sex}/genome.eigenvec'
    output:
        outfile=OUT_DIR + "/REFALT/HAPBASED/{winsize}/{breed}/{sex}/CHR{chr}/result.txt"
    params:
        frq_thresh=0.01,
        npc=10
    resources:
        mem_mb=8000,
        walltime="04:00"
    script:
        "gwas.py"
        
rule merge_results:
    input:
        infiles = expand (OUT_DIR + "/REFALT/HAPBASED/{{winsize}}/{{breed}}/{{sex}}/CHR{chr}/result.txt",chr=chromosomes)
    output:
        outfile=OUT_DIR + "/REFALT/HAPBASED/{winsize}/{breed}/{sex}/result_{frqthresh}.txt"
    script:
        "merge_results.R"

rule manhattan:
    input:
        infile=rules.merge_results.output.outfile
    output:
        plotfile=OUT_DIR + "/REFALT/HAPBASED/{winsize}/{breed}/{sex}/manhattan_{breed}_{sex}.pdf"
    script:
        "manhattan.R"
        
rule all_manhattan:
    input:
        infiles=expand(OUT_DIR + "/REFALT/HAPBASED/{{winsize}}/{breed}/{sex}/result.txt",sex=sexes, breed=['bv','fv','all'])
    output:
        plotfile=OUT_DIR + "/REFALT/HAPBASED/{winsize}/all_manhattan_{winsize}.pdf"
    params:
        breeds = breeds,
        sexes = sexes
    script:
        "all_manhattan.R"
        
        
