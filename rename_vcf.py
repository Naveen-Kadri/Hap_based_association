import gzip

infile  = snakemake.input.vcf
outfile = snakemake.output.vcf
breed = snakemake.wildcards.breed
out = open (outfile, "w")

with gzip.open (infile, "rt") as inf:
    for line in inf:
        if line[0:6] == "#CHROM":
            spl = line.rstrip().split()
            for i, el in enumerate (spl):
                if i >=9:
                    spl[i] = breed + spl[i]
            tw = "\t".join (spl)
            out.write (f'{tw}\n')
        else:
            out.write (f'{line}')
            

    
