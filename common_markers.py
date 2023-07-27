infiles = snakemake.input.infiles
outfile = snakemake.output.outfile
out = open (outfile, 'w')


seen  = {}
n1=0
n2=0
ncommon=0
with open (infiles [0]) as inf:
    for line in inf:
        n1+=1
        spl=line.rstrip().split()
        seen [spl [1]] =1

with open (infiles [1]) as inf:
    for line in inf:
        n2+=1
        spl=line.rstrip().split()
        if spl[1] in seen:
            ncommon+=1
            out.write (f'{spl[1]}\n')
out.close()

print (f'Number of markers in first, second and common file : {n1}\t{n2}\t{ncommon}')
