import gzip
import os
import sys


vcffile = snakemake.input.vcf
phenofile  = snakemake.input.pheno
pcfile = snakemake.input.pc 
window_size = int (snakemake.wildcards.winsize)
npc=int(snakemake.params.npc)
overlap =  int(window_size/3)
frq_thresh = float (snakemake.params.frq_thresh)
outfile = snakemake.output.outfile
tmp = '/cluster/work/pausch/naveen/RECOMBINATION/HAPBASED/' +   ("_".join (snakemake.wildcards))


print('removing the old res file if it exists')
if os.path.isfile (outfile): os.remove (outfile)


pc = {}
print ('reding the pc file')
with open (pcfile) as inf:
    for line in inf:
        spl =line.rstrip().split ()
        npc_used=len(spl)-2
        if npc_used < npc:
            print (f'{npc}s are not available in the pc_file, only {npc_used} will be used')
        tw = "\t".join (spl)
        pc [spl [0]] = spl [2:npc+2]



#hash phenotypes
print ('reding the phenotype file')
pheno=open (phenofile, "rt")
info= {}
for line in pheno:
    myid, fid, pheno, breed=line.rstrip().split("\t")
    myid = breed+myid
    if myid in pc:
        mypc = pc [myid]
        tw = "\t".join (mypc + [ pheno ]  )
        info [myid] = tw
    else:
        print (f'pc info not available for {spl [0] }')
    
def prepare_files () :
    startpos =pos [0]
    endpos= pos [-1]
    ##print (f'{nwin}\t{startpos:,}\t{endpos:,}')
    #first get the frequencies and test only haplotypes above a certain freq.
    freq = dict ()
    n = len(ids)*2
    for myhap1,myhap2 in zip (hap1, hap2):
        freq [myhap1] = freq.get (myhap1,0) +1
        freq [myhap2] = freq.get (myhap2,0) +1

    haplotypes=dict()
    for myhap in freq:
        myfreq=freq.get (myhap)/n
        if myfreq > frq_thresh:
            ##print (f'frequency of hap {myhap} at position {startpos} is {myfreq}')
            haplotypes [myhap] = myfreq
            mypos=open(tmp, "w")
            for myid, myhap1, myhap2 in zip (ids, hap1,hap2):
                if myid in info:
                    mygt=str((myhap == myhap1) + (myhap ==myhap2))
                    to_out = "\t".join (   [myid, mygt, info.get (myid)]  )
                    mypos.write (f"{to_out}\n")
            mypos.close ()
            positions = ":".join (list (map (str, pos)  ))
            Rcommand= " ".join ( [   "Rscript  asso.R", str(startpos), str(endpos), str(myhap), str(myfreq), positions,outfile, tmp, str(npc_used)   ])
            os.system (Rcommand)


hap1 = list ()
hap2 = list()
pos =list () 
nvar = 0
nwin = 0

inf=gzip.open (vcffile, "rt")
for line in inf:
    if line [0:2] != "##":
        if line [0:6] == "#CHROM":
            spl=line.rstrip ().split ()
            ids = spl [9:]
            nhap = len (ids) *2
            nani = len (ids)
        else:
            spl=line.rstrip ().split ()
            if "," not in spl[3] and "," not in spl[4]:       
                gts = spl [9:]
                nvar+=1
                ##if nwin > 10 : break
                pos.append (int(spl [1]))
                alleles = [spl[3], spl[4]]
                for i, gt in enumerate (gts):  ## genotype for a particular variant..
                    if nvar == 1:  #this will run only once ! ..end of it the hap1 and hap2 are lists with a base (haplength is now equal to 1)
                        hap1.append( alleles [ int(gt[0]) ]  )
                        hap2.append( alleles [ int(gt[2]) ]  )
                    else:
                        hap1 [i]= hap1 [i] + alleles [ int(gt[0]) ]
                        hap2 [i]= hap2 [i] + alleles [ int(gt[2]) ]

                if nvar==window_size:##here the length of pos and length of hap all should match!
                    ##print (f'Length of lists.. pos, hap1, hap2 are {len(pos)}, {len (hap1[0]) }, {len (hap2[0])  }, {pos}')
                    nwin +=1
                    prepare_files ()
                    
                    nvar = overlap  #already you have window_size snps !
                    hap1 = [hap[-overlap:]  for hap in hap1]
                    hap2 = [hap[-overlap:]  for hap in hap2]
                    pos  = pos [-overlap:]


                        
print ('removing the temp file')
os.remove(tmp)
