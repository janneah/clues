
import allel
import pandas as pd
import sys

vcffile = allel.read_vcf(sys.argv[1])
sitesfile = open(sys.argv[2], 'w')

df = pd.DataFrame(vcffile['calldata/GT'][:, :, 0],
            columns=vcffile['samples'],
            index=vcffile['variants/POS'])

def get_bases(sr, ref, alt):
    return [x and a or r for x, r, a in zip(sr, ref, alt)]

base_table = df.apply(get_bases, 
                      args=(vcffile['variants/REF'],
                            vcffile['variants/ALT'][:, 0]), axis=0)

#bases_region = base_table.loc[(base_table.index > 100000) & (base_table.index < 110000), ['HG00096', 'HG00097', 'HG00099', 'HG00100']]

#First line
sample_name_str = "\t".join(base_table.columns.values)
title = "NAMES"
print(title + "\t" + sample_name_str, file=sitesfile)

#Second line
title = "REGION"
chrom = str(vcffile['variants/CHROM'][0])
start_position = str(min(vcffile['variants/POS']))
end_position = str(max(vcffile['variants/POS']))
print(title + "\t" + chrom + "\t" + start_position + "\t" + end_position, file=sitesfile)

#The rest
for index, row in base_table.iterrows():
    print(str(index) + "\t" + row.apply(str).str.cat(sep=''), file=sitesfile)

sitesfile.close()