import sys
import gzip

_, vcf_file_name, first_non_par_pos = sys.argv

print('writing haploids', file=sys.stderr)
with gzip.open(vcf_file_name, 'rt') as f:
    for line in f:
        if not line.startswith('#'):
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            if int(POS) < int(first_non_par_pos):
                continue
            haploid_calls = []
            for call in calls:
                for x in call.split('|'):
                    haploid_calls.append(x)
            start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
            start_line.extend(haploid_calls)
            print('\t'.join(start_line))

