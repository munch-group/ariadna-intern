import sys
import gzip

_, vcf_file_name = sys.argv

males = []
females = []

print('sexing samples', file=sys.stderr)

def open_file(file_name, mode):
    try:
        with gzip.open(file_name, mode) as f:
            return f.readlines()
    except OSError:
        with open(file_name, mode) as f:
            return f.readlines()

lines = open_file(vcf_file_name, 'rt')

first_non_par_pos = float('inf')  # Initialize with a large value

for line in lines:
    if line.startswith('#'):
        if line.startswith('#CHROM'):
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
    else:
        fields = line.split()
        if len(fields) < 9:
            print(f"Skipping line with insufficient fields: {line}", file=sys.stderr)
            continue

        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = fields

        if any('|' not in call for call in calls):
            first_non_par_pos = int(POS)
            print('first non-par pos', first_non_par_pos, file=sys.stderr)
            with open('sexes.txt', 'w') as sex_file:
                for i, call in enumerate(calls):
                    if '|' in call:
                        females.append(all_samples[i])
                        print(all_samples[i], 'F', file=sex_file)
                    else:
                        males.append(all_samples[i])
                        print(all_samples[i], 'M', file=sex_file)
            assert len(all_samples) == len(males) + len(females)
            break

print('writing haploids', file=sys.stderr)

for line in lines:
    if line.startswith('#'):
        if line.startswith('#CHROM'):
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
            ids = []
            for sample in all_samples:
                if sample in males:
                    ids.append(sample)
                else:
                    ids.append(sample + "_1")
                    ids.append(sample + "_2")

            start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
            start_line.extend(ids)
            print('\t'.join(start_line))
        else:
            print(line, end='')
    else:
        fields = line.split()
        if len(fields) < 9 or int(fields[1]) < first_non_par_pos:
            continue
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = fields
        haploid_calls = []
        for call in calls:
            for x in call.split('|'):
                haploid_calls.append(x)
        if len(haploid_calls) > 2 * len(females) + len(males):
            print('aborting at par2:', POS, file=sys.stderr)
            break
        start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
        start_line.extend(haploid_calls)
        print('\t'.join(start_line))


