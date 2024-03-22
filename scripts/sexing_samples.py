import sys
import gzip

_, vcf_file_name = sys.argv

males = []
females = []
first_non_par_pos = None

print('sexing samples', file=sys.stderr)
with gzip.open(vcf_file_name, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
        else:
            try:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            except ValueError as e:
                print(f"Skipping line: {line.strip()}")
                continue
            
            if any('|' not in call for call in calls):
                first_non_par_pos = int(POS)
                print(f'first_non_par_pos: {first_non_par_pos}', file=sys.stderr)
                with open('sexes.txt', 'w') as sex_file:
                    for i, call in enumerate(calls):
                        if '|' in call:
                            females.append(all_samples[i])
                            print(all_samples[i], 'F', file=sex_file)
                        else:
                            males.append(all_samples[i])
                            print(all_samples[i], 'M', file=sex_file)
                assert len(all_samples) == len(males) + len(females)
                print(first_non_par_pos)
                break

