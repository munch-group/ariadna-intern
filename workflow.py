### RELATE WORFKLOW ###

from gwf import Workflow

gwf = Workflow(defaults={'account': 'ari-intern'})
import sys, os, re

# directories
output_dir = '/home/ari/ari-intern/people/ari/ariadna-intern/steps'
data_big = '/home/ari/ari-intern/data'
script_dir = '/home/ari/ari-intern/people/ari/ariadna-intern/scripts'
data_dir = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000genome'


# map of recombination rate across the X chromosome made by DECODE genetics
def decode_genetic_maps(decode_hg38_sexavg_per_gen, genetic_map_chrX):
    inputs = [decode_hg38_sexavg_per_gen]
    outputs = [genetic_map_chrX]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    cat {decode_hg38_sexavg_per_gen} | tail -n +2 | grep chrX | cut -f 2,4,5 | (echo pos COMBINED_rate Genetic_Map ; cat - ; ) > {genetic_map_chrX}
    '''
    return inputs, outputs, options, spec

decode_hg38_sexavg_per_gen=f'{data_big}/decode_hg38_sexavg_per_gen.tsv'
genetic_map_chrX=f'{output_dir}/genetic_map_chrX.tsv'

gwf.target_from_template(f'decode_genetic_maps',
    decode_genetic_maps(decode_hg38_sexavg_per_gen, genetic_map_chrX))


# turn diploid females (XX) into two individual haplotypes (haploid individuals) like males
def female_haploid(haploid_vcf, chrX_filtered_eagle2_phased, phased_haplotypes_1000g):
    inputs = [haploid_vcf, chrX_filtered_eagle2_phased]
    outputs = [phased_haplotypes_1000g]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    python {haploid_vcf} {chrX_filtered_eagle2_phased} | gzip > {phased_haplotypes_1000g}
    '''
    return inputs, outputs, options, spec

haploid_vcf=f'{script_dir}/haploid_vcf.py'
chrX_filtered_eagle2_phased=f'{data_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX_filtered_eagle2-phased_v2.vcf.gz'
phased_haplotypes_1000g=f'{output_dir}/1000g_phased_haplotypes.vcf.gz'

gwf.target_from_template(f'female_haploid',
    female_haploid(haploid_vcf, chrX_filtered_eagle2_phased, phased_haplotypes_1000g))


# construct files with haplotype IDs
def haplotype_id(phased_haplotypes_1000g, phased_haplotypes_1000g_id):
    inputs = [phased_haplotypes_1000g]
    outputs = [phased_haplotypes_1000g_id]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    bcftools query -l {phased_haplotypes_1000g} > {phased_haplotypes_1000g_id}
    '''
    return inputs, outputs, options, spec

phased_haplotypes_1000g=f'{output_dir}/1000g_phased_haplotypes.vcf.gz'
phased_haplotypes_1000g_id=f'{output_dir}/1000g_phased_haplotypes_ids.txt'

gwf.target_from_template(f'haplotype_id',
    haplotype_id(phased_haplotypes_1000g, phased_haplotypes_1000g_id))


# construct populations labels mapping each haplotype to a population
def pop_labels(make_poplabels, phased_haplotypes_1000g_id, high_coverage_seq_index, related_high_coverage_seq_index, phased_haplotypes_1000g_poplabels):
    inputs = [make_poplabels, phased_haplotypes_1000g_id, high_coverage_seq_index, related_high_coverage_seq_index]
    outputs = [phased_haplotypes_1000g_poplabels]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    python {make_poplabels} {phased_haplotypes_1000g_id} {high_coverage_seq_index} {related_high_coverage_seq_index} > {phased_haplotypes_1000g_poplabels}
    '''
    return inputs, outputs, options, spec

make_poplabels=f'{script_dir}/make_poplabels.py'
phased_haplotypes_1000g_id=f'{output_dir}/1000g_phased_haplotypes_ids.txt'
high_coverage_seq_index=f'{data_dir}/1000G_2504_high_coverage.sequence.index'
related_high_coverage_seq_index=f'{data_dir}/1000G_698_related_high_coverage.sequence.index'
phased_haplotypes_1000g_poplabels=f'{output_dir}/1000g_phased_haplotypes_poplabels.txt'

gwf.target_from_template(f'pop_labels',
    pop_labels(make_poplabels, phased_haplotypes_1000g_id, high_coverage_seq_index, related_high_coverage_seq_index, phased_haplotypes_1000g_poplabels))


# convert X chromosome VCF for all samples to haps/sample format (format required by RELATE)
def convert_vcf(RelateFileFormats, phased_haplotypes_1000_haps, phased_haplotypes_1000_sample, phased_haplotypes_1000g_poplabels):
    inputs = [RelateFileFormats, phased_haplotypes_1000g_poplabels]
    outputs = [phased_haplotypes_1000_haps, phased_haplotypes_1000_sample]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    {RelateFileFormats} --mode ConvertFromVcf --haps {phased_haplotypes_1000_haps} --sample {phased_haplotypes_1000_sample} -i 1000g_phased_haplotypes --poplabels {phased_haplotypes_1000g_poplabels}
    '''
    return inputs, outputs, options, spec

RelateFileFormats=f'{data_big}/RelateFileFormats'
phased_haplotypes_1000g_poplabels=f'{output_dir}/1000g_phased_haplotypes_poplabels.txt'
phased_haplotypes_1000_haps=f'{output_dir}/1000g_phased_haplotypes.haps'
phased_haplotypes_1000_sample=f'{output_dir}/1000g_phased_haplotypes.sample'

gwf.target_from_template(f'convert_vcf',
    convert_vcf(RelateFileFormats, phased_haplotypes_1000_haps, phased_haplotypes_1000_sample, phased_haplotypes_1000g_poplabels))


# exclude related individuals to avoid biases
def exclude_related(related_high_coverage_seq_index, related_ids_1000g):
    inputs = [related_high_coverage_seq_index]
    outputs = [related_ids_1000g]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    grep -v '#' {related_high_coverage_seq_index} | cut -f 10 > {related_ids_1000g}
    '''
    return inputs, outputs, options, spec

related_high_coverage_seq_index=f'{data_dir}/1000G_698_related_high_coverage.sequence.index'
related_ids_1000g=f'{output_dir}/1000g_related_ids.txt'

gwf.target_from_template(f'exclude_related',
    exclude_related(related_high_coverage_seq_index, related_ids_1000g))


# analyze only individuals from the African LWK population --> find IDs of haplotypes from all other populations so we can exclude them
