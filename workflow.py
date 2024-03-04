### RELATE WORFKLOW ###

from gwf import Workflow
import sys, os, re
from collections import defaultdict
from pathlib import Path
import pandas as pd
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect

gwf = Workflow(defaults={'account': 'ari-intern'})


# directories
out_dir = '/home/ari/ari-intern/people/ari/ariadna-intern/steps'
data_big = '/home/ari/ari-intern/data'
script_dir = '/home/ari/ari-intern/people/ari/ariadna-intern/scripts'
data_dir = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000genome'

# function that modifies file path
def modify_path(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


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
genetic_map_chrX=f'{out_dir}/genetic_map_chrX.tsv'

gwf.target_from_template(f'decode_genetic_maps',
    decode_genetic_maps(decode_hg38_sexavg_per_gen, genetic_map_chrX))


# turn diploid females (XX) into two individual haplotypes (haploid individuals) like males
def female_haploid(haploid_vcf, chrX_filtered_eagle2_phased, phased_haplotypes):
    inputs = [haploid_vcf, chrX_filtered_eagle2_phased]
    outputs = [phased_haplotypes]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    python {haploid_vcf} {chrX_filtered_eagle2_phased} | gzip > {phased_haplotypes}
    '''
    return inputs, outputs, options, spec

haploid_vcf=f'{script_dir}/haploid_vcf.py'
chrX_filtered_eagle2_phased=f'{data_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz'
phased_haplotypes=f'{out_dir}/1000g_phased_haplotypes.vcf.gz'

gwf.target_from_template(f'female_haploid',
    female_haploid(haploid_vcf, chrX_filtered_eagle2_phased, phased_haplotypes))


# construct files with haplotype IDs
def haplotype_id(phased_haplotypes, phased_haplotypes_id):
    inputs = [phased_haplotypes]
    outputs = [phased_haplotypes_id]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    # conda install -c bioconda bcftools
    # conda install openssl   ## to install libcrypto.so.1.0.0 library
    bcftools query -l {phased_haplotypes} > {phased_haplotypes_id}
    '''
    return inputs, outputs, options, spec

phased_haplotypes=f'{out_dir}/1000g_phased_haplotypes.vcf.gz'
phased_haplotypes_id=f'{out_dir}/1000g_phased_haplotypes_ids.txt'

gwf.target_from_template(f'haplotype_id', haplotype_id(phased_haplotypes, phased_haplotypes_id))


# construct populations labels mapping each haplotype to a population
# (group haplotypes according to the population to which the individuals carrying those haplotypes belong)
def pop_labels(make_poplabels, phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index, phased_haplotypes_poplabels):
    inputs = [make_poplabels, phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index]
    outputs = [phased_haplotypes_poplabels]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    python {make_poplabels} {phased_haplotypes_id} {high_coverage_seq_index} {related_high_coverage_seq_index} > {phased_haplotypes_poplabels} 
    '''
    return inputs, outputs, options, spec

make_poplabels=f'{script_dir}/make_poplabels.py'
phased_haplotypes_id=f'{out_dir}/1000g_phased_haplotypes_ids.txt'
high_coverage_seq_index=f'{data_dir}/seq_index/1000G_2504_high_coverage.sequence.index'
related_high_coverage_seq_index=f'{data_dir}/seq_index/1000G_698_related_high_coverage.sequence.index'
phased_haplotypes_poplabels=f'{out_dir}/1000g_phased_haplotypes_poplabels.txt'

gwf.target_from_template(f'pop_labels',
    pop_labels(make_poplabels, phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index, phased_haplotypes_poplabels))


# convert X chromosome VCF for all samples to haps/sample (format required by RELATE)
def convert_vcf(RelateFileFormats, phased_haplotypes_haps, phased_haplotypes_sample, phased_haplotypes, phased_haplotypes_poplabels):
    inputs = [RelateFileFormats, phased_haplotypes_poplabels, phased_haplotypes]
    outputs = [phased_haplotypes_haps, phased_haplotypes_sample]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    {RelateFileFormats} --mode ConvertFromVcf --haps {phased_haplotypes_haps} --sample {phased_haplotypes_sample} -i {phased_haplotypes} --poplabels {phased_haplotypes_poplabels}
    '''
    return inputs, outputs, options, spec

RelateFileFormats='/home/ari/ari-intern/people/ari/relate/bin/RelateFileFormats'
phased_haplotypes_poplabels=f'{out_dir}/1000g_phased_haplotypes_poplabels.txt'
phased_haplotypes=f'{out_dir}/1000g_phased_haplotypes.vcf.gz'
phased_haplotypes_haps=f'{out_dir}/1000g_phased_haplotypes.haps'
phased_haplotypes_sample=f'{out_dir}/1000g_phased_haplotypes.sample'

gwf.target_from_template(f'convert_vcf',
    convert_vcf(RelateFileFormats, phased_haplotypes_haps, phased_haplotypes_sample, phased_haplotypes, phased_haplotypes_poplabels))



## start with specific population ##

# exclude related individuals to avoid biases arising from shared genetic material
def exclude_related(path, population):
    output_dir = f'{out_dir}/{population}/excluded'
    output_path = modify_path(path, parent=output_dir, suffix='_related.txt')
    inputs = {'path' : path}
    outputs = {'path' : output_path}
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v '#' {path} | cut -f 10 > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# find IDs of haplotypes from all other populations so we can exclude them
def other_ppl(path, population):
    output_dir = f'{out_dir}/{population}/excluded'
    output_path = modify_path(path, parent=output_dir, suffix='_non_ppl.txt')
    inputs = {'path' : path}
    outputs = {'path' : output_path}
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v {population} {path} | cut -f 1 -d ' ' > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# combine excluded files: both related and non population individuals
def combine_files(path, related=None):
    output_path = modify_path(path, base='', suffix='excluded_combined.txt')
    out_dir = modify_path(output_path, base='', suffix='')
    inputs = {'path': path, 'related': related}
    outputs = {'path': output_path}
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    mkdir -p {out_dir}
    cat {path} {related} | sort | uniq > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct a list of excluded individuals
def excluded_list(path, haplotype_id=None):
    output_path = modify_path(path, base='', suffix='excluded_list.txt')
    out_dir = modify_path(output_path, base='', suffix='')
    inputs = {'path': path, 'haplotype_id': haplotype_id}
    outputs = {'path': output_path}
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    mkdir -p {out_dir}
    grep -f {path} {haplotype_id} > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct a list of only individuals from the population of interest
def only_ppl(path, poplabels=None):
    output_dir = f'{out_dir}/{population}/included'
    output_path = modify_path(path, parent=output_dir, suffix='_only_ppl.txt')
    inputs = {'path': path, 'poplabels': poplabels}
    outputs = {'path': output_path}
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v -f {path} {poplabels} > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




population = 'LWK'
input_related = [(f'{data_dir}/seq_index/1000G_698_related_high_coverage.sequence.index', population)]
related_target = gwf.map(exclude_related, input_related)
related = related_target.outputs[0]  # list

input_other_ppl = [(f'{out_dir}/1000g_phased_haplotypes_poplabels.txt', population)]
other_ppl_target = gwf.map(other_ppl, input_other_ppl)

combine_target = gwf.map(combine_files, other_ppl_target.outputs, extra = {'related':related})

haplotype_id = f'{out_dir}/1000g_phased_haplotypes_ids.txt'
exclude_list_target = gwf.map(excluded_list, combine_target.outputs, extra = {'haplotype_id':haplotype_id})

poplabels = f'{out_dir}/1000g_phased_haplotypes_poplabels.txt'
include_list = gwf.map(only_ppl, exclude_list_target.outputs, extra = {'poplabels':poplabels})