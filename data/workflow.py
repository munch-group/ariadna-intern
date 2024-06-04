
### WORKFLOW TO DOWNLOAD ALL DATA ###

from gwf import Workflow
gwf = Workflow(defaults={'account': 'ari-intern'})

def download_file(url, output_file, output_dir):
    return f"""
    mkdir -p {output_dir}  # create dir if not existent
    wget {url} -O {output_dir}/{output_file}
    """

files_to_download = [
    # ancestral sequence: human genome as it were at the common ancestor of all humans
    {'url': 'http://ftp.ensembl.org/pub/release-109/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz',
    'output_file': 'homo_sapiens_ancestor_GRCh38.tar.gz'},
    
    # file with SNPs for individuals from the 1000 genomes project
    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz',
    'output_file': 'CCDG_14151_B01_GRM_WGS_2020-08-05_chrX_filtered_eagle2-phased_v2.vcf.gz'},
    
    # index for Variant Call Format (VCF) file for more efficient retrieval of information from a large dataset
    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08_05_chrX_filtered.eagle2-phased.v2.vcf.gz.tbi',
    'output_file': 'CCDG_14151_B01_GRM_WGS_2020-08-05_chrX_filtered.eagle2-phased.v2.vcf.gz.tbi'},

    # documentation files for the phases version of the 1000 genomes data set
    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv',
    'output_file': 'phased-manifest_July2021.tsv'},

    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/README_SNV_INDEL_phasing_111822.pdf',
    'output_file': 'README_SNV_INDEL_phasing_111822.pdf'},
    
    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/1000G_NYGC_phasing_CHANGE_LOG.pdf',
    'output_file': '1000G_NYGC_phasing_CHANGE_LOG.pdf'},

    # sequence index files
    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index',
    'output_file': '1000G_2504_high_coverage.sequence.index'},

    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index',
    'output_file': '1000G_698_related_high_coverage.sequence.index'},

    # strict mask of unreliably called bases
    {'url': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.chrX.mask.fasta.gz',
    'output_file': '20160622.chrX.mask.fasta.gz'},
]

data_dir = 'steps/1000genome'  # specify output directory
options = {'memory': '1g', 'walltime': '00:10:00'}   # specify resource options for each target

for file in files_to_download:
    # replacements done to ensure that the target names are valid
    target_name = f"download_{file['output_file'].replace('.', '_').replace('-', '_').replace(' ', '_').replace(',', '_')}"
    
    # extracting the output file name from the file information dictionary
    output_file = file['output_file']

    # create target for each file to be downloaded
    gwf.target(target_name, inputs=[], outputs=[f"{data_dir}/{output_file}"]) << download_file(file['url'], output_file, data_dir)

if __name__ == '__main__':
    # generate a list of target names
    target_names = [f"download_{file['output_file'].replace('.', '_').replace('-', '_').replace(' ', '_').replace(',', '_')}" for file in files_to_download]

    # create a main target than depends on all download targets
    gwf.target('main', target_names)


