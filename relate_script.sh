#!/bin/bash
mkdir -p ./steps/LWK/relate/run_relate
cd ./steps/LWK/relate/run_relate
/home/ari/ari-intern/people/ari/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample /home/ari/ari-intern/people/ari/ariadna-intern/steps/LWK/relate/1000g_ppl_phased_haplotypes.sample.gz --haps /home/ari/ari-intern/people/ari/ariadna-intern/steps/LWK/relate/1000g_ppl_phased_haplotypes.haps.gz --map . --annot /home/ari/ari-intern/people/ari/ariadna-intern/steps/LWK/relate/1000g_ppl_phased_haplotypes.annot --dist /home/ari/ari-intern/people/ari/ariadna-intern/steps/LWK/relate/1000g_ppl_phased_haplotypes.dist.gz --memory 20 -o ./steps/LWK/relate/run_relate/1000g_ppl_phased_haplotypes
