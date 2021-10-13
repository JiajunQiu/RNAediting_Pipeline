# RNAediting_Pipeline
# Dependency
Install the following softwares and packages before you run the pipeline
bamUtil: https://github.com/statgen/bamUtil
gatk3: conda install gatk
samtools: http://www.htslib.org/
bedtools: https://bedtools.readthedocs.io/en/latest/
scipy: https://www.scipy.org/
statsmodels: https://www.statsmodels.org/stable/index.html

# Usage
python3 RNAediting_Pipeline -i your_bam_file -g your_human_genome_fa_file -s your_vcf_file_for_known_SNPs
