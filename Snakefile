#!/usr/bin/env snakemake

"""
Snakemake Rules
================

Rules for Snakemake.
"""

__version__ = '0.1.0'


rule all_general:
    input:
        expand(
            'reference_data/{build}-all_samples/unprocessed-chr{chr}-n_entries.txt',
            build=config['chromosome_builds'],
            chr=set(list(range(1, 23)) + ['X', 'Y', 'MT'])
        ),
        expand(
            '{build}/{pop}/chr{chr}.bed',
            build=config['chromosome_builds'],
            pop=config['super_populations'],
            chr=set(list(range(1, 23)) + ['X', 'Y', 'MT'])
        ),
        # expand(
        #     '{build}/{pop}/autosome.bed',
        #     build=config['chromosome_builds'],
        #     pop=config['super_populations']
        # )


rule count_entries_unprocessed_vcf:
    """
    Count the number of entires in an unprocessed vcf file
    """
    input:
        vcf='reference_data/{build}-all_samples/unprocessed-chr{chr}.vcf.gz'
    output:
        'reference_data/{build}-all_samples/unprocessed-chr{chr}-n_entries.txt'
    shell:
        'OUT_DIR=$(dirname {input.vcf}); '
	'gunzip -c {input.vcf} | '
	    'grep -v \'#\' | '
	    'wc -l > '
	    '$OUT_DIR/$(basename {input.vcf} .vcf.gz)-n_entries.txt; '


rule vcf_add_rsid_bcftools:
    """
    Replaces existing ID with dbSNP rsid using bcftools.
    In where the ID is missing, it sill be set to '%CHROM\_%POS\_%REF\_%FIRST_ALT'

    Information on dbSNP can be found here:
    https://www.ncbi.nlm.nih.gov/dbvar/content/org_summary/

        For instance, the latest dpSNP vcf at time of writing is here:
        ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/

    If there is an error with the vcf header, it can be edited like so:
    bcftools view -h genotypes-old.vcf.gz > hdr.txt
        edit hdr.txt
            (e.g, ##FILTER=<ID=GENOTYPED,Description="Site was genotyped">)
    bcftools reheader --header hdr.txt genotypes-old.vcf.gz > genotypes.vcf.gz
    """
    input:
        dbsnp='reference_data/{build}-dbsnp.vcf.gz',
        dbsnp_tbi='reference_data/{build}-dbsnp.vcf.gz.tbi',
        file_2_annotate='reference_data/{build}-all_samples/unprocessed-chr{chr}.vcf.gz',
        file_2_annotate_tbi='reference_data/{build}-all_samples/unprocessed-chr{chr}.vcf.gz.tbi'
    output:
        'reference_data/{build}-all_samples/processed-chr{chr}.vcf.gz',
        'reference_data/{build}-all_samples/processed-chr{chr}.vcf.gz.tbi',
        'reference_data/{build}-all_samples/processed-chr{chr}-n_entries.txt'
    shell:
        'OUT_FILE=$(echo {input.file_2_annotate} | sed s/unprocessed/processed/); '
        'OUT_DIR=$(dirname {input.file_2_annotate}); '
        'bcftools annotate '
            '--annotations {input.dbsnp} '
            '--columns ID '
	    '{input.file_2_annotate} | '
        'bcftools annotate '
            '--set-id +\'%CHROM\_%POS\_%REF\_%FIRST_ALT\' '
            '--output $OUT_FILE '
            '--output-type z '
            '-; '
        'tabix -f $OUT_FILE; '
        'gunzip -c $OUT_FILE | '
            'grep -v \'#\' | '
            'wc -l > $OUT_DIR/$(basename $OUT_FILE .vcf.gz)-n_entries.txt; '


rule make_sample_subset:
    """
    Makes a list of sample ids for a specific population.
    Also sub-filters this list to remove closely related participants.
    """
    input:
        sample_info='reference_data/sample_info.tsv',
	related_info='reference_data/sample_info-related.tsv'
    output:
        '{build}/{pop}/samples.txt',
	'{build}/{pop}/samples-unrelated.txt'
    shell:
        'mkdir -p {wildcards.build}/{wildcards.pop}; '
        'cat {input.sample_info} | '
            'grep {wildcards.pop} | '
            'cut -f 1 > {wildcards.build}/{wildcards.pop}/samples.txt; '
        'cat {input.related_info} | '
            'awk -F " " \'NR>1{{print $1}}\' | '
            'grep -v -x -f - {wildcards.build}/{wildcards.pop}/samples.txt > '
            '{wildcards.build}/{wildcards.pop}/samples-unrelated.txt; '


rule subset_vcf:
    """
    Subsets vcf file to a specific sample set.

    WARNING: --force-samples is used which means an error will not be raised
             if <samples_to_keep> contains a sample that is missing from
            <file_2_subset>. Such an instance occurs for EUR:chrY.
    """
    input:
        samples_to_keep='{build}/{pop}/samples-unrelated.txt',
        file_2_subset='reference_data/{build}-all_samples/processed-chr{chr}.vcf.gz'
    output:
        '{build}/{pop}/chr{chr}.vcf.gz',
        '{build}/{pop}/chr{chr}.vcf.gz.tbi',
        '{build}/{pop}/chr{chr}-n_entries.txt'
    shell:
        'mkdir -p {wildcards.build}/{wildcards.pop}; '
        'bcftools view '
            '--samples-file {input.samples_to_keep} '
            '--force-samples '
            '--output-file {wildcards.build}/{wildcards.pop}/chr{wildcards.chr}.vcf.gz '
            '--output-type z '
            '{input.file_2_subset}; '
        'tabix -f {wildcards.build}/{wildcards.pop}/chr{wildcards.chr}.vcf.gz; '
        'gunzip -c {wildcards.build}/{wildcards.pop}/chr{wildcards.chr}.vcf.gz | '
            'grep -v \'#\' | '
            'wc -l > {wildcards.build}/{wildcards.pop}/chr{wildcards.chr}-n_entries.txt; '


rule combine_vcf:
    """
    Combines vcf files.
    """
    input:
        expand(
            '{{build}}/{{pop}}/chr{chr}.vcf.gz',
            chr=range(1, 23)
        )
    output:
        '{build}/{pop}/autosome.vcf.gz',
        '{build}/{pop}/autosome.vcf.gz.tbi',
        '{build}/{pop}/autosome-n_entries.txt'
    shell:
        'bcftools concat '
            '--ligate '
            '--output {wildcards.build}/{wildcards.pop}/autosome.vcf.gz '
            '--output-type z '
            '{input}; '
        'tabix -f {wildcards.build}/{wildcards.pop}/autosome.vcf.gz; '
        'gunzip -c {wildcards.build}/{wildcards.pop}/autosome.vcf.gz | '
            'grep -v \'#\' | '
            'wc -l > {wildcards.build}/{wildcards.pop}/autosome-n_entries.txt'


rule vcf2plink:
    """
    Converts vcf file to plink.
    """
    input:
        vcf='{build}/{pop}/{file}.vcf.gz'
    output:
        '{build}/{pop}/{file}.bed',
        '{build}/{pop}/{file}.bim',
        '{build}/{pop}/{file}.fam'
    shell:
        'OUT_DIR=$(dirname {input.vcf}); '
        'plink2 '
            '--vcf {input.vcf} '
            '--keep-allele-order '
            '--make-bed '
            '--const-fid 0 '
            '--out $OUT_DIR/$(basename {input.vcf} .vcf.gz); '
        'rm $OUT_DIR/$(basename {input.vcf} .vcf.gz).nosex; '
