# Reference data staging

Download and unpack

> All reference data is retrieved here, exclude unused files as desired; using GRCh38_hmf below

```bash
base_url=https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/genomes

fps='
genomes/GRCh37_hmf/Homo_sapiens.GRCh37.GATK.illumina.fasta
genomes/GRCh37_hmf/bwa_index/0.7.17-r1188.tar.gz
genomes/GRCh37_hmf/bwa_index/2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta.0123
genomes/GRCh37_hmf/bwa_index/2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta.bwt.2bit.64
genomes/GRCh37_hmf/bwa_index_image/0.7.17-r1188/Homo_sapiens.GRCh37.GATK.illumina.fasta.img
genomes/GRCh37_hmf/gridss_index/2.13.2/Homo_sapiens.GRCh37.GATK.illumina.fasta.gridsscache
genomes/GRCh37_hmf/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict
genomes/GRCh37_hmf/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai
genomes/GRCh37_hmf/star_index/gencode_19/2.7.3a.tar.gz
genomes/GRCh38_hmf/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
genomes/GRCh38_hmf/bwa_index/0.7.17-r1188.tar.gz
genomes/GRCh38_hmf/bwa_index/2.2.1/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.0123
genomes/GRCh38_hmf/bwa_index/2.2.1/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwt.2bit.64
genomes/GRCh38_hmf/bwa_index_image/0.7.17-r1188/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img
genomes/GRCh38_hmf/gridss_index/2.13.2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gridsscache
genomes/GRCh38_hmf/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict
genomes/GRCh38_hmf/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
genomes/GRCh38_hmf/star_index/gencode_38/2.7.3a.tar.gz
hmf_reference_data/hmftools/5.34_37--2.tar.gz
hmf_reference_data/hmftools/5.34_38--2.tar.gz
hmf_reference_data/panels/tso500_5.34_37--1.tar.gz
hmf_reference_data/panels/tso500_5.34_38--1.tar.gz
virusbreakend/virusbreakenddb_20210401.tar.gz
'

parallel -j4 wget -c -x -nH -P reference_data/ ${base_url}/{} ::: ${fps}
find reference_data/ -name '*.tar.gz' | parallel -j0 'cd {//} && tar -xzvf {/}'
```

Create Nextflow config file for local reference data

```bash
cat <<EOF > refdata.local.config
params {
    genomes {
        'GRCh37_hmf' {
            fasta           = "$(pwd)/genomes/GRCh37_hmf/Homo_sapiens.GRCh37.GATK.illumina.fasta"
            fai             = "$(pwd)/genomes/GRCh37_hmf/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai"
            dict            = "$(pwd)/genomes/GRCh37_hmf/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict"
            bwa_index       = "$(pwd)/genomes/GRCh37_hmf/bwa_index/0.7.17-r1188.tar.gz"
            bwa_index_bseq  = "$(pwd)/genomes/GRCh37_hmf/bwa_index/2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta.0123"
            bwa_index_biidx = "$(pwd)/genomes/GRCh37_hmf/bwa_index/2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta.bwt.2bit.64"
            bwa_index_image = "$(pwd)/genomes/GRCh37_hmf/bwa_index_image/0.7.17-r1188/Homo_sapiens.GRCh37.GATK.illumina.fasta.img"
            gridss_index    = "$(pwd)/genomes/GRCh37_hmf/gridss_index/2.13.2/Homo_sapiens.GRCh37.GATK.illumina.fasta.gridsscache"
            star_index      = "$(pwd)/genomes/GRCh37_hmf/star_index/gencode_19/2.7.3a.tar.gz"
        }
        'GRCh38_hmf' {
            fasta           = "$(pwd)/reference_data/genomes/GRCh38_hmf/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
            fai             = "$(pwd)/reference_data/genomes/GRCh38_hmf/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
            dict            = "$(pwd)/reference_data/genomes/GRCh38_hmf/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict"
            bwa_index       = "$(pwd)/reference_data/genomes/GRCh38_hmf/bwa_index/0.7.17-r1188/"
            bwa_index_bseq  = "$(pwd)/reference_data/genomes/GRCh38_hmf/bwa_index/2.2.1/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.0123"
            bwa_index_biidx = "$(pwd)/reference_data/genomes/GRCh38_hmf/bwa_index/2.2.1/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwt.2bit.64"
            bwa_index_image = "$(pwd)/reference_data/genomes/GRCh38_hmf/bwa_index_image/0.7.17-r1188/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img"
            gridss_index    = "$(pwd)/reference_data/genomes/GRCh38_hmf/gridss_index/2.13.2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gridsscache"
            star_index      = "$(pwd)/reference_data/genomes/GRCh38_hmf/star_index/gencode_38/2.7.3a/"
        }
    }

    ref_data_hmf_data_path = "$(pwd)/reference_data/hmf_reference_data/hmftools/5.34_38--2/"
    ref_data_panel_data_path = "$(pwd)/reference_data/hmf_reference_data/panels/tso500_5.34_38--1/"
    ref_data_virusbreakenddb_path = "$(pwd)/reference_data/virusbreakend/virusbreakenddb_20210401/"
}
EOF
```

Run oncoanalyser with local reference data

> Assumes existing samplesheet at `samplesheet.csv`

```bash
nextflow run oncoanalyser/main.nf \
  \
  -config refdata.local.config \
  -profile docker \
  -revision v0.3.1 \
  \
  --mode targeted \
  --panel tso500 \
  --genome GRCh38_hmf \
  \
  --input samplesheet.csv \
  --outdir output/
```
