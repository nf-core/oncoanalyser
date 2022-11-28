# Examples

## Full

| Filetype  | Keyword   | Description                                                               | Type     |
| ---       | ---       | ---                                                                       | ---      |
| BAM (WGS) | `bam_wgs` | WGS read alignments                                                       | Required |
| BAM (WTS) | `bam_wts` | WTS read alignments                                                       | Optional |
| SV VCF    | `vcf`     | SV VCF produced by an external caller [_used to filter reads for GRIDSS_] | Optional |

```text
id        subject_name   sample_name          sample_type  filetype  filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam_wgs   /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam_wts   /path/to/tumor_bam/sample_one_tumor.rna.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        vcf_sv    /path/to/tumor_sv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam_wgs   /path/to/normal_bam/sample_one_normal.bam
STWO-1    SUBJECT_TWO    SAMPLE_TWO_TUMOR-1   tumor        bam_wgs   /path/to/tumor_bam/sample_two_tumor_one.bam
STWO-1    SUBJECT_TWO    SAMPLE_TWO_NORMAL    normal       bam_wgs   /path/to/normal_bam/sample_two_normal.bam
STWO-2    SUBJECT_TWO    SAMPLE_TWO_TUMOR-2   tumor        bam_wgs   /path/to/tumor_bam/sample_two_tumor_two.bam
STWO-2    SUBJECT_TWO    SAMPLE_TWO_NORMAL    normal       bam_wgs   /path/to/normal_bam/sample_two_normal.bam
STRHEE-1  SUBJECT_THREE  SAMPLE_THREE_TUMOR   tumor        bam_wgs   /path/to/tumor_bam/sample_three_tumor.bam
STRHEE-1  SUBJECT_THREE  SAMPLE_THREE_TUMOR   tumor        vcf_sv    /path/to/tumor_sv_vcf/sample_three_tumor.vcf.gz
STRHEE-1  SUBJECT_THREE  SAMPLE_THREE_NORMAL  normal       bam_wgs   /path/to/normal_bam/sample_three_normal.bam
```

## GRIDSS

See [Full section](#full)

## PURPLE

| Filetype                      | Keyword              | Description                 | Type     |
| ---                           | ---                  | ---                         | ---      |
| AMBER directory               | `amber_dir`          | AMBER output directory      | Required |
| COBALT directory              | `cobalt_dir`         | COBALT output directory     | Required |
| GRIPSS SV VCF (hard filtered) | `vcf_sv_gripss_hard` | Hard filtered GRIPSS SV VCF | Required |
| GRIPSS SV VCF (soft filtered) | `vcf_sv_gripss_soft` | Soft filtered GRIPSS SV VCF | Required |
| SNV/MNV and INDEL VCF         | `vcf_smlv`           | Small SNV/MNV VCF           | Optional |

```text
id        subject_name   sample_name          sample_type  filetype             filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        amber_dir            /path/to/amber_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        cobalt_dir           /path/to/cobalt_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        vcf_sv_gripss_hard   /path/to/tumor_gripss_hard_sv/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        vcf_sv_gripss_soft   /path/to/tumor_gripss_soft_sv/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        vcf_smlv             /path/to/tumor_smlv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       vcf_smlv             /path/to/normal_smlv_vcf/sample_one_normal.vcf.gz
```

## LINX

| Filetype                      | Keyword              | Description                                   | Type     |
| ---                           | ---                  | ---                                           | ---      |
| PURPLE directory              | `purple_dir`         | PURPLE output directory [_LINX somatic_]      | Required |
| GRIPSS SV VCF (hard filtered) | `vcf_sv_gripss_hard` | Hard filtered GRIPSS SV VCF [_LINX germline_] | Required |

```text
id        subject_name   sample_name          sample_type  filetype             filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        purple_dir           /path/to/purple_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       vcf_sv_gripss_hard   /path/to/normal_gripss_hard_sv/sample_one_normal.vcf.gz
```

## GRIDSS-PURPLE-LINX

| Filetype              | Keyword    | Description                                                               | Type     |
| ---                   | ---        | ---                                                                       | ---      |
| BAM (WGS)             | `bam`      | WGS read alignments                                                       | Required |
| SV VCF                | `vcf_sv`   | SV VCF produced by an external caller [_used to filter reads for GRIDSS_] | Optional |
| SNV/MNV and INDEL VCF | `vcf_smlv` | Small SNV/MNV VCF                                                         | Optional |

```text
id        subject_name   sample_name          sample_type  filetype  filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam_wgs   /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        vcf_sv    /path/to/tumor_sv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam_wgs   /path/to/normal_bam/sample_one_normal.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        vcf_smlv  /path/to/tumor_smlv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       vcf_smlv  /path/to/normal_smlv_vcf/sample_one_normal.vcf.gz
```


## LILAC

| Filetype              | Keyword      | Description             | Type     |
| ---                   | ---          | ---                     | ---      |
| BAM (WGS)             | `bam_wgs`    | WGS read alignments     | Required |
| PURPLE directory      | `purple_dir` | PURPLE output directory | Required |

```text
id        subject_name   sample_name          sample_type  filetype     filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam_wgs      /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        purple_dir   /path/to/purple_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam_wgs      /path/to/normal_bam/sample_one_normal.bam
```

## TEAL

| Filetype              | Keyword      | Description             | Type     |
| ---                   | ---          | ---                     | ---      |
| BAM (WGS)             | `bam_wgs`    | WGS read alignments     | Required |
| COBALT directory      | `cobalt_dir` | COBALT output directory | Required |
| PURPLE directory      | `purple_dir` | PURPLE output directory | Required |

```text
id        subject_name   sample_name          sample_type  filetype     filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam_wgs      /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        purple_dir   /path/to/purple_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        cobalt_dir   /path/to/cobalt_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam_wgs      /path/to/normal_bam/sample_one_normal.bam
```
