#!/usr/bin/env python3
import argparse
import collections
import enum
import pathlib
import sys


HEADER_COLS = (
    'id',
    'subject_name',
    'sample_name',
    'sample_type',
    'filetype',
    'filepath',
)


class FileType(enum.Enum):

    BAM_WGS = 'bam_wgs'
    BAM_WTS = 'bam_wts'
    VCF_SV = 'vcf_sv'
    SMLV_VCF = 'vcf_smlv'
    VCF_SV_GRIPSS_SOFT = 'vcf_sv_gripss_soft'
    VCF_SV_GRIPSS_HARD = 'vcf_sv_gripss_hard'
    AMBER_DIR = 'amber_dir'
    COBALT_DIR = 'cobalt_dir'
    PURPLE_DIR = 'purple_dir'

    def __repr__(self):
        return self.value


class SampleType(enum.Enum):

    TUMOR = 'tumor'
    NORMAL = 'normal'

    def __repr__(self):
        return self.value


FILETYPES_EXPECTED = {
    'full': {
        'required': [
            (SampleType.TUMOR, FileType.BAM_WGS),
            (SampleType.NORMAL, FileType.BAM_WGS),
        ],
        'optional': [
            (SampleType.TUMOR, FileType.BAM_WTS),
            (SampleType.TUMOR, FileType.VCF_SV),
            (SampleType.NORMAL, FileType.VCF_SV),
        ],
    },
    'gridss_purple_linx': {
        'required': [
            (SampleType.TUMOR, FileType.BAM_WGS),
            (SampleType.NORMAL, FileType.BAM_WGS),
        ],
        'optional': [
            (SampleType.TUMOR, FileType.VCF_SV),
            (SampleType.NORMAL, FileType.VCF_SV),
            (SampleType.TUMOR, FileType.SMLV_VCF),
            (SampleType.NORMAL, FileType.SMLV_VCF)
        ],
    },
    'gridss': {
        'required': [
            (SampleType.TUMOR, FileType.BAM_WGS),
            (SampleType.NORMAL, FileType.BAM_WGS),
        ],
        'optional': [
            (SampleType.TUMOR, FileType.VCF_SV),
            (SampleType.NORMAL, FileType.VCF_SV),
        ],
    },
    'purple': {
        'required': [
            (SampleType.TUMOR, FileType.AMBER_DIR),
            (SampleType.TUMOR, FileType.COBALT_DIR),
            (SampleType.TUMOR, FileType.VCF_SV_GRIPSS_HARD),
            (SampleType.TUMOR, FileType.VCF_SV_GRIPSS_SOFT),
        ],
        'optional': [
            (SampleType.TUMOR, FileType.SMLV_VCF),
            (SampleType.NORMAL, FileType.SMLV_VCF),
        ],
    },
    # NOTE(SW): expectation that we always at least want to run LINX somatic
    'linx': {
        'required': [
            (SampleType.TUMOR, FileType.PURPLE_DIR),
        ],
        'optional': [
            (SampleType.NORMAL, FileType.VCF_SV_GRIPSS_HARD),
        ],
    },
    'lilac': {
        'required': [
            (SampleType.TUMOR, FileType.BAM_WGS),
            (SampleType.NORMAL, FileType.BAM_WGS),
            (SampleType.TUMOR, FileType.PURPLE_DIR),
        ],
        'optional': [],
    },
    'teal': {
        'required': [
            (SampleType.TUMOR, FileType.BAM_WGS),
            (SampleType.NORMAL, FileType.BAM_WGS),
            (SampleType.TUMOR, FileType.COBALT_DIR),
            (SampleType.TUMOR, FileType.PURPLE_DIR),
        ],
        'optional': [],
    },
}


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_fp',
        required=True,
        type=pathlib.Path,
        help='Input samplesheet filepath',
    )
    parser.add_argument(
        '--mode',
        required=True,
        type=str,
        choices=FILETYPES_EXPECTED.keys(),
        help='Pipeline execution mode',
    )
    args = parser.parse_args()
    if not args.input_fp.exists():
        parser.error(f'Input file {args.input_fp} does not exist')
    return args


def main():
    # Get commandline arguments
    args = get_arguments()

    # Read data
    runs = collections.defaultdict(list)
    with args.input_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        check_header(header_tokens, HEADER_COLS)
        for lts in line_token_gen:
            record = {k: v for k, v in zip(header_tokens, lts)}
            record['sample_type_enum'] = get_enum(record, 'sample_type', SampleType)
            record['filetype_enum'] = get_enum(record, 'filetype', FileType)
            runs[record['id']].append(record)

    # Validate data
    for rid, records in runs.items():
        input_types = {(r['sample_type_enum'], r['filetype_enum']) for r in records}
        check_input_types(
            input_types,
            FILETYPES_EXPECTED[args.mode]['required'],
            FILETYPES_EXPECTED[args.mode]['optional'],
        )
        check_subject_names(rid, records)
        check_bam_sample_names(records, FileType.BAM_WGS)
        check_bam_sample_names(records, FileType.BAM_WTS)
        check_sample_type_sample_names(rid, records)


def get_enum(record, name, enum_class):
    value = record[name]
    try:
        return enum_class(value)
    except ValueError:
        print(f'Got invalid {name} for {record["id"]}: {value}', file=sys.stderr)
        sys.exit(1)


def check_subject_names(rid, records):
    subject_names = {r['subject_name'] for r in records}
    if len(subject_names) > 1:
        subject_names_str = ', '.join(subject_names)
        msg = f'Got multiple subject names for {rid}: {subject_names_str}'
        print(msg, file=sys.stderr)
        sys.exit(1)


def check_bam_sample_names(records, bam_filetype):
    record_tumor_bam = None
    record_normal_bam = None
    for record in records:
        if record['filetype_enum'] != bam_filetype:
            continue
        if record['sample_type_enum'] == SampleType.TUMOR:
            record_tumor_bam = record
        elif record['sample_type_enum'] == SampleType.NORMAL:
            record_normal_bam = record

    if record_tumor_bam is None or record_normal_bam is None:
        return

    if record_tumor_bam['sample_name'] == record_normal_bam['sample_name']:
        sample_name = record_tumor_bam['sample_name']
        print(f'Got identical sample names for \'{repr(bam_filetype)}\' BAMs: {sample_name}', file=sys.stderr)
        sys.exit(1)


def check_sample_type_sample_names(rid, records):
    sample_names_tumor = set()
    sample_names_normal = set()
    for record in records:
        if record['sample_type_enum'] == SampleType.TUMOR:
            sample_names_tumor.add(record['sample_name'])
        elif record['sample_type_enum'] == SampleType.NORMAL:
            sample_names_normal.add(record['sample_name'])
        else:
            assert False
    if len(sample_names_tumor) > 1:
        sn_tumor_str = ', '.join(sample_names_tumor)
        print(f'Got mismatch sample names for {rid} tumor: {sn_tumor_str}', file=sys.stderr)
        sys.exit(1)
    if len(sample_names_normal) > 1:
        sn_normal_str = ', '.join(sample_names_normal)
        print(f'Got mismatch sample names for {rid} normal: {sn_normal_str}', file=sys.stderr)
        sys.exit(1)


def check_header(header_tokens, required_cols):
    check_set_difference(header_tokens, required_cols, 'Found unknown header column')
    check_set_difference(required_cols, header_tokens, 'Missing required column')


def check_input_types(input_types, required, optional):
    check_set_difference(input_types, required + optional, 'Found unknown input type')
    check_set_difference(required, input_types, 'Missing required input type')


def check_set_difference(a, b, message_base):
    d = set(a).difference(b)
    if d :
        d_str = ', '.join(str(e) for e in d)
        plurality = 's' if len(d) > 1 else ''
        print(f'{message_base}{plurality}: {d_str}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
