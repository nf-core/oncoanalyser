#!/usr/bin/env python3
"""Generate the software versions tables in README.md from module environment.yml files.

Usage:
    python3 .github/scripts/update_software_versions.py          # rewrite README.md in place
    python3 .github/scripts/update_software_versions.py --check  # exit 1 if README.md is out of date
"""

import argparse
import pathlib
import re
import sys

REPO_DIR = pathlib.Path(__file__).resolve().parents[2]
README_FP = REPO_DIR / 'README.md'
MODULES_DIR = REPO_DIR / 'modules'

MARKER_START = '<!-- software-versions:start -->'
MARKER_END = '<!-- software-versions:end -->'

# Display name, module environment.yml (relative to modules/), conda package name
WIGITS_TOOLS = [
    ('AMBER', 'local/amber', 'hmftools-amber'),
    ('BamTools', 'local/bamtools', 'hmftools-bam-tools'),
    ('CHORD', 'local/chord', 'hmftools-chord'),
    ('CIDER', 'local/cider', 'hmftools-cider'),
    ('COBALT', 'local/cobalt/run', 'hmftools-cobalt'),
    ('CUPPA', 'local/cuppa', 'hmftools-cuppa'),
    ('ESVEE', 'local/esvee', 'hmftools-esvee'),
    ('FastqTools', 'local/fastqtools', 'hmftools-fastq-tools'),
    ('ISOFOX', 'local/isofox/run', 'hmftools-isofox'),
    ('LILAC', 'local/lilac', 'hmftools-lilac'),
    ('LINX', 'local/linx/somatic', 'hmftools-linx'),
    ('NEO', 'local/neo/finder', 'hmftools-neo'),
    ('ORANGE', 'local/orange', 'hmftools-orange'),
    ('PAVE', 'local/pave/somatic', 'hmftools-pave'),
    ('PEACH', 'local/peach', 'hmftools-peach'),
    ('PURPLE', 'local/purple', 'hmftools-purple'),
    ('QSEE', 'local/qsee', 'hmftools-qsee'),
    ('REDUX', 'local/redux', 'hmftools-redux'),
    ('SAGE', 'local/sage/somatic', 'hmftools-sage'),
    ('SIGS', 'local/sigs', 'hmftools-sigs'),
    ('TEAL', 'local/teal/pipeline', 'hmftools-teal'),
    ('VirusInterpreter', 'local/virusinterpreter', 'hmftools-virus-interpreter'),
    ('WISP', 'local/wisp', 'hmftools-wisp'),
]

# Versions are taken from the modules used during analysis rather than reference data creation (e.g. BWA-MEM2 is read
# from the alignment module, not the nf-core index module)
OTHER_TOOLS = [
    ('BWA-MEM2', 'local/bwa-mem2/mem', 'bwa-mem2'),
    ('fastp', 'local/fastp/split', 'fastp'),
    ('GATK4 (Picard MarkDuplicates)', 'nf-core/gatk4/markduplicates', 'gatk4'),
    ('GRIDSS (VIRUSbreakend)', 'local/virusbreakend', 'gridss'),
    ('linxreport', 'local/linxreport', 'r-linxreport'),
    ('MultiQC', 'local/multiqc', 'multiqc'),
    ('RepeatMasker', 'local/virusbreakend', 'repeatmasker'),
    ('sambamba', 'local/bwa-mem2/mem', 'sambamba'),
    ('samtools', 'local/bwa-mem2/mem', 'samtools'),
    ('STAR', 'local/star/align', 'star'),
]


def parse_dependencies(env_fp):
    # Minimal parser for the flat dependency lists used in module environment.yml files, avoids requiring PyYAML
    dependencies = dict()
    in_dependencies = False
    for line in env_fp.read_text().splitlines():
        if line.startswith('dependencies:'):
            in_dependencies = True
            continue
        if not in_dependencies:
            continue
        if not line.startswith(' '):
            break
        if not (m := re.match(r'^\s+-\s+(?:[\w.-]+::)?([\w.-]+)\s*(?:[=<>]=?)?([^=\s]*)', line)):
            continue
        name, version = m.groups()
        dependencies[name] = version
    return dependencies


def get_version(module, package):
    env_fp = MODULES_DIR / module / 'environment.yml'
    if not env_fp.exists():
        sys.exit(f'error: {env_fp.relative_to(REPO_DIR)} does not exist')
    if not (version := parse_dependencies(env_fp).get(package)):
        sys.exit(f'error: no pinned version for {package} in {env_fp.relative_to(REPO_DIR)}')
    return version


def check_wigits_coverage():
    # Fail if a WiGiTS tool is added to the pipeline without being added to the table above
    listed = {package for _, _, package in WIGITS_TOOLS}
    for env_fp in sorted(MODULES_DIR.glob('**/environment.yml')):
        for package in parse_dependencies(env_fp):
            if package.startswith('hmftools-') and package not in listed:
                sys.exit(f'error: {package} in {env_fp.relative_to(REPO_DIR)} is not listed in {pathlib.Path(__file__).name}')


def format_name(name):
    # Code format the tool name only, leaving any parenthesised note as plain text
    tool, _, note = name.partition(' (')
    return f'`{tool}` ({note}' if note else f'`{tool}`'


def render_table(header, tools):
    rows = [(format_name(name), get_version(module, package)) for name, module, package in tools]
    name_width = max(len(header), *(len(name) for name, _ in rows))
    version_width = max(len('Version'), *(len(version) for _, version in rows))
    lines = [
        f'| {header:<{name_width}} | {"Version":<{version_width}} |',
        f'| {"-" * name_width} | {"-" * version_width} |',
    ]
    lines.extend(f'| {name:<{name_width}} | {version:<{version_width}} |' for name, version in rows)
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--check', action='store_true', help='exit 1 if README.md is out of date')
    args = parser.parse_args()

    check_wigits_coverage()

    tables = '\n\n'.join([
        render_table('WiGiTS tool', WIGITS_TOOLS),
        render_table('Other tool', OTHER_TOOLS),
    ])

    readme = README_FP.read_text()
    pattern = re.compile(f'{re.escape(MARKER_START)}.*?{re.escape(MARKER_END)}', re.DOTALL)
    if not pattern.search(readme):
        sys.exit(f'error: software versions markers not found in {README_FP.name}')
    readme_new = pattern.sub(lambda _: f'{MARKER_START}\n\n{tables}\n\n{MARKER_END}', readme)

    if args.check:
        if readme_new != readme:
            sys.exit(f'error: software versions in {README_FP.name} are out of date, run {pathlib.Path(__file__).name}')
        return

    README_FP.write_text(readme_new)


if __name__ == '__main__':
    main()
