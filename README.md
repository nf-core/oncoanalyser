&nbsp;
&nbsp;
&nbsp;
<p align="center">
🚧🚨 <em>Under development</em> 🚨🚧
</p>

# HMF tools

A comprehensive cancer WGS/WTS analysis and reporting pipeline from the Hartwig Medical Foundation.

## Table of contents

* [Quick start](#quick-start)
* [Pipeline modes](#pipeline-modes)
* [Pipeline tests](#pipeline-tests)
* [License](#license)

## Quick start

1. Install dependencies.
> Assumes `Conda` is installed.
```bash
conda create -y \
  -n hmftools \
  -c bioconda \
  -c conda-forge \
  -c defaults \
  'bcftools >=1.15' \
  'dvc >=2.12.0' \
  'git >=2.37.0' \
  'nextflow >=22.04.0'
conda activate hmftools
```

2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Download the pipeline and reference data.
```bash
# Download pipeline and reference data
git clone https://github.com/umccr/hmftools && cd ./hmftools/hmftools/

git clone https://github.com/umccr/reference_data -b dev reference_data_gitrepo/ && cd reference_data_gitrepo/
dvc pull reference_data/{genomes,hmftools}/ -r storage-s3 && cd ../
ln -s reference_data_gitrepo/reference_data/ reference_data
```

4. Run the pipeline with your data.
```bash
nextflow run ./main.nf -profile docker --input <SAMPLESHEET_FP> --outdir <OUTPUT_DIRECTORY>
```

## Pipeline modes

Several modes of execution are offered and can be accessed using the `--mode` argument.

| Name                  | Description                                   |
| ---                   | ---                                           |
| `full`                | Run all processes in the pipeline             |
| `gridss-purple-linx`  | Run AMBER, COBALT, GRIDSS, GRIPSS, and PURPLE |
| `gridss`              | Run only GRIDSS and GRIPSS                    |
| `purple`              | Run only PURPLE                               |
| `linx`                | Run only LINX                                 |
| `lilac`               | Run only LILAC                                |
| `teal`                | Run only TEAL                                 |

See [examples/README.md](examples/README.md) for example samplesheet required for each mode.

## Pipeline tests

The internal pipeline logic can be tested using a stub run. No actual outputs are generated in a stub run and
only the `stub` definition block of a process is executed. Hence, it completes in a short amount of time but does not
test validity of actual processes beyond declared inputs and outputs.

```bash
nextflow run main.nf --input ./assets/samplesheet.tsv --outdir output_stub/ --max_memory '1.GB' -stub
```

A more comprehensive test that involves both the process `script` definition block and internal pipeline logic can be
run with the supplied downsampled test dataset. This of course takes longer (on the order of minutes) but replicates a
full run and generates 'real' process outputs.

```bash
nextflow run ./main.nf -profile docker,test --outdir output_test/ --max_memory '6.GB'
```

## License

Software and code in this repository are under [GNU General Public License
v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) unless otherwise indicated.
