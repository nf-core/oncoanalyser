//
// Tool output filename accessors for the nf-core/oncoanalyser pipeline
//

def getVcf(sample_id, dir, suffix) {
    if (! dir) { return null }
    return dir.resolve("${sample_id}.${suffix}.vcf.gz")
}

def getVcfTbi(vcf) {
    return vcf ? nextflow.Nextflow.file("${vcf}.tbi") : null
}

def getPurpleSomaticVcf(sample_id, purple_dir) {
    return getVcf(sample_id, purple_dir, 'purple.somatic')
}

def getPurpleGermlineVcf(sample_id, purple_dir) {
    return getVcf(sample_id, purple_dir, 'purple.germline')
}

def getPurpleSvVcf(sample_id, purple_dir) {
    return getVcf(sample_id, purple_dir, 'purple.sv')
}

def getPurpleSvGermlineVcf(sample_id, purple_dir) {
    return getVcf(sample_id, purple_dir, 'purple.sv.germline')
}

def getSageSomaticVcf(sample_id, sage_dir) {
    return getVcf(sample_id, sage_dir, 'sage.somatic')
}

def getSageGermlineVcf(sample_id, sage_dir) {
    return getVcf(sample_id, sage_dir, 'sage.germline')
}

def getSageAppendVcf(sample_id, sage_append_dir) {
    return getVcf(sample_id, sage_append_dir, 'sage.append')
}

def getCircosPlot(sample_id, purple_dir) {
    if (! purple_dir) { return null }
    return purple_dir.resolve("plot/${sample_id}.circos.png")
}

def getCuppaVisData(sample_id, cuppa_dir) {
    if (! cuppa_dir) { return null }
    return cuppa_dir.resolve("${sample_id}.cuppa.vis_data.tsv")
}
