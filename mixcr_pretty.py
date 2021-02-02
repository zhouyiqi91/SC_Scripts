
import sys
import os

def mixcr(outdir, sample, input_file, thread, species):
    report = f"{outdir}/{sample}_align.txt"
    not_align_fq = f"{outdir}/not_align.fq"
    read2_vdjca = f"{outdir}/read2.vdjca"
    alignments = f"{outdir}/{sample}_alignments.txt"

    cmd = f"""
mixcr align \
--force-overwrite \
--species {species} \
-t {thread} \
--not-aligned-R1 {not_align_fq} \
--report {report} \
-OallowPartialAlignments=true \
-OvParameters.geneFeatureToAlign=VTranscriptWithP \
{input_file} \
{read2_vdjca}
mixcr exportAlignments \
{read2_vdjca} {alignments} \
-readIds --force-overwrite -vGene -dGene -jGene -cGene \
-nFeature CDR3 -aaFeature CDR3\n
mixcr exportAlignmentsPretty {read2_vdjca} > {sample}_pretty.txt
"""

    print(cmd)
    os.system(cmd)
    return alignments


input_file = sys.argv[1]
outdir = './'
sample = os.path.basename(input_file)
thread = 2
species = 'hs'
mixcr(outdir, sample, input_file, thread, species)