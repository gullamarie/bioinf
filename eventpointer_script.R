# Setup

.libPaths("/usr/local/lib/R/site-library")
.libPaths()
print("Done setting libpath.")

print("Loading Eventpointer")
library("EventPointer")
print("Done loading Evenpointer")

# Number of cores to use
NUM_CORES <- 10
RESULTS_DIR <- "/directory/eventpointer/samples/results"

# Run EP
# Read list of BAM filenames from a file
samples <- scan("/directory/eventpointer/samples/results/index_bams/all_filenames.tsv", what="list", sep="\n")
path_to_samples <- "/directory/bam_files"
path_to_gtf <- "/directory/reference_annotation/gencode.v19.annotation.gtf"

print("PrepareBam START")
sg_feature_counts <- PrepareBam_EP(Samples=samples, SamplePath=path_to_samples, Ref_Transc="GTF", fileTransc=path_to_gtf, cores=NUM_CORES)
# Save preliminary results
save.image(file.path(RESULTS_DIR, "preparebam_finished.RData"))
print("PrepareBam STOP")

# Run EventDetection function
print("EventDetection START")
data(sg_feature_counts)
# Directory to write the EventsFound_RNASeq.txt file
txt_path <- RESULTS_DIR
AllEvents_RNASeq <- EventDetection(sg_feature_counts, cores=NUM_CORES, Path=txt_path)
# Save preliminary results
save.image(file.path(RESULTS_DIR, "eventdetection_finished.RData"))
print("EventDetection STOP")


