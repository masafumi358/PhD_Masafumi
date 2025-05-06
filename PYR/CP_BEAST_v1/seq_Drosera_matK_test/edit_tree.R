#set working directory
wd = "~/Github/PhD_Masafumi/PYR/CP_BEAST_v1/seq_Drosera_matK_test/"
setwd(wd)

# Load ape package
library(ape)

# Read your BEAST MCC tree
trees <- read.nexus("Drosera_beast_matK_msl.trees")
tr <- trees[[1]]  # Pick one sample
tr$edge.length <- tr$edge.length*10

# Set up PDF output
pdf("Drosera_beast_tree_simple.pdf", width = 7.5, height = 9)

# Plot the tree with clearer spacing
plot(tr,
     cex = 0.7,              # Tip label size
     label.offset = 0.01,    # Space between tip and label
     no.margin = FALSE,      # Adds margin
     edge.width = 1.2        # Thicker branches
)

# Add scale bar
axisPhylo()

# Close PDF device
dev.off()
system("open Drosera_beast_tree_simple.pdf")

