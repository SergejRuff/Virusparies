# Virusparies

VirusHunterGatherer is a data-driven tool designed for high-throughput virus discovery.
The process involves Virushunter conducting sensitive homology-based detection of
viral sequence reads in unprocessed data. This approach identifies the most conserved regions
of a virus, which serve as seeds in the Virusgatherer step for a seed-based assembly
of full-length viral genome sequences.

The Virusparies package provides a set of plotting functions tailored for visualizing
VirushunterGatherer Hittables. The name draws inspiration from the hunter-gatherer metaphor,
with "paries" derived from Latin meaning "wall". It symbolizes the parietal art left by
ancient hunters and gatherers on walls, summarizing their stories and beliefs.

# Installation

## First install the "remotes" package

install.packages("remotes")

## Then install the Virusparies package

remotes::install_github("SergejRuff/Virusparies")


# Contributions

Sergej Ruff formulated the idea behind Virusparies and was responsible for its implementation.

Chris Lauber and Li Chuin Chong from Twincore - Centre for Experimental and Clinical Infection Research provided ideas for improvements. 
Chris Lauber is the main developer behind the VirusHuntergatherer software and the Group Leader of the Computational Virology working group at the Institute for Experimental Virology, TWINCORE.
