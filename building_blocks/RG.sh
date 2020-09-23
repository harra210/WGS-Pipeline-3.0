#!/bin/bash
#
# Read group fields to be inserted for GATK at either BWA or Picard
# ID = Read group identifier, which will be used as the sample name
# LB = Library group identifier, which is the library that each sample was generated on
# PL = Platform identifier, which will normally be Illumina unless otherwise specified
# SM = Sample, which will be the samples name
# PU = Platform unit, which is the flowcell barcode
