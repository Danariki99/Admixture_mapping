import pickle
import numpy as np
import pandas as pd
import os
from snputils.genobj.snp.snpobj import SNPObject
from snputils.io.snp.read.vcf import VCFReader
from snputils.io.snp.write.vcf import VCFWriter

print('reading file')

reader = VCFReader("/private/groups/ioannidislab/smeriglio/trial/merged_vcfs/ancestry_EUR.vcf")

snpobj = reader.read()

print('finished reading file')

snpobj.calldata_gt = np.logical_not(snpobj.calldata_gt).astype(int)
new_filename = '/private/groups/ioannidislab/smeriglio/trial/merged_vcfs/negative_ancestry_EUR_merged.vcf'

vcf_writer = VCFWriter(snpobj, new_filename)
vcf_writer.write()

print('finished writing file')