#!/usr/bin/env python3
import ssam, inspect

print("=== SSAMAnalysis.find_localmax ===")
print(inspect.signature(ssam.SSAMAnalysis.find_localmax))
print(inspect.getdoc(ssam.SSAMAnalysis.find_localmax))

print("\n=== SSAMAnalysis.cluster_vectors ===")
print(inspect.signature(ssam.SSAMAnalysis.cluster_vectors))

print("\n=== SSAMAnalysis.map_celltypes ===")
print(inspect.signature(ssam.SSAMAnalysis.map_celltypes))

print("\n=== SSAMDataset attributes ===")
print([m for m in dir(ssam.SSAMDataset) if not m.startswith('_')])
