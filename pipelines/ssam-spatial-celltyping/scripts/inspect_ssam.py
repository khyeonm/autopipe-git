import ssam, inspect, os

d = os.path.dirname(ssam.__file__)
out = []

out.append('=== SSAM package location ===')
out.append(d)
out.append('')

out.append('=== Files in ssam package ===')
out += os.listdir(d)
out.append('')

out.append('=== ssam exports (dir) ===')
out.append(str(dir(ssam)))
out.append('')

for name in ['SSAMDataset', 'SSAMAnalysis']:
    cls = getattr(ssam, name, None)
    if cls:
        out.append('=== ' + name + '.__init__ ===')
        try:
            out.append(inspect.getsource(cls.__init__))
        except Exception as e:
            out.append('(source unavailable: ' + str(e) + ')')
        out.append('')

ana = ssam.SSAMAnalysis
for method in ['run_kde', 'normalize_vectors', 'threshold_vectors', 'map_celltypes', 'find_localmax', 'run_umap']:
    fn = getattr(ana, method, None)
    if fn:
        out.append('=== SSAMAnalysis.' + method + ' ===')
        try:
            out.append(inspect.getsource(fn))
        except Exception as e:
            out.append('(source unavailable: ' + str(e) + ')')
        out.append('')

with open('/output/ssam_source.txt', 'w') as f:
    f.write('\n'.join(out))
print('Done')

