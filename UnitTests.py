from NequickG import SlantRayAnalysis
import numpy as np

h1, lat1, lon1, h2, lat2, lon2 = (150, 40, 0, 20000, 55, 0)
sra = SlantRayAnalysis(150, 40, 0, 20000, 55, 0)
print sra.rp, sra.latp, sra.lonp
s1, s2 = sra.ray_endpoints()
assert(s1 <s2 )
hs, lats, lons = sra.ray_coords(np.linspace(7000, 15000, 10))

print lats
print lons

assert not np.any(lats > lat2)
assert not np.any(lats < lat1)
assert not np.any(lons > lon2)
assert not np.any(lons < lon1)
assert not np.any(hs > h2)
assert not np.any(hs < h1)


assert (s2 ** 2 + sra.rp ** 2 - (20000 + 6371.2) ** 2) < 0.1
assert (s1 ** 2 + sra.rp ** 2 - (150 + 6371.2) ** 2) < 0.1


