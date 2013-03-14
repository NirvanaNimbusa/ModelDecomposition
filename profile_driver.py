#!/usr/bin/env python

import cProfile,pstats,detile_image

cProfile.run("""detile_image.main('ldconfig.csv', 'xelnagatiles.jpg', '78', '100','2','2')""",'detiling_profile')

print "-----------------------------------------------------------"

p = pstats.Stats('detiling_profile')
p.strip_dirs().sort_stats('cum','time').print_stats()
p.print_callees()
