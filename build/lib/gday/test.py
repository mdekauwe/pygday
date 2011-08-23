import pstats
p = pstats.Stats('x')
p.sort_stats('name')
p.print_stats()