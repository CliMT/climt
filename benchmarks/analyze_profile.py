import pstats
import sys

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'profile.stats'

p = pstats.Stats(filename)
print("Top 30 functions by cumulative time:")
p.strip_dirs().sort_stats('cumulative').print_stats(30)

print("\nTop 30 functions by total time (internal time):")
p.strip_dirs().sort_stats('tottime').print_stats(30)
