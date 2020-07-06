import re

names = ['Vaccinia virus (strain Western Reserve) (VACV) (Vaccinia virus (strain WR))',
         'Horsepox virus (HSPV)', 'Rabbitpox virus', 'Vaccinia virus', 'Cowpox virus (CPV)']

rgx = re.compile(r'(\w+\svirus)')
for name in names:
    m = re.match(rgx, name)
    print(m.group(0))