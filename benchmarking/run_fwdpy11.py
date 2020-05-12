import sys

import fwdpy11

N = int(sys.argv[1])
seed = int(sys.argv[2])

pdict = {
    "nregions": [],
    "sregions": [],
    "recregions": [],
    "rates": [0, 0, 0],
    "gvalue": fwdpy11.Multiplicative(2.0),
    "simlen": 5 * N,
}
mp = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation(N, 1.0)
rng = fwdpy11.GSLrng(seed)
fwdpy11.evolvets(rng, pop, mp, 100, suppress_table_indexing=True)
ts = pop.dump_tables_to_tskit()
ts.dump("fwdpy11.trees")
