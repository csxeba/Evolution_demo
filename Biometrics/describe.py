def heterozygosity(pop):
    "Averaged heterozygosity on all locuses"
    holder = [heterozygote(ind) for ind in pop.individuals]
    return round(sum(holder) / len(holder), 2)


def heterozygote(ind):
    holder = []
    for locusA, locusB in zip(ind.chromosomeA, ind.chromosomeB):
        holder.append(int(locusA != locusB))
    return sum(holder) / len(holder)


def mutants(pop):
    holder = [int(ind.mutant) for ind in pop.individuals]
    return round(sum(holder) / len(holder), 2)


def gene_loss(pop):
    """Equals the number of different genes divided by
    the initial number of individual genes"""
    lost = round((len(gene_pool(pop)) / len(pop.init_gene_pool)) * 100, 2)
    lost = 100 - lost  # Possibly not?
    return lost


def gene_pool(pop):
    holder = [tuple(x.chromosomeA) for x in pop.individuals] + \
             [tuple(x.chromosomeB) for x in pop.individuals]
    return set(holder)


def fitness(pop):
    return pop.fitness()


def topn(pop, n, best=True):
    distances = [avgDistance(ind, pop) for ind in pop.individuals]
    distance_tr = sorted(distances, reverse=not bool(best))[0]
    if best:
        inds = [ind for ind in pop.individuals
                if avgDistance(ind, pop) <= distance_tr]
    else:
        inds = [ind for ind in pop.individuals
                if avgDistance(ind, pop) >= distance_tr]
    return inds[:n]


def best(pop):
    ind = topn(pop, 1)[0]
    return round(avgDistance(ind, pop), 2)


def distance(chromosome, target, squared=False):
    "Calculate the (sqared) Euclidean distance of two vectors"
    euc_sq = []

    for a, b in zip(chromosome, target):
        euc_sq.append((a - b) ** 2)

    sq = sum(euc_sq)

    if squared:
        return sq

    import math
    return math.sqrt(sq)


def avgDistance(ind, pop):
    return (distance(ind.phenotype, pop.sTarget) + \
            distance(ind.phenotype, pop.rTarget)) / 2


def describe(pop, show=None):
    "Print out useful information about a population"
    distances = [avgDistance(ind, pop) for ind in pop.individuals]
    print("------------------------")
    size = len(pop.individuals)
    if (not isinstance(show, int)) or (show < 0) or (show > size):
        show = 0
    for ind in pop.individuals[:show]:
        print("Ahem")
        print(ind)
    print("Size:\t", size, sep="\t")
    print("Best fitness:", min(fitnesses), sep="\t")
    print("Avg fitness:", pop.fitness(), sep="\t")
    print("-- Population genetics descriptors --")
    print("Heterozygosity:", heterozygosity(pop), sep="\t")
    print("Gene loss:", str(round(gene_loss(pop), 2)) + "%", "\n", sep="\t")
