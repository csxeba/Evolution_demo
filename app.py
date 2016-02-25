from tkinter import *

import Biometrics.describe as desc
from ColorEvolution import *


class App(Tk):
    def __init__(self):
        Tk.__init__(self)

        self.title("Evolúciós algoritmus -- szimuláció")

        self.calc_genotypes = BooleanVar(value=True)
        self.show_genotypes = BooleanVar(value=True)

        self.epochs = IntVar(value=0)
        self.size = IntVar(value=0)
        self.avgFitness = DoubleVar(value=0.0)
        self.bestFitness = DoubleVar(value=0.0)
        self.heterozygosity = DoubleVar(value=0.0)
        self.mutants = DoubleVar(value=0.0)

        self.log = []

        #### HEADER INITIALIZATION

        # MAIN STORY FRAME
        story = Frame(self, bd=4, relief=RAISED)
        story.pack()

        #### SCALES' INITIALIZATION
        leftFrame = Frame(story, bd=4, relief=RAISED)
        leftFrame.pack(side=LEFT)

        controls = Frame(leftFrame, bd=4, relief=RAISED)
        controls.pack()

        self.chains = (
            "Szelekciós R (lókusz 1)",
            "Szelekciós G (lókusz 2)",
            "Szelekciós B (lókusz 3)",
            "Maximális populációméret",
            "Szelektálódó egyedek (%)",
            "Reprodukálódó egyedek (%)",
            "Sz. nyomás csökk. körönként (%)",
            "Crossing over ráta (%)",
            "Mutációs ráta (%)")
        self.popvars = {c: None for c in self.chains}

        maxes = (255, 255, 255, 1000, 100, 100, 100, 100, 100)

        self.values = {c: (IntVar(), m) for c, m in zip(self.chains, maxes)}

        for c in self.chains:
            ##            rown = counter % 3
            ##            coln = int(str(counter/3)[0])
            v = self.values[c][0]
            m = self.values[c][1]
            Scale(controls, from_=0, to=m, sliderlength=30, length=200,
                  label=c, font=("Times", 8), variable=v,
                  orient=HORIZONTAL, command=self._update_vals
                  ).pack()

        # /SCALES' INITALIZATION

        rightFrame = Frame(story, bd=4, relief=RAISED)
        rightFrame.pack(side=LEFT)

        descriptors = Frame(rightFrame, bd=4, relief=RAISED)
        descriptors.pack()

        w1 = w2 = 12

        for i, c, v in zip(
                range(6),
                ("Korszak:", "Méret:", "Átl. fitnesz:",
                 "Legjobb egyed:", "Heterozigozitás:", "Mutáns arány:"),
                (self.epochs, self.size, self.avgFitness,
                 self.bestFitness, self.heterozygosity, self.mutants)):
            rown = i % 3
            coln = i % 2
            px = 15

            f = Frame(descriptors, bd=2, relief=RAISED)
            f.grid(row=rown, column=coln, padx=px)

            Label(f, text=c, width=w1, bd=2, relief=RIDGE
                  ).pack()
            Label(f, textvar=v, width=w2, bd=2, relief=RIDGE
                  ).pack()

        genPhen = Frame(rightFrame, bd=4, relief=RAISED)
        genPhen.pack()

        calcPheno = Frame(genPhen, bd=4, relief=RAISED)
        calcPheno.pack()

        Label(calcPheno, text="Fitnesz kiszámítása:").pack()

        px = 8

        Radiobutton(calcPheno, text="genotípus alapján", var=self.calc_genotypes,
                    value=True, indicatoron=0,
                    command=lambda: self._update_vals(None)
                    ).pack(side=LEFT, padx=px)
        Radiobutton(calcPheno, text="fenotípus alapján", var=self.calc_genotypes,
                    value=False, indicatoron=0,
                    command=lambda: self._update_vals(None)
                    ).pack(side=LEFT, padx=px)

        showPheno = Frame(genPhen, bd=4, relief=RAISED)
        showPheno.pack()

        Label(showPheno, text="Egyedek megjelenítése:").pack()

        Button(showPheno, text="genotípus alapján",
               command=lambda: self.change_view(True)
               ).pack(side=LEFT, padx=px)
        Button(showPheno, text="fenotípus alapján",
               command=lambda: self.change_view(False)
               ).pack(side=LEFT, padx=px)

        canframe = Frame(rightFrame, bd=4, relief=RAISED)
        canframe.pack(side=LEFT)

        self.canvas = Canvas(canframe, width=240, height=240)
        self.canvas.pack(padx=2, pady=2)

        # Buttons frame at the bottom
        buttons = Frame(self, bd=4, relief=RAISED)
        buttons.pack(fill=X)

        for c, func in zip(("+1 korszak", "+10 korszak",
                            "Szelekció", "Szaporítás", "Mutáció"),
                           (self.epoch, lambda: self.epoch(10),
                            self.selection, self.reproduction,
                            self.mutation)):
            Button(buttons, text=c, command=func, width=12
                   ).pack(side=LEFT, padx=2)

        self._init_pop()
        self.protocol("WM_DELETE_WINDOW", self.dumpLog)

    def _init_pop(self):
        self.population = Population(95)
        p = self.population
        hypers = [p.sTarget[0], p.sTarget[1], p.sTarget[2],
                  p.limit, p.diers_perc, p.reproducers_perc,
                  p.diers_decr_rate, p.crossing_over_rate,
                  p.mutation_rate]
        hypers = [int(hyper) if hyper >= 1 else int(hyper * 100)
                  for hyper in hypers]
        for c, h in zip(self.chains, hypers):
            self.values[c][0].set(h)

        self.update_canvas()

    def _update_vals(self, x):
        p = self.population
        hypers = [p.sTarget[0], p.sTarget[1], p.sTarget[2],
                  p.limit, p.diers_perc, p.reproducers_perc,
                  p.diers_decr_rate, p.crossing_over_rate,
                  p.mutation_rate]
        for i, c in enumerate(self.chains[:3]):
            p.sTarget[i] = self.values[c][0].get()

        p.limit = self.values["Maximális populációméret"][0].get()
        p.diers_perc = self.values["Szelektálódó egyedek (%)"][0].get()
        p.reproducers_perc = self.values["Reprodukálódó egyedek (%)"][0].get()
        p.crossing_over_rate = self.values["Crossing over ráta (%)"][0].get()
        p.mutation_rate = self.values["Mutációs ráta (%)"][0].get()
        p.fitness_by_gene = self.calc_genotypes

    def epoch(self, epochs=1):
        for _ in range(epochs):
            mutation(self.population)
            self.population.reproduction(0)
            self.population.selection(0)
            while len(self.population.individuals) >= self.population.limit:
                self.population.selection(0)
            self.update_canvas()
            self.epochs.set(self.epochs.get() + 1)

            self.log.append((
                self.size.get(),
                self.avgFitness.get(),
                self.bestFitness.get(),
                self.heterozygosity.get(),
                self.mutants.get(),
                desc.topn(self.population, 1)[0].phenotype))

    def selection(self):
        self.population.selection(0)
        self.update_canvas()

    def reproduction(self):
        self.population.reproduction(0)
        self.update_canvas()

    def mutation(self):
        mutation(self.population)
        self.update_canvas()

    def update_canvas(self):
        self.canvas.delete("all")
        inds = []
        for i, ind in enumerate(self.population.individuals):
            inds.append(IndRepr(self.canvas, ind, i))
            inds[-1].show(self.show_genotypes.get())
        self.size.set(len(self.population.individuals))
        self.avgFitness.set(self.population.fitness())
        self.bestFitness.set(desc.best(self.population))
        self.heterozygosity.set(desc.heterozygosity(self.population))
        self.mutants.set(desc.mutants(self.population))

    def change_view(self, x):
        self.show_genotypes.set(x)
        self.update_canvas()

    def dumpLog(self):
        chain = ''
        headers = ["Epoch", "Size", "avgFitness", "bestFitness", "Heterozygosity",
                   "Mutants", "BestColor"]
        for h in headers:
            chain = chain + h + "\t"
        chain = chain[:-1] + "\n"

        for epoch, record in enumerate(self.log):
            epoch += 1
            chain = chain + str(epoch) + "\t"
            for field in record:
                chain = chain + str(field) + "\t"
            chain = chain[:-1] + "\n"

        fl = open("run.log", mode="w")
        fl.write(chain)
        fl.close()

        self.destroy()


class IndRepr:
    def __init__(self, master, ind, no):
        self.master = master
        rown = math.trunc(no / 10)
        coln = no % 10
        size = 24
        x, y = 12 + 24 * rown, 12 + 24 * coln
        x00, y00 = x - 10, y - 10
        x01, y01 = x, y + 10
        x10, y10 = x, y - 10
        x11, y11 = x + 10, y + 10
        self.coords = ((x00, y00, x01, y01),
                       (x10, y10, x11, y11))

        self.colorA = '#%02x%02x%02x' % (ind.chromosomeA[0], ind.chromosomeA[1],
                                         ind.chromosomeA[2])
        self.colorB = '#%02x%02x%02x' % (ind.chromosomeB[0], ind.chromosomeB[1],
                                         ind.chromosomeB[2])
        self.colorPh = '#%02x%02x%02x' % (ind.phenotype[0], ind.phenotype[1],
                                          ind.phenotype[2])

    def show(self, geno):
        if geno:
            c1 = self.colorA
            c2 = self.colorB
            self.master.create_rectangle(*self.coords[0], fill=c1)
            self.master.create_rectangle(*self.coords[1], fill=c2)

        else:
            c = self.colorPh
            self.master.create_rectangle(self.coords[0][0], self.coords[0][1],
                                         self.coords[1][2], self.coords[1][3],
                                         fill=self.colorPh)


if __name__ == '__main__':
    app = App()
    app.mainloop()
