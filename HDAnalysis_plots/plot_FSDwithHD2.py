import matplotlib.pyplot as plt
import numpy
import os
import time

def plotFSDwithHD2(familySizeList1,maximumXFS,minimumXFS, quant,
                   title_file1, subtitle, pdf, relative=False, diff = True):
    if diff is False:
        colors = ["#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4"]
        labels = ["HD=1", "HD=2", "HD=3", "HD=4", "HD=5-8","HD>8"]
    else:
        colors = ["#93A6AB", "#403C14", "#731E41", "#BAB591", "#085B6F", "#E8AA35", "#726C66"]
        if relative is True:
            labels = ["d=0", "d=0.1", "d=0.2", "d=0.3", "d=0.4", "d=0.5-0.8", "d>0.8"]
        else:
            labels = ["d=0","d=1", "d=2", "d=3", "d=4", "d=5-8","d>8"]

    fig = plt.figure(figsize=(6, 7))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.1)
    p1 = numpy.bincount(numpy.concatenate((familySizeList1)))
    maximumY = numpy.amax(p1)

    if len(range(minimumXFS, maximumXFS)) == 0:
        range1 = range(minimumXFS - 1, minimumXFS + 2)
    else:
        range1 = range(0, maximumXFS + 2)
    counts = plt.hist(familySizeList1, label=labels,
                      color=colors, stacked=True,
                      rwidth=0.8,alpha=1, align="left",
                      edgecolor="None",bins=range1)
    plt.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(1.45, 1))

    plt.title(title_file1, fontsize=12)
    plt.suptitle(subtitle, y=1, x=0.5, fontsize=14)
    plt.xlabel("No. of Family Members", fontsize=12)
    plt.ylabel("Absolute Frequency", fontsize=12)

    ticks = numpy.arange(0, maximumXFS + 1, 1)
    ticks1 = map(str, ticks)
    if maximumXFS >= 20:
        ticks1[len(ticks1) - 1] = ">=20"
    plt.xticks(numpy.array(ticks), ticks1)
    [l.set_visible(False) for (i, l) in enumerate(ax.get_xticklabels()) if i % 5 != 0]

    plt.xlim((0, maximumXFS + 1))
    if len(numpy.concatenate(familySizeList1)) != 0:
        plt.ylim((0, max(numpy.bincount(numpy.concatenate(familySizeList1))) * 1.1))

    plt.ylim((0, maximumY * 1.2))
    legend = "\nmax. family size: \nabsolute frequency: \nrelative frequency: "
    plt.text(0.15, -0.08, legend, size=12, transform=plt.gcf().transFigure)

    count = numpy.bincount(quant)  # original counts
    legend1 = "{}\n{}\n{:.5f}" \
        .format(max(quant), count[len(count) - 1], float(count[len(count) - 1]) / sum(count))
    plt.text(0.5, -0.08, legend1, size=12, transform=plt.gcf().transFigure)
    legend3 = "singletons\n{:,}\n{:.5f}".format(int(counts[0][len(counts[0]) - 1][1]),
                                                float(counts[0][len(counts[0]) - 1][1]) / sum(
                                                    counts[0][len(counts[0]) - 1]))
    plt.text(0.7, -0.08, legend3, transform=plt.gcf().transFigure, size=12)

    plt.grid(b=True, which='major', color='#424242', linestyle=':')

    #plt.savefig("{}_FSDwithHD_{}.png".format(title_savedFile, name), bbox_inches="tight")
    pdf.savefig(fig, bbox_inches="tight")

    plt.close("all")