import matplotlib.pyplot as plt
import numpy
from collections import Counter
import os
import time

def plotHDwithFSD(list1,maximumX,minimumX, subtitle, lenTags, title_file1,pdf,
                   xlabel,relative=False):
    if relative is True:
        step = 0.1
    else:
        step = 1

    fig = plt.figure(figsize=(6, 8))
    plt.subplots_adjust(bottom=0.1)
    con_list1 = numpy.concatenate(list1)
    p1 = numpy.array([v for k, v in sorted(Counter(con_list1).iteritems())])
    maximumY = numpy.amax(p1)

    if relative is True:  # relative difference
        bin1 = numpy.arange(-1, maximumX + 0.2, 0.1)
    else:
        bin1 = maximumX + 1

    counts = plt.hist(list1, bins=bin1, edgecolor='black', linewidth=1,
                      label=["FS=1", "FS=2", "FS=3", "FS=4", "FS=5-10",
                             "FS>10"], rwidth=0.8,
                      color=["#808080", "#FFFFCC", "#FFBF00", "#DF0101", "#0431B4", "#86B404"],
                      stacked=True, alpha=1,
                      align="left",
                      range=(0, maximumX + 1))
    plt.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(1.45, 1))
    bins = counts[1]  # width of bins
    counts = numpy.array(map(int, counts[0][5]))
    plt.suptitle(subtitle, y=1, x=0.5, fontsize=14)
    plt.title(title_file1, fontsize=12)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel("Absolute Frequency", fontsize=12)

    plt.grid(b=True, which='major', color='#424242', linestyle=':')
    plt.axis((minimumX - step, maximumX + step, 0, numpy.amax(counts) + sum(counts) * 0.1))
    plt.xticks(numpy.arange(0, maximumX + step, step))

    plt.ylim((0, maximumY * 1.2))

    bin_centers = -0.4 * numpy.diff(bins) + bins[:-1]
    for x_label, label in zip(counts, bin_centers):  # labels for values
        if x_label == 0:
            continue
        else:
            plt.annotate("{:,}\n{:.3f}".format(x_label, float(x_label) / sum(counts), 1),
                         xy=(label, x_label + len(con_list1) * 0.01),
                         xycoords="data", color="#000066",fontsize=10)

    legend = "sample size= {:,} against {:,}".format(sum(counts), lenTags)
    plt.text(0.14, -0.01, legend, size=12, transform=plt.gcf().transFigure)

    pdf.savefig(fig, bbox_inches="tight")
#        plt.savefig("{}_HDwithFSD_{}.png".format(title_savedFile, name), bbox_inches="tight")

  #  print('File as "{}_HDwithFSD_{}.png" saved in your home directory!'.format(title_savedFile, name))
    plt.close("all")

    plt.clf()
