import matplotlib.pyplot as plt
import numpy
import os
import time

def plotHDwithinSeq_Sum2(sum1, sum2,min_value, lenTags, title_file1, pdf):
    fig = plt.figure(figsize=(6, 8))
    plt.subplots_adjust(bottom=0.1)

    ham = [numpy.array(min_value), sum1, sum2]  # new hd within tags

    maximumX = numpy.amax(numpy.concatenate(ham))
    minimumX = numpy.amin(numpy.concatenate(ham))
    maximumY = numpy.amax(numpy.concatenate(map(lambda (x): numpy.bincount(x), ham)))

    if len(range(minimumX, maximumX)) == 0:
        range1 = minimumX
    else:
        range1 = range(minimumX, maximumX + 2)

    counts = plt.hist(ham, align="left", rwidth=0.8, stacked=False,
                      label=["HD of whole tag", "tag1 - a\nvs. tag2 - a", "tag1 - b\nvs. tag2 - b"],
                      bins=range1, color=["#585858", "#58ACFA", "#FA5858"], edgecolor='black', linewidth=1)
    plt.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(1.55, 1))
    plt.suptitle('Hamming distances within tags', fontsize=14)
    plt.title(title_file1, fontsize=12)
    plt.xlabel("Hamming Distance", fontsize=12)
    plt.ylabel("Absolute Frequency", fontsize=12)
    plt.grid(b=True, which='major', color='#424242', linestyle=':')


    plt.axis((minimumX - 1, maximumX + 1, 0, maximumY * 1.1))
    plt.xticks(numpy.arange(minimumX - 1, maximumX + 1, 1.0))
    plt.ylim((0, maximumY * 1.1))

    legend = "sample size= {:,} against {:,}".format(len(ham[0]), lenTags, lenTags)
    plt.text(0.14, -0.01, legend, size=12, transform=plt.gcf().transFigure)
    pdf.savefig(fig, bbox_inches="tight")

#        plt.savefig("{}_within Seq_{}.png".format(title_savedFile, name), bbox_inches="tight")

    plt.close("all")
    plt.clf()



