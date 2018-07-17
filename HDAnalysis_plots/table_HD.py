import numpy
from collections import Counter
def createTableHD(list1, row_label):
    selfAB = numpy.concatenate(list1)
    uniqueHD = numpy.unique(selfAB)
    nr = numpy.arange(0, len(uniqueHD), 1)
    count = numpy.zeros((len(uniqueHD), 6))
    state = 1
    for i in list1:
        counts = list(Counter(i).items())
        hd = [item[0] for item in counts]
        c = [item[1] for item in counts]
        table = numpy.column_stack((hd, c))
        if len(table) == 0:
            state = state + 1
            continue
        else:
            if state == 1:
                for i, l  in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 0] = j[1]
            if state == 2:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 1] = j[1]

            if state == 3:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 2] = j[1]

            if state == 4:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 3] = j[1]

            if state == 5:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 4] = j[1]

            if state == 6:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 5] = j[1]
            state = state + 1

        sumRow = count.sum(axis=1)
        sumCol = count.sum(axis=0)
        first = ["{}{}".format(row_label,i) for i in uniqueHD]
        final = numpy.column_stack((first, count, sumRow))

    return (final, sumCol)

def createTableHDwithTags(list1):
    selfAB = numpy.concatenate(list1)
    uniqueHD = numpy.unique(selfAB)
    nr = numpy.arange(0, len(uniqueHD), 1)
    count = numpy.zeros((len(uniqueHD), 3))

    state = 1
    for i in list1:
        counts = list(Counter(i).items())
        hd = [item[0] for item in counts]
        c = [item[1] for item in counts]
        table = numpy.column_stack((hd, c))
        if len(table) == 0:
            state = state + 1
            continue
        else:
            if state == 1:
                for i, l  in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 0] = j[1]
            if state == 2:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 1] = j[1]

            if state == 3:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 2] = j[1]
            state = state + 1

        sumRow = count.sum(axis=1)
        sumCol = count.sum(axis=0)
        first = ["HD={}".format(i) for i in uniqueHD]
        final = numpy.column_stack((first, count, sumRow))

    return (final, sumCol)


def createFileHD(summary, sumCol, overallSum, output_file, name,sep):
    output_file.write(name)
    output_file.write("\n")
    output_file.write("{}FS=1{}FS=2{}FS=3{}FS=4{}FS=5-10{}FS>10{}sum{}\n".format(sep,sep,sep,sep,sep,sep,sep,sep))
    for item in summary:
        for nr in item:
            if "HD" not in nr and "diff" not in nr:
                nr = nr.astype(float)
                nr = nr.astype(int)
            output_file.write("{}{}".format(nr,sep))
        output_file.write("\n")
    output_file.write("sum{}".format(sep))
    sumCol = map(int, sumCol)
    for el in sumCol:
        output_file.write("{}{}".format(el,sep))
    output_file.write("{}{}".format(overallSum.astype(int),sep))
    output_file.write("\n\n")

def createFileHDwithinTag(summary, sumCol, overallSum, output_file, name,sep):
    output_file.write(name)
    output_file.write("\n")
    output_file.write("{}HD of whole tag;tag1-half1 vs. tag2-half1{}tag1-half2 vs. tag2-half2{}sum{}\n".format(sep,sep,sep,sep))
    for item in summary:
        for nr in item:
            if "HD" not in nr:
                nr = nr.astype(float)
                nr = nr.astype(int)
            output_file.write("{}{}".format(nr,sep))
        output_file.write("\n")
    output_file.write("sum{}".format(sep))
    sumCol = map(int, sumCol)
    for el in sumCol:
        output_file.write("{}{}".format(el,sep))
    output_file.write("{}{}".format(overallSum.astype(int),sep))
    output_file.write("\n\n")
