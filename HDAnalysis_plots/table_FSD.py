import numpy
from collections import Counter

def createTableFSD2(list1, diff=True):
    selfAB = numpy.concatenate(list1)
    uniqueFS = numpy.unique(selfAB)
    nr = numpy.arange(0, len(uniqueFS), 1)
    if diff is False:
        count = numpy.zeros((len(uniqueFS), 6))
    else:
        count = numpy.zeros((len(uniqueFS), 7))

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
                for i, l  in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 0] = j[1]
            if state == 2:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 1] = j[1]

            if state == 3:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 2] = j[1]

            if state == 4:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 3] = j[1]

            if state == 5:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 4] = j[1]

            if state == 6:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 5] = j[1]

            if state == 7:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 6] = j[1]
            state = state + 1

        sumRow = count.sum(axis=1)
        sumCol = count.sum(axis=0)

    uniqueFS = uniqueFS.astype(str)
    if uniqueFS[len(uniqueFS) - 1] == "20":
        uniqueFS[len(uniqueFS) - 1] = ">20"

    first = ["FS={}".format(i) for i in uniqueFS]
    final = numpy.column_stack((first, count, sumRow))

    return (final, sumCol)

def createFileFSD2(summary, sumCol, overallSum, output_file, name,sep, rel=False, diff=True):
    output_file.write(name)
    output_file.write("\n")
    if diff is False:
        output_file.write("{}HD=1{}HD=2{}HD=3{}HD=4{}HD=5-8{}HD>8{}sum{}\n".format(sep,sep,sep,sep,sep,sep,sep,sep))
    else:
        if rel is False:
            output_file.write("{}diff=0{}diff=1{}diff=2{}diff=3{}diff=4{}diff=5-8{}diff>8{}sum{}\n".format(sep,sep,sep,sep,sep,sep,sep,sep,sep))
        else:
            output_file.write("{}diff=0{}diff=0.1{}diff=0.2{}diff=0.3{}diff=0.4{}diff=0.5-0.8{}diff>0.8{}sum{}\n".format(sep,sep,sep,sep,sep,sep,sep,sep,sep))

    for item in summary:
        for nr in item:
            if "FS" not in nr and "diff" not in nr:
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

