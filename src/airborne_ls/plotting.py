"""
Some simple plotting tools. Currently implements only ASCII-based plotting, for
very simple plots without introducing a dependency on matplotlib or similar.
"""
import numpy as np


def asciiBoxplot(values, title, leftMargin=5, width=80):
    """
    Print an ASCII character boxplot of the given values.

    The axis values are unlikely to work well if everything is either very large or
    very small, since these give rise to many digits, which won't fit neatly along
    the scale. With this in mind, labels are only shown as integers, so don't bother if
    the data is all < 1.

    Parameters:
      values (numpy array): Data whose distribution is to be plotted
      title (str): If not None, a title to print above the plot
      leftMargin (int): Number of characters blank at the left of the plot
      width (int): Number of characters occupied by the plot itself
    """
    values = np.array(values)
    median = np.median(values)
    lq = np.percentile(values, 25)
    uq = np.percentile(values, 75)
    lowWhisker = np.percentile(values, 5)
    highWhisker = np.percentile(values, 95)
    minVal = values.min()
    maxVal = values.max()

    # Re-scale the critical values to the plot coordinates (i.e. characters across the screen)
    critVals = np.array([minVal, lowWhisker, lq, median, uq, highWhisker, maxVal])
    critVals_char = (leftMargin + width * (critVals - minVal) / (maxVal - minVal)).astype(np.int32)

    if title is not None:
        titleLen = len(title)
        padding = leftMargin + (width - titleLen) // 2
        titleLine = ' ' * padding + title
        print(titleLine)

    # Construct the boxplot itself. Min & max are 'o', whiskers are '-',
    # box is '=' and median is '|'
    plotLine = np.array([' '] * (width + leftMargin + 1))
    plotLine[critVals_char[0]] = 'o'
    plotLine[critVals_char[1]:critVals_char[2]] = '-'
    plotLine[critVals_char[2]:critVals_char[4]] = '='
    plotLine[critVals_char[3]] = '|'
    plotLine[critVals_char[4]:critVals_char[5]] = '-'
    plotLine[critVals_char[6]] = 'o'
    plotLineStr = ''.join(plotLine)
    print(plotLineStr)

    # Construct the value labels line
    labelsLine = np.array([' '] * (width + 2 * leftMargin))
    # We label the min & max, the quartiles, and the median
    for i in [0, 2, 3, 4, 6]:
        val = critVals[i]
        pos = critVals_char[i]
        valStr = str(int(round(val)))
        posLeft = pos - len(valStr) // 2
        labelsLine[posLeft:posLeft + len(valStr)] = [c for c in valStr]
    labelsLineStr = ''.join(labelsLine)
    print(labelsLineStr)
