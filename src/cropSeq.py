#! /usr/bin/env python

def determineStartEnd(sequence, mingap):
    '''
    Determines the start and the end of a sequence

    Parameters
    ----------
    sequence: array
        sequence

    Returns
    -------
    todo
        value for start and end of sequence
    '''
    start = 0
    end = 0

    start = findValue(sequence, mingap)
    end = len(sequence) - findValue(sequence[::-1])

    # todo: look into this!
    if start > end:
        return (0, 0)
    return(start, end)

def findValue(sequence, mingap=10):
    '''
    Determines the start of the given sequence

    Parameters
    ----------
    sequence: array
        sequence

    Returns
    -------
    int
        value for start or end of sequence
    '''

    position = 0
    boundary1 = 50
    boundary2 = 80
    boundary3 = 20

    # todo: make boundarys parameters!

    gaps = countGaps(sequence)

    if len(gaps) < 11:
        return(gaps[0] + 1)

    if len(gaps) <= 80:
        boundary1 = 10
        boundary2 = 19
        boundary3 = 10


    # this pattern doesn't indicate an incomplete sequence, set start to 0
    if gaps[boundary1] < boundary3:
        return 0

    # for more fluctuation within the sequence, meaning we observe a few nt within many gaps -> indicates incomplete sequence
    for n in range(0, boundary2):
        if gaps[n+1] - gaps[n] > mingap:
            position = n + 1 + gaps[n+1]
    if position > 0:
        return position

    # if none of above, take the first nt/aa as start of the sequence
    return(gaps[0] + 1)

def countGaps(sequence):
    '''
    Counts the gaps in a sequence for each non-gap position

    Parameters
    ----------

    sequence: array

    Returns
    -------
    todo
        array of ints of the length of the sequence without gaps

    '''

    gapNumbers = []
    gapCounter = 0

    for element in sequence:
        if element == '-':
            gapCounter += 1
        else:
            gapNumbers.append(gapCounter)

    return gapNumbers

