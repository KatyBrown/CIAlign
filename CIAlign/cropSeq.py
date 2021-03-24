#! /usr/bin/env python


def determineStartEnd(sequence, name, log, mingap_perc=0.05, redefine_perc=0.1):
    '''
    Determines the start and the end of a sequence by calling a subroutine

    Parameters
    ----------
    sequence: np.array
        sequence

    name: string
        name of the sequence

    log: logging.Logger
        An open log file object

    mingap_perc: float
        proportion of the sequence length (excluding gaps) that is the
        threshold for change in gap numbers within the first 10% of non-gap
        positions of the sequence \
        (default: 0.05)

    redefine_perc: float
        proportion of the sequence length (excluding gaps) that we check for
        change in gap numbers to redefine start/end
        (default: 0.1)

    Returns
    -------
    start: int
        redefined start of the sequence

    end: int
        redefined end of the seuquence

    '''

    start = 0
    end = 0
    gaps = countGaps(sequence)
    # make sure the mingap is not lower than 2, and if it is give a warning
    mingap = int(mingap_perc*len(gaps))
    # print("mingap:", mingap)
    if mingap < 2 and len(gaps) > 20:
        log.warning("Given the length of sequence %s, the mingap_perc threshold is too low and the change in gap numbers was replaced by 2." % name)
        print("Given the length of sequence %s, the mingap_perc threshold is too low and the change in gap numbers was replaced by 2.\n" % name)

    start = findValue(sequence, mingap_perc, redefine_perc)
    # put in reverse for end
    end = len(sequence) - findValue(sequence[::-1], mingap_perc, redefine_perc)

    if start > end:
        return (0, 0)
    return(start, end)


def findValue(sequence, mingap_perc=0.05, redefine_perc=0.1):
    '''
    Determines the start of the given sequence

    Parameters
    ----------
    sequence: numpy array
        sequence

    name: string
        name of the sequence

    log: logging.Logger
        An open log file object

    mingap_perc: float
        proportion of the sequence length (excluding gaps) that is the
        threshold for change in gap numbers within the first checkfor_perc %
        of non-gap positions of the sequence
        (default: 0.05)

    redefine_perc: float
        proportion of the sequence length (excluding gaps) that we check for
        change in gap numbers to redefine start/end
        (default: 0.1)

    Returns
    -------
    int
        value for start or end (when put in reverse) of sequence
    '''

    position = 0

    gaps = countGaps(sequence)

    seq_length = len(gaps)
    boundary1 = int(0.1*seq_length)
    # threshold for how many non-gap positions we look at
    boundary2 = int(redefine_perc*seq_length)
    boundary3 = int(0.01*seq_length)
    # the threshold for the change in gap numbers
    mingap = int(mingap_perc*seq_length)
    # make sure mingap is not lower than 2 and if it is, set it to 2
    if mingap < 2:
        mingap = 2

    # for very short sequences it is not desirable to redefine
    if seq_length < 21:
        return(gaps[0] + 1)

    # this pattern doesn't indicate an incomplete sequence, set start to 0
    if gaps[boundary1] < boundary3:
        return 0

    # for more fluctuation within the sequence, meaning we observe a few nt
    # within many gaps -> indicates incomplete sequence
    for n in range(0, boundary2):
        if gaps[n+1] - gaps[n] >= mingap:
            position = n + 1 + gaps[n+1]
    if position > 0:
        return position

    # if none of above, take the first nt/aa as start of the sequence
    return(gaps[0] + 1)


def countGaps(sequence):
    '''
    Counts the gaps in a sequence preceding each non-gap position

    Parameters
    ----------

    sequence: numpy array
        sequence

    Returns
    -------
    gapNumbers: list
        list of ints of the length of the sequence without gaps

    '''

    gapNumbers = []
    gapCounter = 0

    for element in sequence:
        if element == '-':
            gapCounter += 1
        else:
            gapNumbers.append(gapCounter)

    return gapNumbers
