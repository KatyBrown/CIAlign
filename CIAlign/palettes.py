#!/usr/bin/env python3


def base():
    '''
    Returns the hexadecimal values for black and white

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for black and white
    '''

    return {'black': '#000000',
            'white': '#FFFFFF'}


def CBSafe():
    '''
    Returns the hexadecimal values for a colour blind safe colour palette

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for the colours used
        in the CIAlign mini alignments
    '''

    b = base()
    b.update({'yellow_nt': "#c9c433",
              'green_nt': "#56ae6c",
              'red_nt': "#a22c49",
              'blue_nt': "#0038a2",
              'grey_nt': "#6979d3",
              'red_aa': "#a22c49",
              'yellow_aa': "#c9c433",
              'blue_aa':  "#0038a2",
              'orange_aa': "#e57700",
              'midblue_aa': "#589aab",
              'cyan_aa': "#50d3cb",
              'lightgrey_aa': '#eae2ea',
              'green_aa': "#56ae6c",
              'darkgrey_aa': "#888988",
              'purple_aa': '#89236a',
              'paleblue_aa': '#e669ca',
              'peach_aa': "#ffc4a9",
              'tan_aa': "#936e23",
              'remove_insertions': "#9db341",
              'remove_divergent': "#7066bc",
              'crop_ends': '#020545',
              'remove_gap_only': '#f9c1d2',
              'remove_short': "#c85133",
              'crop_divergent': '#ff00d1',
              'crop_ltrs': '#1e7507',
              'crop_tirs': '#ffb900',
              'user': '#fff4a3'})
    return (b)


def Bright():
    '''
    Returns the hexadecimal values for a colour blind safe colour palette

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for the colours used
        in the CIAlign mini alignments
    '''
    b = base()
    b.update({'yellow_nt': "#ffd500",
              'red_nt': "#64bc3c",
              'green_nt': "#f20707",
              'blue_nt': "#0907f2",
              'grey_nt': "#c7d1d0",
              'red_aa': "#f20707",
              'yellow_aa': "#ffd500",
              'blue_aa':  "#0907f2",
              'orange_aa': "#f4aa03",
              'midblue_aa': "#03b5f4",
              'cyan_aa': "#03f4dd",
              'lightgrey_aa': '#f0f4f4',
              'green_aa': "#64bc3c",
              'darkgrey_aa': "#808080",
              'purple_aa': '#6f0cac',
              'paleblue_aa': '#cde3f8',
              'peach_aa': "#f8e7cd",
              'tan_aa': "#d2a867",
              'remove_insertions': "#9db341",
              'remove_divergent': "#7066bc",
              'crop_ends': '#020545',
              'remove_gaponly': '#f9c1d2',
              'remove_short': "#c85133",
              'crop_divergent': '#ff00d1',
              'crop_ltrs': '#1e7507',
              'crop_tirs': '#ffb900',
              'user': '#fff4a3'})
    return (b)
