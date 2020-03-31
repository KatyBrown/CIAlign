#! /usr/bin/env python

from time import time
import psutil
import os

def start_mem_time(name):

    # create dict
    start_values = {}

    start_values["name"] = name

    #Init setzen
    start_values["t"] = -time()
    start_values["base_mem"] = memory_usage_psutil()

    return start_values

def end_mem_time(start_values, file):
    # start_values is a dict

    # Div ist damit automatisch bestimmt
    total_time = start_values["t"] + time()

    new_mem = memory_usage_psutil()
    alloc_mem = new_mem - start_values["base_mem"]

    # write to file
    out = open(file, "a")
    out.write("Function: " + start_values["name"] + "\n")
    out.write("time: " + str(total_time) + " s\n")
    out.write("memory: " + str(alloc_mem) + " MB\n\n")
    out.close()


def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)
    return mem
