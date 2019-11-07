import cutadapt
import cutadapt.align
import collections
from .demultiplex import Demultiplexer
import subprocess
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pypipegraph as ppg
import numpy as np

AdapterMatch = collections.namedtuple(
        "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
) 


def locate(adapter_sequence, maximal_error_rate = 1, start = True):
    if start:
        where = (
            cutadapt.align.START_WITHIN_SEQ1
            | cutadapt.align.START_WITHIN_SEQ2
            | cutadapt.align.STOP_WITHIN_SEQ2
        )
    else:
        where = (
            cutadapt.align.STOP_WITHIN_SEQ1
            | cutadapt.align.START_WITHIN_SEQ2
            | cutadapt.align.STOP_WITHIN_SEQ2
        )
    adapter = cutadapt.align.Aligner(
        adapter_sequence,
        maximal_error_rate / len(adapter_sequence),
        where,
        wildcard_ref=True,
        wildcard_query=False,
        )
    def match(seq):
        alignment = adapter.locate(seq)
        if alignment is None:
            return None
        _match = AdapterMatch(*alignment)
        print(_match)
        return _match
    return match

def count_raw_input_reads(gz_filename1):
    if gz_filename1.endswith('.gz'):
        p1 = subprocess.Popen(['gunzip', '-c', gz_filename1], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
        x = int(p2.communicate()[0][:-1])
        x = x/4
        return x
    else:
        t = subprocess.check_output(['wc', '-l', gz_filename1])
        if t == b"":
            return 0
        else:
            x = int(t.split()[0])
            x = x/4
            return x

def plot_read_counts(demultiplexer_or_list, outfile, dependencies = [], log = True):
    outfile_df = outfile + ".tsv"
    demultiplexer = demultiplexer_or_list
    deps = []
    if isinstance(demultiplexer_or_list, Demultiplexer):
        demultiplexer = [demultiplexer]
    for dm in demultiplexer:
        deps.append(dm.do_demultiplex())
        deps.append(dm.input_sample.prepare_input())
        for sample in dm.get_samples().values():
            deps.append(sample.prepare_input())
    def __count():
        samples_to_plot = []
        color = []
        for dm in demultiplexer:
            samples_to_plot.append(dm.input_sample)
            color.append("g")
            dm_samples = dm.get_samples()
            for dm_sample_name in dm_samples:
                samples_to_plot.append(dm_samples[dm_sample_name])
                color.append("c")
        tmp = {'Sample' : [], 'Count' : []}
        for sample in samples_to_plot:
            read_count = count_raw_input_reads(str(sample.get_aligner_input_filenames()[0]))
            tmp["Sample"].append(sample.name)
            tmp['Count'].append(read_count)
        df = pd.DataFrame(tmp)
        df["Color"] = color
        df.to_csv(outfile_df, sep = "\t", index = False)

        fig = plt.figure(figsize=(12, 12))
        x = np.arange(len(samples_to_plot))
        plt.bar(x, df['Count'].values, width = .8, color = df['Color'].values)
        labels = [str(x.strip()) for x in df["Sample"].values]
        plt.xticks(ticks = x, labels = labels, rotation = 75, ha = "right")
        if log:
            plt.yscale('log')
        plt.title("Demultiplexed read counts")
        plt.tight_layout()
        fig.savefig(outfile)

    return ppg.MultiFileGeneratingJob([outfile, outfile_df], __count).depends_on(dependencies + deps)