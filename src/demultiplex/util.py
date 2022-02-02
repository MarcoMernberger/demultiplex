import cutadapt
import cutadapt.align
import collections
import subprocess
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pypipegraph as ppg
import numpy as np
from mplots import MPPlotJob
from .demultiplex import Demultiplexer

AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


def locate(adapter_sequence, maximal_error_rate=1, start=True):
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
    if gz_filename1.endswith(".gz"):
        p1 = subprocess.Popen(["gunzip", "-c", gz_filename1], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
        x = int(p2.communicate()[0][:-1])
        x = x / 4
        return x
    else:
        t = subprocess.check_output(["wc", "-l", gz_filename1])
        if t == b"":
            return 0
        else:
            x = int(t.split()[0])
            x = x / 4
            return x


def plot_read_counts(demultiplexer_or_list, outfile, dependencies=[], log=True):
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
        tmp = {"Sample": [], "Count": []}
        for sample in samples_to_plot:
            read_count = count_raw_input_reads(
                str(sample.get_aligner_input_filenames()[0])
            )
            tmp["Sample"].append(sample.name)
            tmp["Count"].append(read_count)
        df = pd.DataFrame(tmp)
        df["Color"] = color
        df.to_csv(outfile_df, sep="\t", index=False)

        fig = plt.figure(figsize=(12, 12))
        x = np.arange(len(samples_to_plot))
        plt.bar(x, df["Count"].values, width=0.8, color=df["Color"].values)
        labels = [str(x.strip()) for x in df["Sample"].values]
        plt.xticks(ticks=x, labels=labels, rotation=75, ha="right")
        if log:
            plt.yscale("log")
        plt.title("Demultiplexed read counts")
        plt.tight_layout()
        fig.savefig(outfile)

    return ppg.MultiFileGeneratingJob([outfile, outfile_df], __count).depends_on(
        dependencies + deps
    )


def plot_read_counts_percentage(
    demultiplexer_or_list, outfile, dependencies=[], log=True
):
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
        samples_to_plot = {}
        color = []
        tmp = {"Sample": [], "Percentage Reads": [], "Color": []}
        for dm in demultiplexer:
            samples_to_plot[dm.input_sample] = [dm.input_sample]
            read_count_100 = count_raw_input_reads(
                str(dm.input_sample.get_aligner_input_filenames()[0])
            )
            tmp["Color"].append("g")
            tmp["Sample"].append(dm.input_sample.name)
            tmp["Percentage Reads"].append(100.0)
            dm_samples = dm.get_samples()
            for dm_sample_name in dm_samples:
                sample = dm_samples[dm_sample_name]
                read_count = count_raw_input_reads(
                    str(sample.get_aligner_input_filenames()[0])
                )
                perc = 100 * (float(read_count) / read_count_100)
                tmp["Color"].append("c")
                tmp["Sample"].append(sample.name)
                tmp["Percentage Reads"].append(perc)
        df = pd.DataFrame(tmp)
        return df

    def __plot(df):
        fig = plt.figure(figsize=(12, 12))
        x = np.arange(df.shape[0])
        plt.barh(x, df["Percentage Reads"].values, height=0.8, color=df["Color"].values)
        for i, v in enumerate(df["Percentage Reads"].values):
            plt.gca().text(
                v + 2,
                i - 0.3,
                f"{v:.2f}%",
                color=df["Color"].values[i],
            )
        labels = [str(x.strip()) for x in df["Sample"].values]
        plt.yticks(ticks=x, labels=labels, ha="right")
        plt.xlim([0, 115])
        plt.title("Demultiplexed read count percentages")
        plt.tight_layout()
        return fig

    return MPPlotJob(outfile, __count, __plot).depends_on(dependencies + deps)