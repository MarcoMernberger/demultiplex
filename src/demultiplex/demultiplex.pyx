"""
Created on August 02, 2019

@author: mernberger
"""
import cutadapt
import pypipegraph as ppg
import itertools
import gzip
import os
from pathlib import Path
import mbf_align
from cpython cimport bool
from mbf_align import build_fastq_strategy
import pandas as pd 
import collections
import shutil
import functools

cdef:
    struct Read:
        char* seq
        char* name
        char* qual

import cutadapt.align
from collections import namedtuple

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans

rev_comp_table = maketrans(
    b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR"
)

AdapterMatch = collections.namedtuple(
        "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
) 

def locate(adapter_sequence, maximal_error_rate = 1, start = True):
    if start:
        where = (
            #cutadapt.align.START_WITHIN_SEQ1
            cutadapt.align.START_WITHIN_SEQ2
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

def cutadapt_locate(char* seq, char* adapter, float max_error_rate, bool start):
    """use cutadapt to find an adapter (or barcode) in read, returns the start position (0 based) or -1"""
    cdef int pos
    cdef int rstart
    cdef int rstop
    #can we find the exact sequence?
    pos = seq.find(adapter)
    if pos >= 0:
        return pos
    #if not, use cutadapt to find the best match
    elif max_error_rate == 0:
        return -1
    else:
        locate_adapter = locate(adapter, max_error_rate, start)        
        match = locate_adapter(seq)
        if match is None:
            return -1 
        if start:
            return match.rstop
        else:
            return match.rstart

def hamming(char* seq1, char* seq2):
    cdef int l
    l = len(seq1)
    if l != len(seq2):
        raise ValueError("Strings Must be of the same length")
    cdef dist
    dist = 0
    cdef int ii
    for ii in range(l):
        if seq1[ii] != seq2[ii]:
            dist += 1
    return dist


def hamming_ignore_N(seq1, seq2, int max_errors):
    cdef int l
    l = len(seq1)
    if l != len(seq2):
        raise ValueError("Strings Must be of the same length")
    cdef dist
    dist = 0
    cdef int ii
    for ii in range(l):
        if seq1[ii] == 'N' or seq2[ii] == 'N':
            continue
        elif seq1[ii] != seq2[ii]:
            dist += 1
            if dist > max_errors:
                return dist
    return dist

def with_primers(adapter_sequence_begin, adapter_sequence_end):
    adapter_sequence_begin_reverse = adapter_sequence_begin[::-1].translate(rev_comp_table)
    adapter_sequence_end_reverse = adapter_sequence_end[::-1].translate(rev_comp_table)
    def filter(func):
        @functools.wraps(func)
        def wrapper(*args):
            value = func(*args)
            return value
        return wrapper
    return filter

#######################################################
def get_decision_callback(adapter_sequence_begin, adapter_sequence_end, max_error_rate, paired, forward_only):
    adapter_sequence_begin_reverse = adapter_sequence_begin[::-1].translate(rev_comp_table)
    adapter_sequence_end_reverse = adapter_sequence_end[::-1].translate(rev_comp_table)
    
#    if start and stop given:
#        if require_both_primers:
#            if require_both_primer_in_both_reads:
#                accept if all indices are not -1
#            else:
#                accept if start_fwd and stop_fwd are not -1, trim both reverses if present
#        else:
#            accept if start_fwd is not -1, trim all others if present
#    else:
#        accept if start_fwd is in read_1, trim before start_rev in seq2


    #when do we accept?
    #both reads are relevant  --> find start_fwd in seq1 and stop_rev in seq2, trim stop_fwd in seq1 and start_rev in seq2
    #                         --> if not accept missing stop, discard if stop_fwd is missing in seq1 and start_rev in seq2
    #                         -->>> need to find them stop_rev, start_rev, stop_fwd, start_fwd
    #only start is relevant   --> find start_fwd in seq1, trim start_rev in seq2
    #                         --> if stop is given, find and trim
    #                         -->>> find start_fwd, start_rev, stop_fwd and stop_r4everse only if given

    #we dont need this
    def __filter_func_start_stop_both_relevant_and_present(tup):
        seq1, qual1, name1 = tup[0]
        seq2, qual2, name2 = tup[1]
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_begin, max_error_rate, True)        
        if read1_i1 == -1 and read2_i1 == -1:
        # forward adapter nowhere to be found, discard
            return False, None, None
        elif read1_i1 == -1 and read2_i1 != -1:
            #adapter in read 2, switch
            seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            read1_i1 = read2_i1
        elif read1_i1 != -1 and read2_i1 != -1:
            raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True)
        if read2_i1 ==-1:
            #stop primer not found, discard
            return False, None, None
        read1_i2 = cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False)
        if read1_i2 == -1:
            return False, None, None
        read2_i2 = cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False)
        if read2_i2 == -1:
            return False, None, None
        return True, (read1_i1,  read1_i2), (read1_i1,  read1_i2)

    #forward reverse and paired
    def __filter_func_paired_both_relevant(tup):
        seq1, qual1, name1 = tup[0]
        seq2, qual2, name2 = tup[1]
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_begin, max_error_rate, True)        
        if read1_i1 == -1 and read2_i1 == -1:
        # forward adapter nowhere to be found, discard
            return False, None, None
        elif read1_i1 == -1 and read2_i1 != -1:
            #adapter in read 2, switch
            seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            read1_i1 = read2_i1
        elif read1_i1 != -1 and read2_i1 != -1:
            raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True)
        if read2_i1 ==-1:
            #stop primer not found, discard
            return False, None, None
        #accept, both primer found
        read1_i2 = max(0, cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False))
        read2_i2 = max(0, cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False))
        return True, (read1_i1,  read1_i2), (read1_i1,  read1_i2)

    #forward only and paired
    def __filter_func_paired_start_relevant(tup):            
        seq1, qual1, name1 = tup[0]
        seq2, qual2, name2 = tup[1]
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_begin, max_error_rate, True)        
        if read1_i1 == -1 and read2_i1 == -1:
        # forward adapter nowhere to be found, discard
            return False, None, None
        elif read1_i1 == -1 and read2_i1 != -1:
            #adapter in read 2, switch
            seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            read1_i1 = read2_i1
        elif read1_i1 != -1 and read2_i1 != -1:
            raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        #accept already
        read2_i1 = max(0, cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True))
        read1_i2 = max(0, cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False))
        read2_i2 = max(0, cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False))
        return True, (read1_i1,  read1_i2), (read1_i1,  read1_i2)

    #forward reverse and single
    def __filter_func_single_both_relevant(read):
        seq1, qual1, name1 = read
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        if read1_i1 == -1:
            return False, None
        read1_i2 = max(0, cutadapt_locate(seq1, adapter_sequence_end, max_error_rate, False))
        return True, (read1_i1,  read1_i2)

    #forward only and single
    def __filter_func_single_start_relevant(read):
        seq1, qual1, name1 = read
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        if read1_i1 == -1:
            return False, None
        return True, (read1_i1, len(seq1))

    if paired:
        if forward_only:
            filter_func = __filter_func_paired_start_relevant
        else:
            filter_func = __filter_func_paired_both_relevant
    else:
        if forward_only:
            filter_func = __filter_func_single_start_relevant
        else:
            filter_func = __filter_func_single_both_relevant
 
    return filter_func

def _demultiplex_hamming(barcode_df, int max_errors):
    """creates a list of quality filters for the PE lane"""
    ret = {}
    forward_rows =[]
    reverse_rows = []
    cdef int len_fwd
    cdef int len_rev
    for row in barcode_df.iterrows():
        key = row[1]['key']
        fwd = row[1]['forward']
        rev = row[1]['reverse']
        fwd_cmpl = fwd[::-1].translate(rev_comp_table)
        rev_cmpl = rev[::-1].translate(rev_comp_table)
        fwd = fwd.upper()
        rev = rev.upper()
        fwd_cmpl = fwd_cmpl.upper()
        rev_cmpl = rev_cmpl.upper()
        len_fwd = len(fwd)
        len_rev = len(rev)
        def qf(quality, sequence,
               fwd = fwd, rev = rev, fwd_cmpl = fwd_cmpl,
               rev_cmpl = rev_cmpl, len_rev = len_rev, len_fwd = len_fwd):
            cdef int dist_start
            cdef int dist_end
            cdef int len_seq

            len_seq = len(sequence)
            #fwd + rev_complement
            dist_start = hamming_ignore_N(sequence[:len_fwd], fwd, max_errors)
            if dist_start <= max_errors:
                dist_end  = hamming_ignore_N(sequence[-len_rev:], rev_cmpl, max_errors)
                if dist_end <= max_errors:
                    return (len_fwd, len_seq - len_rev)
                else:
                    return (len_fwd, len_seq)
            #rev + fwd complement
            else:
                dist_start = hamming_ignore_N(sequence[:len_rev], rev, max_errors)
                if dist_start <= max_errors:
                    dist_end = hamming_ignore_N(sequence[-len_fwd:], fwd_cmpl, max_errors)
                    if dist_end <= max_errors:
                        return (len_rev, len_seq - len_fwd)
                    else:
                        return (len_rev, len_seq)
                else:
                    return (-1, -1)
        ret[key] = qf
    return ret


def read_fastq_iterator_retrieve_index(file_object):
    """retrieve the index-read from name, append to the read and keep track of the indices in seq and quality
    Yield (seq, name, quality)
    """
    row1 = file_object.readline() # name
    row2 = file_object.readline() # seq
    row3 = file_object.readline() # ignored
    row4 = file_object.readline() # quality
    while row1:
        index = row1.split(':')[-1][:-1]
        seq = index+row2[:-1]
        quality = (row4[:-1])
        name = row1[1:-1]
        yield (seq, name, quality)
        row1 = file_object.readline()
        row2 = file_object.readline()
        row3 = file_object.readline()
        row4 = file_object.readline()

def get_fastq_reads(filename, total_index_length = False):
    """get fastq reads - either straight as seq,name, quality or with the index read stiched back onto the seq,
    if @total_index_length is set"""
#    f = chipseq._common.BlockedFileAdaptor(filename)
#    if total_index_length:
#        return read_fastq_iterator_retrieve_index(f)
#    else:
#        return chipseq._common.read_fastq_iterator(f)
    return None

def get_decision_callback(adapter_sequence_begin, adapter_sequence_end, paired, forward_only, max_error_rate):
    adapter_sequence_begin_reverse = adapter_sequence_begin[::-1].translate(rev_comp_table)
    adapter_sequence_end_reverse = adapter_sequence_end[::-1].translate(rev_comp_table)

    #we dont need this
    def __filter_func_start_stop_both_relevant_and_present(tup):
        seq1, qual1, name1 = tup[0]
        seq2, qual2, name2 = tup[1]
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_begin, max_error_rate, True)        
        if read1_i1 == -1 and read2_i1 == -1:
        # forward adapter nowhere to be found, discard
            return False, None, None
        elif read1_i1 == -1 and read2_i1 != -1:
            #adapter in read 2, switch
            seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            read1_i1 = read2_i1
        elif read1_i1 != -1 and read2_i1 != -1:
            raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True)
        if read2_i1 ==-1:
            #stop primer not found, discard
            return False, None, None
        read1_i2 = cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False)
        if read1_i2 == -1:
            return False, None, None
        read2_i2 = cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False)
        if read2_i2 == -1:
            return False, None, None
        return True, ((read1_i1,  read1_i2), (read1_i1,  read1_i2))

    #forward reverse and paired
    def __filter_func_paired_both_relevant(tup):
        seq1, qual1, name1 = tup[0]
        seq2, qual2, name2 = tup[1]
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_begin, max_error_rate, True)        
        if read1_i1 == -1 and read2_i1 == -1:
        # forward adapter nowhere to be found, discard
            return False, None, None
        elif read1_i1 == -1 and read2_i1 != -1:
            #adapter in read 2, switch
            seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            read1_i1 = read2_i1
        elif read1_i1 != -1 and read2_i1 != -1:
            raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True)
        if read2_i1 ==-1:
            #stop primer not found, discard
            return False, None, None
        #accept, both primer found
        read1_i2 = max(0, cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False))
        read2_i2 = max(0, cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False))
        return True, ((read1_i1,  read1_i2), (read1_i1,  read1_i2))

    #forward only and paired
    def __filter_func_paired_start_relevant(tup):            
        seq1, qual1, name1 = tup[0]
        seq2, qual2, name2 = tup[1]
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_begin, max_error_rate, True)        
        if read1_i1 == -1 and read2_i1 == -1:
        # forward adapter nowhere to be found, discard
            return False, None, None
        elif read1_i1 == -1 and read2_i1 != -1:
            #adapter in read 2, switch
            seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            read1_i1 = read2_i1
        elif read1_i1 != -1 and read2_i1 != -1:
            raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        #accept already
        read2_i1 = max(0, cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True))
        read1_i2 = max(0, cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False))
        read2_i2 = max(0, cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False))
        return True, ((read1_i1,  read1_i2), (read1_i1,  read1_i2))

    #forward reverse and single
    def __filter_func_single_both_relevant(read):
        seq1, qual1, name1 = read
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        if read1_i1 == -1:
            return False, None
        read1_i2 = max(0, cutadapt_locate(seq1, adapter_sequence_end, max_error_rate, False))
        return True, ((read1_i1,  read1_i2))

    #forward only and single
    def __filter_func_single_start_relevant(read):
        seq1, qual1, name1 = read
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        if read1_i1 == -1:
            return False, None
        return True, ((read1_i1, len(seq1)))

    if paired:
        if forward_only:
            filter_func = __filter_func_paired_start_relevant
        else:
            filter_func = __filter_func_paired_both_relevant
    else:
        if forward_only:
            filter_func = __filter_func_single_start_relevant
        else:
            filter_func = __filter_func_single_both_relevant
 
    return filter_func


def with_primers(adapter_sequence_begin, adapter_sequence_end):
    adapter_sequence_begin_reverse = adapter_sequence_begin[::-1].translate(rev_comp_table)
    adapter_sequence_end_reverse = adapter_sequence_end[::-1].translate(rev_comp_table)
    def filter(func):
        @functools.wraps(func)
        def wrapper(*args):
            value = func(*args)
            return value
        return wrapper
    return filter

class Demultiplexer:

    def __init__(self, name, barcode_df_or_file, input_strategy, outputdir = None, max_error_rate = 0):
        self.name = name
        if outputdir is None:
            self.result_dir = Path('cache') / self.name
        else:
            self.result_dir = outputdir
            if isinstance(self.result_dir, str):
                self.result_dir = Path(self.result_dir)
        self.result_dir.mkdir(parents = True, exist_ok = True)
        if isinstance(barcode_df_or_file, pd.DataFrame):
            self.barcodes = barcode_df_or_file
        elif isinstance(barcode_df_or_file, str):
            self.barcodes = pd.read_csv(barcode_df_or_file, sep = '\t')
        elif isinstance(barcode_df_or_file, Path):
            self.barcodes = pd.read_csv(str(barcode_df_or_file.resolve()), sep = '\t')
        else:
            raise ValueError("No flanking barcodes supplied. please enter a barcode dataframe or a file path.")
        self.input_strategy = build_fastq_strategy(input_strategy)
        self.maximal_error_rate = max_error_rate
        self.input_files = self.input_strategy()
        self.decision_callbacks = self.initialize_decision_callbacks() #_demultiplex_cutadapt(barcode_df, max_error_rate)

    def initialize_decision_callbacks(self):
        callbacks = {}
        paired = len(self.input_files[0]) == 2
        codes = set()
        forward_only = False
        for code in self.barcodes['reverse'].values:
            if not isinstance(code, str) or code == '':
                forward_only = True
            codes.add(code)
        if len(codes) < len(self.barcodes):
            forward_only = True
        for row in self.barcodes.iterrows():
            key = row[1]['key']
            fwd = row[1]['forward'].encode()
            rev = row[1]['reverse'].encode()
            callbacks[key] = get_decision_callback(fwd, rev, self.maximal_error_rate, paired, forward_only)
#            callbacks[key] = with_primers(fwd, rev)(self.get_decision_callback(paired, forward_only, self.maximal_error_rate)) 

        return callbacks

    def __get_iterator(self, paired):
        myiterator =mbf_align.fastq2.get_iterator("fastq")
        def _iterreads_paired(filetuple):
            for readX in zip(
                myiterator(str(filetuple[0]), reverse_reads = False), myiterator(str(filetuple[1]), reverse_reads = False)
                ):
                    yield readX
        def _iterreads_single(filetuple):
            for readX in myiterator(str(filetuple[0]), reverse_reads = False):
                yield readX
        if paired:
            return _iterreads_paired
        else:
            return _iterreads_single

    def __write_readX(self, readX, fhs):
        for i, read in enumerate(readX):
            fhs[i].write(('@%s\n%s\n+\n%s\n' % read))

    def _decide_on_barcode(self, readX):
        if len(readX) == 2:
            readA, readB = readX
            #check if reads belong together
            if not readA[2].split()[0] == readB[2].split()[0]:
                raise ValueError('Read pairs not matching: %s %s' % (readA[2], readB[2]))
        for key in self.decision_callbacks:
            qf = self.decision_callbacks[key]
            accept, trims, reads = qf(readX)
            if accept:
                trimmed_reads = []
                for i, read in enumerate(reads):
                    trim = trims[i]
                    trimmed = (read[2].decode(), read[0][trim[0]:trim[1]].decode(), read[1][trim[0]:trim[1]].decode())
                    trimmed_reads.append(trimmed)
                return key, tuple(trimmed_reads)
        return 'discarded', readX

    def do_demultiplex(self, dependencies = None):
        if dependencies is not None:
            deps = dependencies
        else:
            deps = []
        paired = len(self.input_files[0]) == 2
        output_files = {}
        filenames = []
        seen = set()
        sample_names = []
        for ftup in self.input_files:
            lib_name = ftup[0].name.split("_S")[0]            
            if lib_name in seen:
                continue
            seen.add(lib_name)
            sample_names.extend([f"{lib_name}_{key}" for key in list(self.decision_callbacks.keys())] + [f"{lib_name}_discarded"])
            for fn in ftup:
                deps.append(ppg.FileChecksumInvariant(fn))
        if paired:
            for sample_name in sample_names:
                (self.result_dir / f"{sample_name}").mkdir(exist_ok = True, parents = True)
                output_files[sample_name] = [self.result_dir / f"{sample_name}" / f"{sample_name}_R1_.fastq", self.result_dir / f"{sample_name}" / f"{sample_name}_R2_.fastq"]
                filenames.extend(output_files[sample_name])
        else:
            for sample_name in sample_name:
                output_files[sample_name] = [self.result_dir / f"{sample_name}" / f"{sample_name}_R1_.fastq"]
                filenames.extend(output_files[sample_name])
        def dump():
            #open tmp files to write in 
            tmp_files = {}
            for sample_name in output_files:
                tmp_files[sample_name] = [open(str(f), "w") for f in output_files[sample_name]]
            read_iterator = self.__get_iterator(paired)
            for ftup in self.input_files:
                lib_name = ftup[0].name.split("_S")[0]            

                for readX in read_iterator(ftup):
                    key, trimmedX = self._decide_on_barcode(readX)
                    sample_name = f"{lib_name}_{key}"
                    self.__write_readX(trimmedX, tmp_files[sample_name])
            for key in tmp_files:
                for f in tmp_files[key]:
                    f.close()
            with (self.result_dir / "done.txt").open('w') as done:
                done.write("demultiplexing done")
        
        params = [self.maximal_error_rate]
        deps.append(ppg.ParameterInvariant(self.name + "_do_demultiplex_params", params))
        return ppg.FileGeneratingJob(str(self.result_dir / "done.txt"), dump).depends_on(deps)

    def clean_up(self, folder, paired, dependencies):
        if paired:
            names = [f"{folder.name}_R1_.fastq", f"{folder.name}_R2_.fastq"]
        else:
            names = [f"{folder.name}_R1_.fastq"]
        outfiles = [folder / name for name in names]
        def __dump():
            for file in folder.iterdir():
                if file.name[:-4] == "_tmp":
                    abspath = file.resolve()
                    shutil.move(abspath, abspath[:-4])
        return ppg.MultiFileGeneratingJob(outfiles, __dump).depends_on(dependencies)

    def get_samples(self, add_dependencies = []):
        res = {}
        paired = len(self.input_files[0]) == 2
        pairing = "paired" if paired else "single"
        dependencies = add_dependencies + [self.do_demultiplex(add_dependencies)]
        seen = set()
        for tup in self.input_files:
            library_name = tup[0].name.split("_S")[0]
            if library_name in seen:
                continue
            seen.add(library_name)
            for key in self.decision_callbacks:
                sample_name = f"{library_name}_{key}"
                folder = self.result_dir / sample_name
                job = self.clean_up(folder, paired, dependencies)
                sample = mbf_align.Sample(
                    sample_name,
                    input_strategy = mbf_align.FASTQsFromJob(job),
                    reverse_reads = False,
                    fastq_processor=mbf_align.fastq2.Straight(),
                    pairing=pairing,
                    vid=None
                    )
                res[sample_name] = sample
        return res




''''
    def __init__(self, name, barcode_df_or_file, input_strategy, outputdir = None, max_error_rate = 0):
        self.name = name
        if outputdir is None:
            self.output_dir = Path('cache') / self.name
        if isinstance(barcode_df_or_file, pd.DataFrame):
            self.barcodes = barcode_df_or_file
        elif isinstance(barcode_df_or_file, str):
            self.barcodes = pd.read_csv(barcode_df_or_file, sep = '\t')
        elif isinstance(barcode_df_or_file, Path):
            self.barcodes = pd.read_csv(str(barcode_df_or_file.resolve()), sep = '\t')
        else:
            raise ValueError("No flanking baqrcodes supplied. please enter a barcode dataframe or a file path.")
        self.input_strategy = build_fastq_strategy(input_strategy)


    def discover_input_files(self):
        """
        a) input is a string or path:
            1 - all input files may be in one folder
                --> check R1 R2 and return a file iterator
            2 - folder may contain subfolders with separate lanes
                --> foreach as 1   
        b) input is a list of files:
            --> as 1
        c) input is a list of tuples/lists:
            --> forearch as 1
        """
        if isinstance(self.input_path, str):
            files = os.listdir(self.input_path)
            pairs = []
            for file in files:
                filename = os.path.join(self.input_path, file)
                if '_R1_' in filename:
                    other_filename = filename.replace('_R1_', '_R2_')
                    if os.path.exists(other_filename):
                        pairs.append((filename, other_filename))
                    else:
                        raise ValueError('No file with mate pairs detected. Shoould have been %s.' % other_filename)
                else:
                    continue
            if len(pairs) == 0:
                raise ValueError('No input files detected for %s.' % self.name)
            pairs_dict = {self.name : pairs}
            return pairs_dict
        elif isinstance(self.input_path, dict):
            return self.input_path
        else:
            raise ValueError("input_path must be either a string or a dictionary of keys to list of paired file tuples")

    def do_demultiplex(self, dependencies = None):
        input_pairs = self.input_strategy()
        print(input_pairs)
        raise ValueError()
        output_keys = []
        for file_key in self.input_files:
            output_keys.extend(["{}_{}".format(file_key, key) for key in list(self.decision_callbacks.keys()) + ['discarded']])
        def dump():
            start_time = time.time()
            outfile_time = open(self.name+'.time', 'w')
            outfile_time.write(("start time = %s\n" % start_time))
            tmp_filenames = []
            output_files = {key: self._open_output_files(key, tmp_filenames) for key in output_keys}
            for file_key in self.input_files:
                for readX in self._iterreads(self.input_files[file_key]):
                    key, trimmedX = self._decide_on_barcode(readX)
                    self.write_readX(trimmedX, output_files, "{}_{}".format(file_key, key))
            for key in output_keys:
                for fh in output_files[key]:
                    fh.close()
            for filename in tmp_filenames:
                shutil.move(filename, filename.replace('.fastq_tmp', '.fastq'))
            stop_time = time.time() - start_time
            outfile_time.write(("elapsed time = %s" % stop_time))
            outfile_time.close()
        if dependencies is not None:
            deps = dependencies
        else:
            deps = []
        params = [self.total_index_length, self.method, self.max_error_rate]
        for key in self.input_files:
            list_of_pairs = self.input_files[key]
            for tuple in list_of_pairs:
                for fn in tuple:
                    deps.append(ppg.FileChecksumInvariant(fn))
                    params.append(fn)
                    os.path.join(self.cache_dir, key, '1_R1_.fastq_tmp')
        filenames = [os.path.join(self.cache_dir, key, '1_R1_.fastq') for key in output_keys]
        filenames += [os.path.join(self.cache_dir, key, '1_R2_.fastq') for key in output_keys]
        deps.append(ppg.ParameterInvariant(self.name + "_do_demultiplex_params", params))
        return ppg.MultiFileGeneratingJob(filenames=filenames,
                                          function=dump, empty_ok = True 
                                          ).depends_on(deps)

    def get_samples(self, add_dependencies = []):
        res = {}
        dependencies = add_dependencies + [self.do_demultiplex(add_dependencies)]
        for input_library_name in self.input_files:
            for key in self.decision_callbacks:
                sample_name = f"{input_library_name}_{key}"
                folder = self.result_dir / f"{input_library_name}_{key}"
                sample = mbf_align.Sample(
                    sample_name,
                    input_strategy = mbf_align.FASTQsFromFolder(folder),
                    reverse_reads = False,
                    fastq_processor=mbf_align.Straight(),
                    pairing="pairing",
                    vid=None
                    )
                res[sample_name] = sample
        return res

class PEMagicDemultiplexer:

    def __init__(self, name, input_path_or_dict_of_key_to_files, barcode_df, method, max_error_rate = 0, total_index_length = False):
        """
        This is a demultiplexer class for those times when you need to do some complicated index juggling.
        @name name of the demultiplexer, this will become a suffix to the lane names
        @input_path path to fastqs
        @method can be 'cut_adapt' or 'hamming'
        They are called in order on each read(pair) and decide whether it belongs to their respective lane.

        @total_index_length if this is a number, it signals that we need to append the index read to the read and
        then take this number of bases as actual index ... this happens when the actual index differs
        from the index length indicated in the sequencing run.
        """
        self.method = method
        self.max_error_rate = max_error_rate
        if self.method == 'cutadapt':
            self.decision_callbacks = _demultiplex_cutadapt(barcode_df, max_error_rate)
        elif self.method == 'hamming':
            self.decision_callbacks = _demultiplex_hamming(barcode_df, max_error_rate)
        else:
            valid_methods = 'cut_adapt', 'hamming'
            raise ValueError("Invalid method %s, allowed: %s" % (method, valid_methods))
        self.input_path = input_path_or_dict_of_key_to_files
        self.name = name
        self.total_index_length = total_index_length
        self.cache_dir = os.path.join('cache', 'PEMagicDemultiplexer', self.name)
        self.input_files = self.discover_input_files()
        ppg.assert_uniqueness_of_object(self)

    def discover_input_files(self):
        if isinstance(self.input_path, str):
            files = os.listdir(self.input_path)
            pairs = []
            for file in files:
                filename = os.path.join(self.input_path, file)
                if '_R1_' in filename:
                    other_filename = filename.replace('_R1_', '_R2_')
                    if os.path.exists(other_filename):
                        pairs.append((filename, other_filename))
                    else:
                        raise ValueError('No file with mate pairs detected. Shoould have been %s.' % other_filename)
                else:
                    continue
            if len(pairs) == 0:
                raise ValueError('No input files detected for %s.' % self.name)
            pairs_dict = {self.name : pairs}
            return pairs_dict
        elif isinstance(self.input_path, dict):
            return self.input_path
        else:
            raise ValueError("input_path must be either a string or a dictionary of keys to list of paired file tuples")
    
    def _open_output_files(self, key, tmp_filenames):
        exptools.common.ensure_path(os.path.join(self.cache_dir, key))
        fname1 = os.path.join(self.cache_dir, key, '1_R1_.fastq_tmp')
        fname2 = os.path.join(self.cache_dir, key, '1_R2_.fastq_tmp')
        tmp_filenames.extend([fname1, fname2])
        return open(fname1, 'w'), open(fname2, 'w')

    def _iterreads(self, pairs):
        for pair in pairs: #self.input_files:
            for read_pair in zip(
                get_fastq_reads(pair[0], self.total_index_length),
                get_fastq_reads(pair[1], self.total_index_length)
                ):
                yield read_pair

    def write_readX(self, readX, output_files, key):
        output_files[key][0].write(('@%s\n%s\n+\n%s\n' % readX[0]))
        output_files[key][1].write(('@%s\n%s\n+\n%s\n' % readX[1]))

    def _decide_on_barcode(self, readX):
        readA, readB = readX
        #check if reads belong together
        #now the qualityfilter function also gets the name of the read (qfs is called with read quality, sequence and name)
        if not readA[1].split()[0] == readB[1].split()[0]:
            raise ValueError('Read pairs not matching: %s %s' % (readA[1], readB[1]))
        #read = (seq, name, quality)
        for key in self.decision_callbacks:
            qf = self.decision_callbacks[key]
            accept, trimA, trimB = qf(readA[2], readA[0], readB[2], readB[0])
#            if b"M03491:55:000000000-C92V5:1:2107:17141:1701" in readA[1]:
#                print(readA[0])
#                print(trimA, accept)
#                print((trimA[0] != -1), (trimA[0] < trimA[1]))
#                print(readB[0])
#                print(trimB, forward_found_B, rev_found_B)
#                print((trimB[0] != -1), (trimB[0] < trimB[1]))
#                print((forward_found_A and rev_found_B) or (forward_found_B and rev_found_A))
#                print(key)
            if accept:
                if (trimA[0] == -1) or (trimA[0] >= trimA[1]):
                    #raise ValueError("Something is wrong with trimA, was {}, read name = {}".format(trimA, readA[1]))
                    continue
                if (trimB[0] == -1) or (trimB[0] >= trimB[1]):
                    #raise ValueError("Something is wrong with trimB, was {}, read name = {}".format(trimB, readB[1]))    
                    continue
                trimmedA = (readA[1].decode(), readA[0][trimA[0]:trimA[1]].decode(), readA[2][trimA[0]:trimA[1]].decode())
                trimmedB = (readB[1].decode(), readB[0][trimB[0]:trimB[1]].decode(), readB[2][trimB[0]:trimB[1]].decode())
                return key, (trimmedA, trimmedB)
        return 'discarded', ((readA[1].decode(), readA[0].decode(), readA[2].decode()), (readB[1].decode(), readB[0].decode(), readB[2].decode()))

    def do_demultiplex(self, dependencies = None):
        input_pairs = self.input_strategy()
        print(input_pairs)
        raise ValueError()
        output_keys = []
        for file_key in self.input_files:
            output_keys.extend(["{}_{}".format(file_key, key) for key in list(self.decision_callbacks.keys()) + ['discarded']])
        def dump():
            start_time = time.time()
            outfile_time = open(self.name+'.time', 'w')
            outfile_time.write(("start time = %s\n" % start_time))
            tmp_filenames = []
            output_files = {key: self._open_output_files(key, tmp_filenames) for key in output_keys}
            for file_key in self.input_files:
                for readX in self._iterreads(self.input_files[file_key]):
                    key, trimmedX = self._decide_on_barcode(readX)
                    self.write_readX(trimmedX, output_files, "{}_{}".format(file_key, key))
            for key in output_keys:
                for fh in output_files[key]:
                    fh.close()
            for filename in tmp_filenames:
                shutil.move(filename, filename.replace('.fastq_tmp', '.fastq'))
            stop_time = time.time() - start_time
            outfile_time.write(("elapsed time = %s" % stop_time))
            outfile_time.close()
        if dependencies is not None:
            deps = dependencies
        else:
            deps = []
        params = [self.total_index_length, self.method, self.max_error_rate]
        for key in self.input_files:
            list_of_pairs = self.input_files[key]
            for tuple in list_of_pairs:
                for fn in tuple:
                    deps.append(ppg.FileChecksumInvariant(fn))
                    params.append(fn)
                    os.path.join(self.cache_dir, key, '1_R1_.fastq_tmp')
        filenames = [os.path.join(self.cache_dir, key, '1_R1_.fastq') for key in output_keys]
        filenames += [os.path.join(self.cache_dir, key, '1_R2_.fastq') for key in output_keys]
        deps.append(ppg.ParameterInvariant(self.name + "_do_demultiplex_params", params))
        return ppg.MultiFileGeneratingJob(filenames=filenames,
                                          function=dump, empty_ok = True 
                                          ).depends_on(deps)

    def get_lanes(self, add_dependencies = []):
        res = {}
        dependencies = add_dependencies + [self.do_demultiplex(add_dependencies)]
        for file_key in self.input_files:
            for key in self.decision_callbacks:
                folder = os.path.join(self.cache_dir, "{}_{}".format(file_key, key))
                name = str('%s_%s_%s' % (self.name, file_key, key))
                files = [(os.path.join(folder, '1_R1_.fastq'), os.path.join(folder, '1_R2_.fastq'))]
                lane = chipseq.lanes.PairedEndLane(name, files, dependencies=dependencies, deps_external=True)
                res['%s_%s_%s' % (self.name, file_key, key)] = lane
        return res


'''
'''
class SEMagicDemultiplexer:

    def __init__(self, name, input_path, decision_callbacks, total_index_length = False):
        """
        This is a demultiplexer class for those times when you need to do some complicated index juggling. This time for Single Reads
	use the PEMagicDemultiplexer for PairedEnd data
        @name name of the demultiplexer, this will become a suffix to the lane names
        @input_path path to fastqs folder
        @decision_callbacks a dictionary of quality filters as for Lanes
        @total_index_length if this is a number, it signals that we need to append the index read to the read and
        then take this number of bases as actual index ... this happens when the actual index differs
        from the index length indicated in the sequencing run.
        """
        self.input_path = input_path
        self.decision_callbacks = decision_callbacks
        self.name = name
        self.total_index_length = total_index_length
        self.cache_dir = os.path.join('cache', 'SEMagicDemultiplexer', self.name)
        self.input_files = self.discover_input_files()
        ppg.assert_uniqueness_of_object(self)

    def discover_input_files(self):
        files = os.listdir(self.input_path)
        inputfiles = []
        for file in files:
            filename = os.path.join(self.input_path, file)
            if filename.endswith('.fastq') or filename.endswith('.fastq.gz'):
                inputfiles.append(filename)
        return inputfiles

    def _open_output_files(self, key, tmp_filenames):
        exptools.common.ensure_path(os.path.join(self.cache_dir, key))
        fname1 = os.path.join(self.cache_dir, key, '1_R1_.fastq_tmp')
        tmp_filenames.extend([fname1])
        return open(fname1, 'wb'),

    def _iterreads(self):
        for fn in self.input_files:
            for read in get_fastq_reads(fn, self.total_index_length):
                yield read

    def write_readX(self, readX, output_files, key):
        output_files[key][0].write('@%s\n%s\n+\n%s\n' % readX)

    def _decide_on_barcode(self, readA):
        #check if reads belong together
        #now the qualityfilter function also gets the name of the read (qfs is called with read quality, sequence and name)
        #read = (seq, name, quality)
        for key in self.decision_callbacks:
            qf = self.decision_callbacks[key]
            trimA = qf(readA[2], readA[0], readA[1])
            if (trimA[0] != -1) and (trimA[0] < trimA[1]):
                trimmedA = (readA[1], readA[0][trimA[0]:trimA[1]], readA[2][trimA[0]:trimA[1]])
                #print '-------------readA-----------'
                #print readA
                #raise ValueError()
                return key, trimmedA
        return 'discarded', (readA[1], readA[0], readA[2])

    def do_demultiplex(self, dependencies = None):
        output_keys = list(self.decision_callbacks.keys()) + ['discarded']
        def dump():
            start_time = time.time()
            outfile_time = open(self.name+'.time', 'wb')
            outfile_time.write("start time = %s\n" % start_time)
            tmp_filenames = []
            output_files = {key: self._open_output_files(key, tmp_filenames) for key in output_keys}
            for readX in self._iterreads():
                key, trimmedX = self._decide_on_barcode(readX)
                self.write_readX(trimmedX, output_files, key)
            for key in output_keys:
                for fh in output_files[key]:
                    fh.close()
            for filename in tmp_filenames:
                shutil.move(filename, filename.replace('.fastq_tmp', '.fastq'))
            stop_time = time.time() - start_time
            outfile_time.write("elapsed time = %s" % stop_time)
            outfile_time.close()
        if dependencies is not None:
            deps = dependencies
        else:
            deps = []
        params = [self.total_index_length, self.input_path, self.name]
        for filename in self.input_files:
            deps.append(ppg.FileChecksumInvariant(filename))
            params.append(filename)
        filenames = [os.path.join(self.cache_dir, key, '1_R1_.fastq') for key in self.decision_callbacks]
        deps.append(ppg.ParameterInvariant(":".join(filenames), params))
        return ppg.MultiFileGeneratingJob(filenames,
                                          dump, 
                                          empty_ok = True).depends_on(deps)
'''