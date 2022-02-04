#embed_pos_in_docstring
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
import cutadapt.align
from collections import namedtuple

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans


cdef:
    struct Read:
        char* seq
        char* name
        char* qual


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
        if start:
            return pos + len(adapter)  # we want to trim the adapter
        else:
            return pos  # if its the end adapter, we need to trim at the beginning
    #if not, use cutadapt to find the best match
    elif max_error_rate == 0:
        return -1
    else:
        locate_adapter = locate(adapter, max_error_rate, start)        
        match = locate_adapter(seq)
        if match is None:
            return -1 
        if start:
            return match.rstop  # we want to trim the adapter
        else:
            return match.rstart  # if its the end adapter, we need to trim at the beginning


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
            if read2_i1 < read1_i1:
                #adapter in read 2, something freakish with read12
                seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
                read1_i1 = read2_i1
            read2_i1 = -1
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
        if (read1_i1>=read1_i2) or (read2_i1>=read2_i2):
            return False, None, None
        return True, ((read1_i1,  read1_i2), (read2_i1,  read2_i2))

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
            read1_i1, read2_i1 = read2_i1, read1_i1
        elif read1_i1 != -1 and read2_i1 != -1:
            if read2_i1 < read1_i1:
                #adapter in read 2, something freakish with read12
                seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
                read1_i1, read2_i1 = read2_i1, read1_i1
            #raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        read2_i1 = cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True)
        if read2_i1 ==-1:
            #stop primer not found, discard
            return False, None, None
        #accept, both primer found
        read1_i2 = cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False)
        if read1_i2 < 1:
            read1_i2 = len(seq1)
        read2_i2 = cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False)
        if read2_i2 < 1:
            read2_i2 = len(seq2)
        if (read1_i1>=read1_i2) or (read2_i1>=read2_i2):
            return False, None, None
        return True, ((read1_i1,  read1_i2), (read1_i1,  read1_i2)), ((seq1, qual1, name1), (seq2, qual2, name2))

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
            if read2_i1 < read1_i1:
                #adapter in read 2, something freakish with read12
                seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
                read1_i1 = read2_i1
            read2_i1 = -1
            #raise ValueError("adapter_Sequence_begin found in both reads, this should not happen.")
        else:
            pass
        #accept already
        #read1_i1 = start index for slice on read 1
        #read1_i2 = stop index for slice on read 1
        #read2_i1 = start index for slice on read 2
        #read2_i2 = stop index for slice on read 2
        read1_i2 = cutadapt_locate(seq1, adapter_sequence_end_reverse, max_error_rate, False)
        if read1_i2 < 1: 
            read1_i2 = len(seq1)
        read2_i1 = max(0, cutadapt_locate(seq2, adapter_sequence_end, max_error_rate, True))
        read2_i2 = cutadapt_locate(seq2, adapter_sequence_begin_reverse, max_error_rate, False)
        if read2_i2 < 1:
            read2_i2 = len(seq2)
        if (read1_i1>=read1_i2) or (read2_i1>=read2_i2):
            return False, None, None
        return True, ((read1_i1,  read1_i2), (read1_i1,  read1_i2)), ((seq1, qual1, name1), (seq2, qual2, name2))

    #forward reverse and single
    def __filter_func_single_both_relevant(read):
        seq1, qual1, name1 = read
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        if read1_i1 == -1:
            return False, None, None
        read1_i2 = cutadapt_locate(seq1, adapter_sequence_end, max_error_rate, False)
        if read1_i2 == -1:
            return False, None, None
        if (read1_i1>=read1_i2):
            return False, None, None
        return True, ((read1_i1,  read1_i2)), ((seq1, qual1, name1))

    #forward only and single
    def __filter_func_single_start_relevant(read):
        seq1, qual1, name1 = read
        read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
        if read1_i1 == -1:
            #if not found, try it reverse
            seq1 = seq1[::-1].translate(rev_comp_table)
            qual1 = qual1[::-1]
            read1_i1 = cutadapt_locate(seq1, adapter_sequence_begin, max_error_rate, True)
            if read1_i1 == -1:
                #if still not found: discard
                return False, None, None
        return True, ((read1_i1, len(seq1))), ((seq1, qual1, name1))
    
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

    def __init__(self, sample, barcode_df_or_file, outputdir = None, max_error_rate = 0, filter_func = None):
        self.name = f"DM_{sample.name}"
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
        self.input_sample = sample
        self.maximal_error_rate = max_error_rate
        self.input_files = [self.input_sample.get_aligner_input_filenames()]
        self.pairing = len(self.input_files[0]) == 2
        self.decision_callbacks = self.initialize_decision_callbacks()
        self.filter_func = filter_func

    def initialize_decision_callbacks(self):
        callbacks = {}
        codes = set()
        forward_only = False
        for code in self.barcodes['reverse'].values:
            if not isinstance(code, str) or code == '':
                forward_only = True
            codes.add(code)
        if len(codes) < len(self.barcodes):
            forward_only = True
        for _, row in self.barcodes.iterrows():
            key = row['key']
            fwd = row['forward'].encode()
            rev = row['reverse'].encode()
            callbacks[key] = get_decision_callback(fwd, rev, self.pairing, forward_only, self.maximal_error_rate)

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

    def __decode(self, read):
        return (read[0].decode(), read[1].decode(), read[2].decode())

    def _decide_on_barcode(self, readX):
        if self.pairing:
            readA, readB = readX
            #check if reads belong together
            if not readA[2].split()[0] == readB[2].split()[0]:
                raise ValueError('Read pairs not matching: %s %s' % (readA[2], readB[2]))
        else:
            readX = self.__decode(readX)
        for key in self.decision_callbacks:
            qf = self.decision_callbacks[key]
            accept, trims, reads = qf(readX)
            if accept:
                trimmed_reads = []
                for i, read in enumerate(reads):
                    trim = trims[i]
                    trimmed = self.__decode((read[2], read[0][trim[0]:trim[1]], read[1][trim[0]:trim[1]]))
                    trimmed_reads.append(trimmed)
                return key, tuple(trimmed_reads)
        if self.pairing:
            readX = self.__decode((readX[0][2], readX[0][0], readX[0][1])), self.__decode((readX[1][2], readX[1][0], readX[1][1]))
        else:
            readX = self.__decode(readX[2], readX[0], readX[1])
        return 'discarded', readX

    def do_demultiplex(self, dependencies = []):
        deps = dependencies
        output_files = {}
        filenames = []
        lib_name = self.input_sample.name            
        sample_names = [f"{lib_name}_{key}" for key in list(self.decision_callbacks.keys())] + [f"{lib_name}_discarded"]
        if self.pairing:
            for sample_name in sample_names:
                (self.result_dir / sample_name).mkdir(parents = True, exist_ok = True)
                output_files[sample_name] = [self.result_dir / sample_name / f"{sample_name}_R1_.fastq", self.result_dir / sample_name / f"{sample_name}_R2_.fastq"]
                filenames.extend(output_files[sample_name])
        else:
            for sample_name in sample_name:
                (self.result_dir / sample_name).mkdir(parents = True, exist_ok = True)
                output_files[sample_name] = [self.result_dir / sample_name / f"{sample_name}_R1_.fastq"]
                filenames.extend(output_files[sample_name])
        def dump():
            #open tmp files to write in 
            tmp_files = {}
            for sample_name in output_files:
                tmp_files[sample_name] = [Path(str(f)+".tmp").open("w") for f in output_files[sample_name]]
            read_iterator = self.__get_iterator(self.pairing)
            lib_name = self.input_sample.name            
            for ftup in self.input_files:
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
        deps.append(self.input_sample.prepare_input())
        return ppg.FileGeneratingJob(str(self.result_dir / "done.txt"), dump).depends_on(deps)

    def clean_up(self, lib_name, sample_name, dependencies):
        if self.pairing:
            outfiles = [
                self.result_dir / sample_name / f"{sample_name}_R1_.fastq",
                self.result_dir / sample_name / f"{sample_name}_R2_.fastq"
            ]
        else:
            outfiles = [
                self.result_dir / sample_name / f"{sample_name}_R1_.fastq"
            ]
        folder = self.result_dir / sample_name
        def __dump():
            for file in folder.iterdir():
                if file.suffix == ".tmp":
                    shutil.move(str(file), str(file.parent / file.stem))

        return ppg.MultiFileGeneratingJob(outfiles, __dump).depends_on(dependencies)

    def get_samples(self, add_dependencies = []):
        if hasattr(self, "samples"):
            return self.samples
        else:
            self.register_samples(add_dependencies)
            return self.samples
        

    def register_samples(self, add_dependencies = []):
        res = {}
        pairing = "paired" if self.pairing else "single"
        dependencies = add_dependencies + [self.do_demultiplex(add_dependencies)]
        library_name = self.input_sample.name
        fastq_processor = mbf_align.fastq2.Straight()
        if self.filter_func is not None:
            fastq_processor = mbf_align.fastq2.Paired_Filtered_Trimmed(self.filter_func)
        for key in list(self.decision_callbacks.keys()) + ["discarded"]:
            sample_name = f"{library_name}_{key}"
            folder = self.result_dir / sample_name
            job = self.clean_up(library_name, sample_name, dependencies)
            sample = mbf_align.Sample(
                sample_name,
                input_strategy = mbf_align.FASTQsFromJob(job),
                reverse_reads = False,
                fastq_processor = fastq_processor,
                pairing = pairing,
                vid = None
                )
            res[sample_name] = sample
        self.samples = res

cdef int axpy(int a, int x, int y):
    return a * x + y