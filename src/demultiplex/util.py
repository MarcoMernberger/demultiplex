import cutadapt
import cutadapt.align
import collections

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