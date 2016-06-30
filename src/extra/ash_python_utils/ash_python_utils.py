def CreateLogger(name=""):
    import logging
    import sys
    log = logging.getLogger(name)
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log


def AttachFileLogger(log, log_filename, mode="w"):
    import logging
    log_handler = logging.FileHandler(log_filename, mode=mode)
    log.addHandler(log_handler)
    log.info("Log will be written to " + log_filename)


def linear_search(obj, item, start=0):
    for i in xrange(start, len(obj)):
        if obj[i] == item:
            return i
    return -1


def idFormatByFileName(fname):
    import re
    if re.match(r"^.*\.fa(sta)?(\.gz)?$", fname):
        return "fasta"
    elif re.match(r"^.*\.((fq)|(fastq))(\.gz)?$", fname):
        return "fastq"
    else:
        raise "Unrecognized file type"

assert idFormatByFileName("fq.fa") == "fasta"
assert idFormatByFileName("fq.fa.gz") == "fasta"
assert idFormatByFileName("fq.fasta") == "fasta"
assert idFormatByFileName("fq.fasta.gz") == "fasta"
assert idFormatByFileName("fa.fq") == "fastq"
assert idFormatByFileName("fa.fq.gz") == "fastq"
assert idFormatByFileName("fa.fastq") == "fastq"
assert idFormatByFileName("fa.fastq.gz") == "fastq"


def smart_open(filename, mode="r"):
    import gzip
    import re

    if "w" in mode:
        MODE = "w"
    elif "a" in mode:
        MODE = "a"
    else:
        MODE = "r"

    if re.match(r"^.*\.gz$", filename):
        assert(MODE != "a")
        fh = gzip.open(filename, mode=MODE)
    else:
        fh = open(filename, mode=mode)
    return fh


def md5_file(fname):
    import hashlib

    hash_md5 = hashlib.md5()
    with smart_open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def fastx2fastx(input_file, output_file, quality=50):
    from Bio import SeqIO

    input_format, output_format = idFormatByFileName(input_file), idFormatByFileName(output_file)

    with smart_open(input_file) as fin, smart_open(output_file, "w") as fout:
        if input_format == output_format:
            for line in fin:
                fout.write(line)
        else:
            for record in SeqIO.parse(fin, input_format):
                if output_format == "fastq":
                    record.letter_annotations["phred_quality"] = [quality] * len(record)
                SeqIO.write(record, fout, output_format)


fq2fa = fastx2fastx
fa2fq = fastx2fastx


def mkdir_p(path):
    "From http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python"
    import errno
    import os

    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


class FakeLog:

    def info(self, msg):
        print msg

    def warn(self, msg):
        print msg
