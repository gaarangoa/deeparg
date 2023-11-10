import gzip
import bz2
import zipfile


class Validate():
    def __init__(self):
        self.info = {
            "compressed": False,
            "fileType": None
        }
        self.magic_dict = {
            "\x1f\x8b\x08": "gz",
            "\x42\x5a\x68": "bz2",
            "\x50\x4b\x03\x04": "zip"
        }
        self.max_len = max(len(x) for x in self.magic_dict)

    def compressed(self, filename):
        try:
            f = open(filename)
        except Exception as e:
            f = open(filename+".gz")

        file_start = f.read(self.max_len)
        for magic, filetype in self.magic_dict.items():
            if file_start.startswith(magic):
                return filetype
        return "raw"

    def isFasta(self, line):
        if ">" in line[0]:
            return True
        else:
            return False

    def isFastq(self, line):
        if "@" in line[0]:
            return True
        else:
            return False

    def isGz(self, line):
        if "gz" in line:
            return True
        else:
            return False

    def dataType(self, filename):
        # determine if the file is compressed
        comp = self.compressed(filename)
        # determine if the file is fasta or fastq
        self.info['compressed'] = comp
        if comp == "gz":
            try:
                x = gzip.open(filename)
            except:
                x = gzip.open(filename+".gz")
            x = x.read(10)
        elif comp == "bz2":
            x = bz2.open(filename)
            x = x.read(10)
        elif comp == "zip":
            x = zipfile.open(filename)
            x = x.read(10)
        elif comp == "raw":
            x = open(filename)
            x = x.read(10)
        else:
            self.info['compressed'] = "undefined"

        if self.isFasta(x):
            self.info['fileType'] = "fasta"
            return self.info
        elif self.isFastq(x):
            self.info['fileType'] = "fastq"
            return self.info
        else:
            self.info['fileType'] = "undefined"
            return self.info
