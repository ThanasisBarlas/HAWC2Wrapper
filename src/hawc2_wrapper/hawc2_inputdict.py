"""pure python reader for HAWC2 htc files"""

import numpy as np


def _convert_val(val):
    """ convert a string value read from a file into either int, float or string """

    try:
        return int(val)
    except:
        try:
            return float(val)
        except:
            return val


class Entry(object):
    """class for "name: val" entries in input files"""

    def __init__(self, name, val, comment=''):

        self.name = name
        self.val = val
        self.comment = comment

    def _print(self, tab=0):
        """pretty print the entry"""

        if isinstance(self.val, str):
            return '%s%s' % (tab * ' ', self.name + ' ' + self.val + ' ;\n')
        try:
            return '%s%s' % (tab * ' ', self.name + ' '+ ' '.join(map(str, self.val)) + ' ;\n')
        except:
            return '%s%s' % (tab * ' ', self.name + ' ' + str(self.val) + ' ;\n')

class Section(object):
    """Class for list of Entry objects in a HAWC2 htc file"""

    def __init__(self, name):

        self.name = name
        self.entries = []
        self.comment = ''

    def get_entry(self, name):
        """convenience method to get a specific entry from section"""

        items = []
        for entry in self.entries:
            if entry.name == name:
                if isinstance(entry, Section):
                    return entry
                else:
                    items.append(entry.val)
        if len(items) == 0:
            return None
        if len(items) == 1:
            return items[0]
        else:
            return items

    def _next(self):

        try:
            entry = self.entries[self.pos]
            self.pos += 1
            return entry
        except:
            return False

    def _print(self, tab=0):
        """pretty print section recursively"""

        self.pos = 0
        string = ''

        string += '%sbegin %s' % (tab * ' ', self.name + ' ;\n')
        while self.pos < len(self.entries):
            entry = self._next()
            string += entry._print(tab=tab + 2)
        string += '%send %s' % (tab * ' ', self.name + ' ;\n')

        return string


class HAWC2InputDict(object):
    """
    class for reading a HAWC2 htc file

    file is read into a nested list with Section objects with lists of
    Entry objects with name, val for each input parameter.

    all values are converted to either integer, float or string

    pc, ae, and st files are not read.
    """

    def __init__(self):

        self.body_order = []

    def _print(self):

        string = ''
        for section in self.htc:
            string += section._print()
        return string

    def read(self, filename):

        fid = open(filename, 'r')
        self.lines = fid.readlines()
        self.nl = len(self.lines)
        self.htc = []
        self.pos = 0

        tmp = self._read_section(Section('tmp'))
        self.htc = tmp.entries
        self.inp = tmp

    def _read_section(self, section):
    
        while self.pos < self.nl:
            line = self._next()
            comment = line[line.find(';')+1:].strip()
            param = line[:line.find(';')]
            param = param.strip().split()
            if len(param) == 0:
                continue
            elif param[0] == ';':
                continue
            elif param[0] == 'begin':
                newsection = Section(param[1])
                sec = self._read_section(newsection)
                section.entries.append(sec)
            elif param[0] == 'end':
                return section
            elif param[0] == 'exit':
                return section
            elif param[0] == 'continue_in_file':
                htc = HAWC2InputDict()
                htc.read(param[1])
                self.body_order.extend(htc.body_order)
                section.entries.extend(htc.htc)
            else:
                vals = [_convert_val(val) for val in param[1:]]
                if param[0] == 'name' and section.name == 'main_body':
                    self.body_order.append(param[1])

                if len(vals) == 1:
                    section.entries.append(Entry(param[0],vals[0], comment))
                else:
                    section.entries.append(Entry(param[0], vals, comment))

        return section

    def _next(self):

        try:
            line = self.lines[self.pos]
            self.pos += 1
            return line 
        except:
            return False


def read_hawc2_st_file(filename):
    """
    Reader for a HAWC2 beam structure file.

    The format of the file should be:
    main_s[0] dm[1] x_cg[2] y_cg[3] ri_x[4] ri_y[5] x_sh[6] y_sh[7] E[8] ...
    G[9] I_x[10] I_y[11] K[12] k_x[13] k_y[14] A[15] pitch[16] x_e[17] y_e[18]

    Sub-classes can overwrite this function to change the reader's behaviour.
    """
    fid = open(filename, 'r')
    st_sets = [] 
    line = fid.readline()
    while line:
        line = fid.readline()
        if line.find('$') != -1:
            ni = int(line.split()[1])
            st_data = np.zeros((ni, 19)) 
            for i in range(ni):
                tmp = fid.readline().split()
                st_data[i, :] = [float(tmp[j]) for j in range(19)]

            st = {}
            st['s'] = st_data[:, 0]
            st['dm'] = st_data[:, 1]
            st['x_cg'] = st_data[:, 2]
            st['y_cg'] = st_data[:, 3]
            st['ri_x'] = st_data[:, 4]
            st['ri_y'] = st_data[:, 5]
            st['x_sh'] = st_data[:, 6]
            st['y_sh'] = st_data[:, 7]
            st['E'] = st_data[:, 8]
            st['G'] = st_data[:, 9]
            st['I_x'] = st_data[:, 10]
            st['I_y'] = st_data[:, 11]
            st['K'] = st_data[:, 12]
            st['k_x'] = st_data[:, 13]
            st['k_y'] = st_data[:, 14]
            st['A'] = st_data[:, 15]
            st['pitch'] = st_data[:, 16]
            st['x_e'] = st_data[:, 17]
            st['y_e'] = st_data[:, 18]
            st_sets.append(st)
    fid.close()
    return st_sets

def read_hawc2_pc_file(filename):
    """Read a pc airfoil data file into a dictionary"""

    fid = open(filename, 'r')
    ltmp = fid.readline().strip('\r\n').split()
    nset = int(ltmp[0])
    desc = " ".join(ltmp[1:])
    sets = []
    for i in range(nset):
        pcset = {}
        pcset['polars'] = []
        npo =  int(fid.readline().split()[0])
        rawdata = [row.strip().split('\t') for row in fid]
        rthick = []
        ii = 0
        for n in range(npo):
            polar = {}
            line = rawdata[ii][0].split()
            ni = int(line[1])
            polar['np'] = ni
            polar['rthick'] = float(line[2])
            polar['desc'] = " ".join(line[3:])
            data = np.zeros((ni,4))
            ii += 1
            for i in range(ni):
                dline = rawdata[ii][0]
                data[i, :] = [float(var) for var in dline.split()]
                ii += 1
            polar['aoa'] = data[:, 0]
            polar['cl'] = data[:, 1]
            polar['cd'] = data[:, 2]
            polar['cm'] = data[:, 3]
            rthick.append(polar['rthick'])
            pcset['polars'].append(polar)
        pcset['rthick'] = rthick
        sets.append(pcset)

    fid.close()
    return [desc, sets]

def read_hawc2_ae_file(filename):
    """read blade chord and relative thickness from an ae file"""

    fid = open(filename, 'r')
    blade_ae = {}

    # read header
    line = fid.readline()  # '1 HAWC_AE data\n'
    line = fid.readline()

    # read data
    data = np.loadtxt(fid, usecols=(0, 1, 2, 3))
    blade_ae['s'] = data[:, 0]
    blade_ae['chord'] = data[:, 1]
    blade_ae['rthick'] = data[:, 2]
    blade_ae['aeset'] = data[:, 3]

    fid.close()

    return blade_ae

if __name__ == '__main__':

    r = HAWC2InputReader()
    r.read('hawc2_master.htc')

