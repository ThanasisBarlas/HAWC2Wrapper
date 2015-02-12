
from openmdao.main.api import Assembly

from hawc2_wrapper.hawc2_inputreader import HAWC2InputReader
from hawc2_wrapper.hawc2_inputwriter import HAWC2InputWriter, HAWC2SInputWriter
from openmdao.lib.datatypes.api import List, VarTree, Array, Float, Str, Bool


class ReadWriteH2(Assembly):
    
    def configure(self):

        self.add('reader', HAWC2InputReader())
        self.add('writer', HAWC2InputWriter())
        self.driver.workflow.add(['reader', 'writer'])

        self.reader.htc_master_file = 'hawc2_master_full_K.htc'
        self.writer.from_file = True

        self.connect('reader.vartrees', 'writer.vartrees')

class ReadWriteH2s(Assembly):
    
    
    def configure(self):

        self.add('reader', HAWC2InputReader())
        self.add('writer', HAWC2SInputWriter())
        self.driver.workflow.add(['reader', 'writer'])

        self.reader.htc_master_file = 'hawc2s_master.htc'
        self.writer.from_file = True

        self.connect('reader.vartrees', 'writer.vartrees')

if __name__ == '__main__':

    h2 = ReadWriteH2()
    h2.run()
    # h2s = ReadWriteH2s()
    # h2s.run()