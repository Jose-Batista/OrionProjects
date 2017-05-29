
import pickle
from floe.api import InputPort, OutputPort


class PickleEncoder(object):
    def encode(self, item):
        return pickle.dumps(item)

    def decode(self, data):
        return pickle.loads(data)

class ObjectInputPort(PickleEncoder, InputPort):
    pass

class ObjectOutputPort(PickleEncoder, OutputPort):
    pass

class VirtualScreeningData:
    def __init__(self):
        self.act_list = list()
        self.baitset = list()

class ActiveList:
    def __init__(self):
        self.act_list = list()
        self.fp_list = list()
        self.baitset = list()

