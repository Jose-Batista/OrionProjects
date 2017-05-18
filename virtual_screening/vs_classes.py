
import pickle
from floe.api import InputPort, OutputPort


class PickleEncoder(object):
    def encode(self, item):
        return pickle.dumps(item)

    def decode(self, data):
        return pickle.loads(data)

class MyInputPort(PickleEncoder, InputPort):
    pass

class MyOutputPort(PickleEncoder, OutputPort):
    pass

class VirtualScreeningData:
    def __init__(self):
        self.act_list = list()
        self.baitset = list()
