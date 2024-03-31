from abc import ABC, abstractmethod

class subckt(ABC):

    @abstractmethod
    def __init__(self, **nodes):
        pass

    @abstractmethod
    def set_params(self, **params):
        pass

    @abstractmethod
    def add(self, ckt):
        pass

    @abstractmethod
    def add_small(self, op, ckt_sml, **nodes):
        pass
