import configparser


class BaseReader:
    """Base class for config readers"""
    def __init__(self, config: configparser.ConfigParser, rin):
        self.config = config
        self.rin = rin