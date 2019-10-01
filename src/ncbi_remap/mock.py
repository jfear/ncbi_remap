from snakemake.io import Namedlist, InputFiles, OutputFiles, Wildcards, Params, Log


class MockSnake:
    """Mock the snakemake class used inside of python scripts."""

    def __init__(self, input=None, output=None, params=None, wildcards=None, log="", config=None):
        self.input = InputFiles(self.make_namedlist(input))
        self.output = OutputFiles(self.make_namedlist(output))
        self.params = Params(self.make_namedlist(params))
        self.wildcards = Wildcards(self.make_namedlist(wildcards))
        self.log = Log(log)
        self.config = config or {}
        self.rulename = "mock"

    def make_namedlist(self, item, first_lvl=True):
        if isinstance(item, list):
            return Namedlist(item)
        elif isinstance(item, dict):
            return Namedlist(fromdict={k: self.make_namedlist(v, False) for k, v in item.items()})
        elif isinstance(item, (str, int, float)):
            if first_lvl:
                return Namedlist([item])
            else:
                return item

        return Namedlist([])
