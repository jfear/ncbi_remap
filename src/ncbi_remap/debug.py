import os
from tempfile import NamedTemporaryFile

from larval_gonad.config import read_config
from larval_gonad.mock import MockSnake


def snakemake_debug(**kwargs):
    """Helper to debug scripts run by snakemake.

    This is a helper function for debuging scripts run by snakemake. It will
    set the working direcotory if `workdir` is provided. It will return a
    MockSnake object that can be used in place of a snakemake object.

    Possible kwargs include:
        workdir: str
            If provided will set the working directory to this location.
        input: str, list, dict
            Mock inputs for snakemake rule.
        output: str, list, dict
            Mock outputs for snakemake rule.
        params: str, list, dict
            Mock params for snakemake rule.
        wildcards: str, list, dict
            Mock wildcards for snakemake rule.
        log: str
            Mock log for snakemake rule.
        config: str
            Mock config for snakemake rule.
    
    Returns
    -------
    MockSnake: Mock Snakemake object.

    """
    if "workdir" in kwargs:
        try:
            os.chdir(os.path.join(os.getcwd(), kwargs["workdir"]))
            print(os.getcwd())
        except:
            pass

    if "config" in kwargs:
        config = read_config(kwargs["config"])
    else:
        config = None

    if "output" in kwargs:
        output = kwargs["output"]
    else:
        tmp = NamedTemporaryFile()
        output = tmp.name
        print(output)

    return MockSnake(
        input=kwargs.get("input", None),
        output=output,
        params=kwargs.get("params", None),
        wildcards=kwargs.get("wildcards", None),
        log=kwargs.get("log", ""),
        config=config,
    )
