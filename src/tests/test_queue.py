from pathlib import Path
import pickle
import pytest

from ncbi_remap.queue import Queue

TEST_DATA = Path(Path(__file__).parent / "../test_data/srx2srr.csv").absolute()

@pytest.fixture(scope="session")
def srx_folder(tmpdir_factory):
    tmp = tmpdir_factory.mktemp("SRXs")
    tmp.ensure("SRXs/DRX000774")
    tmp.ensure("SRXs/DRX000775")
    tmp.ensure("SRXs/DRX000998")
    return str(tmp / "SRXs")


@pytest.fixture(scope="session")
def srx_list():
    return ["DRX000774", "DRX000775", "DRX000998"]


@pytest.fixture(scope="session")
def srx_pkl(tmpdir_factory):
    srxs = set(["DRX000774", "DRX000775", "DRX000998"])
    fn = tmpdir_factory.mktemp("data").join("srxs.pkl")
    pickle.dump(srxs, open(fn, "wb"))
    return str(fn)


@pytest.fixture(scope="session")
def srr_folder(tmpdir_factory):
    tmp = tmpdir_factory.mktemp("SRRs")
    tmp.ensure("SRRs/DRR001177")
    tmp.ensure("SRRs/DRR001178")
    tmp.ensure("SRRs/DRR001444")
    return str(tmp / "SRRs")


@pytest.fixture(scope="session")
def srr_list():
    return ["DRR001177", "DRR001178", "DRR001444"]


@pytest.fixture(scope="session")
def srr_pkl(tmpdir_factory):
    srrs = set(["DRR001177", "DRR001178", "DRR001444"])
    fn = tmpdir_factory.mktemp("data").join("srrs.pkl")
    pickle.dump(srrs, open(fn, "wb"))
    return str(fn)


def test_base_queue():
    queue = Queue(srx2srr=TEST_DATA)
    assert queue.n_srxs > 0
    assert queue.n_srrs > 0


def test_srx_target_folder(srx_folder):
    queue = Queue(targets=srx_folder, srx2srr=TEST_DATA)
    assert queue.n_srxs == 3
    assert queue.n_srrs > 0


def test_srx_target_list(srx_list):
    queue = Queue(targets=srx_list, srx2srr=TEST_DATA)
    assert queue.n_srxs == 3
    assert queue.n_srrs > 0


def test_srx_target_pkl(srx_pkl):
    queue = Queue(targets=srx_pkl, srx2srr=TEST_DATA)
    assert queue.n_srxs == 3
    assert queue.n_srrs > 0


def test_srr_target_folder(srr_folder):
    queue = Queue(targets=srr_folder, srx2srr=TEST_DATA)
    assert queue.n_srxs > 0
    assert queue.n_srrs == 3


def test_srr_target_list(srr_list):
    queue = Queue(targets=srr_list, srx2srr=TEST_DATA)
    assert queue.n_srxs > 0
    assert queue.n_srrs == 3


def test_srr_target_pkl(srr_pkl):
    queue = Queue(targets=srr_pkl, srx2srr=TEST_DATA)
    assert queue.n_srxs == 3
    assert queue.n_srrs > 0


def test_srx_problems_folder(srx_folder):
    queue = Queue(problems=srx_folder, srx2srr=TEST_DATA)
    assert queue.n_srxs == len(queue.srx2srr.srx.unique()) - 3


def test_srr_problems_folder(srr_folder):
    queue = Queue(problems=srr_folder, srx2srr=TEST_DATA)
    assert queue.n_srrs == len(queue.srx2srr.srr.unique()) - 3


def test_srx_subset_folder(srx_folder):
    queue = Queue(subset=srx_folder, srx2srr=TEST_DATA)
    assert queue.n_srxs == 3


def test_srr_subset_folder(srr_folder):
    queue = Queue(subset=srr_folder, srx2srr=TEST_DATA)
    assert queue.n_srrs == 3


def test_srx_sample_table():
    queue = Queue(targets=["DRX000774", "DRX000775", "SRX1006261"], size=2, srx2srr=TEST_DATA)
    assert sorted(queue.srxs) == sorted(queue.sample_table.srx.unique().tolist())


def test_get_srx():
    queue = Queue(srx2srr=TEST_DATA)
    srx = queue.get_srx("DRR001177")
    assert srx == "DRX000774"


def test_get_srrs():
    queue = Queue(srx2srr=TEST_DATA)
    srrs = queue.get_srrs("DRX000774")
    assert srrs == ["DRR001177"]

def test_srr_expand():
    srx = "DRX000774"
    srr = "DRR001177"
    queue = Queue(srx2srr=TEST_DATA)
    assert queue.expand("test/{srx}/{srr}", srr) == f"test/{srx}/{srr}"

def test_srx_expand():
    srx = "SRX1006261"
    queue = Queue(srx2srr=TEST_DATA)
    assert len(queue.expand("test/{srx}/{srr}", srx)) == 4
