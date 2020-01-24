from python_code.get_acceptors_and_donors import dna_data_split
import io
from .tests_utils import TEST_DATA, SEQ_FROM_TEST_DATA


def test_dna_data_split():
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_split(f)

    assert sequences == [SEQ_FROM_TEST_DATA]
