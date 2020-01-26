from python_code.get_acceptors_and_donors import dna_data_read, get_acceptors_and_donors, DONOR_SEQ, ACCEPTOR_SEQ
import io
from .tests_utils import TEST_DATA, SEQ_FROM_TEST_DATA


def test_dna_data_split():
    """
    Test of file is read in proper way
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    assert sequences == [SEQ_FROM_TEST_DATA]


def test_get_acceptors_and_donors_1():
    """
    Dummy test which check if read sequences of length 2 ar only donor/acceptor sequence.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    t_d, f_d, t_a, f_a = get_acceptors_and_donors(0, 2, sequences)

    for x in t_d + f_d:
        assert x == DONOR_SEQ

    for x in t_a + f_a:
        assert x == ACCEPTOR_SEQ


def test_get_acceptors_and_donors_2():
    """
    Check if donor/acceptors sequence is in correct position.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    a = 2
    b = 3
    t_d, f_d, t_a, f_a = get_acceptors_and_donors(a, b, sequences)

    print(t_d, f_d, t_a, f_a)
    for x in t_d + f_d:
        assert x[a: a+2] == DONOR_SEQ

    for x in t_a + f_a:
        assert x[a: a+2] == ACCEPTOR_SEQ
